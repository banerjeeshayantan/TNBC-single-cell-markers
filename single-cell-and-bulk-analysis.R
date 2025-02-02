library(SingleCellExperiment)
library(org.Hs.eg.db)
library(tidyverse)
library(SingleR)
library(mlbench)
library(MASS)
library(randomForest)
library(pROC)
library(scater)
library(cowplot)
library(caret)
library(CelliD)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(Seurat)

# Load protein-coding genes to filter out non-protein-coding ones
data("HgProteinCodingGenes")

# Read the expression matrix file
# The dataset contains gene expression values in TPM (Transcripts Per Million)
tnbc <- read.delim("~/TNBC-single-cell/latest-codes/TNBC_expression_matrix_TPM.txt")
dim(tnbc) # Check dimensions: 23368 genes by 5873 cells

# Read the annotation file containing metadata about each cell
tnbc_anno <- read.delim("~/TNBC-single-cell/TNBC-archive/TNBC_anno.txt")

# Create a Seurat object from the raw expression matrix
tnbc_seurat <- CreateSeuratObject(counts = tnbc)

# Apply log transformation to normalize data (log(TPM/10+1))
tnbc_seurat <- SetAssayData(object = tnbc_seurat, slot="data", new.data = log1p(as.matrix(tnbc/10)))

# Remove mid and post-treatment cells, retaining only pre-treatment cells
tnbc_anno <- tnbc_anno[!tnbc_anno$Treatment_status %in% c("Mid", "Post"),]

# Subset the Seurat object to retain only pre-treatment cells
selected_cells <- which(names(tnbc_seurat$orig.ident) %in% tnbc_anno$Sample_id)
tnbc_seurat_pre <- tnbc_seurat[, selected_cells]

# Remove a cell that lacks corresponding phenotyping data
tnbc_seurat_pre <- tnbc_seurat_pre[, !names(tnbc_seurat_pre$orig.ident) == "KTN206_P_CCAAGCCATGG"]

# Remove housekeeping genes from Tirosh et al., which maintain constant expression
housekeeping_genes <- read.delim("~/TNBC-single-cell/TNBC-archive/housekeeping_genes.txt", header = FALSE)
selected_genes <- setdiff(rownames(tnbc_seurat_pre), housekeeping_genes$V1)
tnbc_seurat_pre <- tnbc_seurat_pre[selected_genes, ]

# Keep only protein-coding genes for downstream analysis
tnbc_protein_coding <- tnbc_seurat_pre[rownames(tnbc_seurat_pre) %in% HgProteinCodingGenes,]

dim(tnbc_protein_coding) # Check dimensions: 17406 genes by 2511 cells

# Compute mitochondrial gene expression ratio and remove high mitochondrial cells
tnbc_protein_coding$mitoRatio <- PercentageFeatureSet(object = tnbc_protein_coding, pattern = "^MT") / 100
filtered_seurat <- subset(tnbc_protein_coding, mitoRatio < 0.05) # Retain cells with <5% mitochondrial content

dim(filtered_seurat) # Cells reduced to 2034

# Remove genes expressed in fewer than 10 cells
filtered_seurat <- filtered_seurat[rowSums(filtered_seurat@assays$RNA@counts > 0) >= 10, ]

# Update metadata to match the remaining cells
idx6 <- which(tnbc_anno$Sample_id %in% colnames(filtered_seurat))
tnbc_anno_updated <- tnbc_anno[idx6,]
filtered_seurat@meta.data$Chemo_labels <- tnbc_anno_updated$Labels

# Cell cycle scoring using pre-defined gene lists
load("~/TNBC-single-cell/TNBC-archive/cycle.rda")
seurat_phase <- CellCycleScoring(filtered_seurat, g2m.features = g2m_genes, s.features = s_genes)

# Identify most variable genes for feature selection
seurat_phase <- FindVariableFeatures(seurat_phase, selection.method = "vst", nfeatures = 3000, verbose = FALSE)

# Run PCA to visualize variance explained by different principal components
seurat_phase <- RunPCA(seurat_phase)

# Perform UMAP and t-SNE dimensionality reduction
seurat_phase <- RunUMAP(seurat_phase, dims = 1:40)
seurat_phase <- RunTSNE(seurat_phase, dims = 1:40)

# Integrate the dataset using Seurat's integration method
split_seurat <- SplitObject(seurat_phase, split.by = "Chemo_labels")
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 2000)
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, reduction = "cca", anchor.features = integ_features, dims = 1:40, k.anchor = 20)
seurat_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "LogNormalize", k.weight = 60, features.to.integrate = rownames(seurat_phase))
seurat_integrated <- ScaleData(seurat_integrated)

# Run PCA and visualize clusters after integration
seurat_integrated <- RunPCA(object = seurat_integrated)
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40, reduction = "pca")
seurat_integrated <- RunTSNE(seurat_integrated, dims = 1:40, reduction = "pca")

# Perform clustering to identify distinct cellular populations
seurat_integrated <- FindNeighbors(object = seurat_integrated, dims = 1:40)
seurat_integrated <- FindClusters(object = seurat_integrated, resolution = 0.5)

# Perform differential gene expression (DEG) analysis
seurat_integrated_markers <- FindAllMarkers(seurat_integrated, min.pct = 0.25, logfc.threshold = 1, only.pos = TRUE)
seurat_integrated_markers <- seurat_integrated_markers[seurat_integrated_markers$p_val_adj < 0.05, ]

# Save results for downstream analysis
write.table(seurat_integrated_markers, "~/TNBC-single-cell/latest-codes/Pre_only_TNBC_Markers_FindAllMarkers_Res=0.4_9_clusters.txt", sep="\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# Save integrated Seurat object
saveRDS(seurat_integrated, "~/TNBC-single-cell/latest-codes/Seurat_Integrated_Post_Only_2371_cells.rds")


# Save marker genes identified in FindAllMarkers
write.table(seurat_integrated_markers, "~/TNBC-single-cell/latest-codes/Pre_only_TNBC_Markers_FindAllMarkers_Res=0.4_9_clusters.txt", sep="\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# Subset T and B cell populations based on marker gene expression
DefaultAssay(seurat_integrated) <- "RNA"
seurat_tandb_cell_cd8_cd4_new <- subset(seurat_integrated, subset = (CD8A > 0 | CD8B > 0 | CD4 > 0 | CD19 > 0 | MS4A1 > 0))
DefaultAssay(seurat_tandb_cell_cd8_cd4_new) <- "integrated"

# Identify highly variable genes in the subsetted population
seurat_tandb_cell_cd8_cd4_new <- FindVariableFeatures(seurat_tandb_cell_cd8_cd4_new, 
                                                      selection.method = "vst", 
                                                      nfeatures = 3000, 
                                                      verbose = FALSE)
DefaultAssay(seurat_tandb_cell_cd8_cd4_new) <- "RNA"

# Perform differential expression analysis
seurat_tandb_cell_integrated_markers <- FindAllMarkers(seurat_tandb_cell_cd8_cd4_new, min.pct = 0.25, logfc.threshold = 1, only.pos = TRUE)
seurat_tandb_cell_integrated_markers <- seurat_tandb_cell_integrated_markers[seurat_tandb_cell_integrated_markers$p_val_adj < 0.05,]
seurat_tandb_cell_integrated_markers <- seurat_tandb_cell_integrated_markers[order(seurat_tandb_cell_integrated_markers$avg_log2FC, decreasing = TRUE),]

# Perform Gene Ontology enrichment analysis
all_genes_human <- unique(rownames(seurat_phase))
GO_BP <- enrichGO(gene = seurat_tandb_cell_integrated_markers$gene, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', 
                  ont = "BP", pAdjustMethod = "BH", universe = all_genes_human, pvalueCutoff = 0.01, qvalueCutoff = 0.05)
GO_CC <- enrichGO(gene = seurat_tandb_cell_integrated_markers$gene, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', 
                  ont = "CC", pAdjustMethod = "BH", universe = all_genes_human, pvalueCutoff = 0.01, qvalueCutoff = 0.05)
GO_MF <- enrichGO(gene = seurat_tandb_cell_integrated_markers$gene, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', 
                  ont = "MF", pAdjustMethod = "BH", universe = all_genes_human, pvalueCutoff = 0.01, qvalueCutoff = 0.05)

# Generate dot plots for GO terms
dotplot(GO_BP, x="count", showCategory=20, color = 'p.adjust', font.size = 8) +
  dotplot(GO_CC, x="count", showCategory=20, color = 'p.adjust', font.size = 8) +
  dotplot(GO_MF, x="count", showCategory=20, color = 'p.adjust', font.size = 8)

# Perform PCA visualization
PCAPlot(seurat_tandb_cell_cd8_cd4_new, split.by = "Chemo_labels")  

# Save filtered list of marker genes
write.table(intersect(seurat_tandb_cell_integrated_markers$gene, pooled_sig_gene_list), 
            "~/TNBC-single-cell/TNBC-archive/Top_features/latest_feats_5April_DEGs_extracted_from_TILS.txt", 
            col.names = TRUE, row.names = FALSE, sep="\t", quote = FALSE)

# Save the Seurat object
saveRDS(seurat_t_cell_cd8_cd4_new, "~/TNBC-single-cell/latest-codes/Seurat_T_cells_CD8_CD4.rds") 

#just subsetting the expression values for B cells
mat <- GetAssayData(object = seurat_phase_t_cell, assay = "RNA", slot = "data")
seurat_phase_t_cell$Chemo_labels
mat = as.matrix(mat)
mat[1:10,1:10]
mat_t = as.data.frame(t(mat))
mat_t$Labels = seurat_phase_t_cell$Chemo_labels

d = mat_t[, -c(which(colnames(mat_t) %in% "Labels"))]
yvar = mat_t$Labels
y = factor(yvar)
levels(y) = c("0","1")
getwd()
write.table(mat_t, "TNBC_expression_Bcells_only.txt", col.names = T, row.names = T, sep="\t", quote = F)


# Perform cell type annotation using SingleR
library(celldex)
ref.data <- BlueprintEncodeData()
DefaultAssay(seurat_integrated) <- "integrated"
seurat_integrated_sce <- as.SingleCellExperiment(seurat_integrated)
pred_be <- SingleR(test = seurat_integrated_sce, ref = ref.data, assay.type.test=1, labels = ref.data$label.main)

# Perform classification using glmnet
require(glmnet)
mat <- GetAssayData(object = seurat_phase_t_cell, assay = "RNA", slot = "data")
mat_t <- as.data.frame(t(as.matrix(mat)))
mat_t$Labels <- seurat_phase_t_cell$Chemo_labels
d <- mat_t[, -which(colnames(mat_t) %in% "Labels")]
y <- factor(mat_t$Labels)
levels(y) <- c("0", "1")
glmnet1 <- cv.glmnet(x=as.matrix(d), y=y, type.measure='class', nfolds=5, alpha=1, family="binomial")
c <- coef(glmnet1, s='lambda.min', exact=TRUE)
variables <- row.names(c)[which(c != 0)][-1]

# Generate violin plots for differential gene expression
df_long <- pivot_longer(seurat_phase_data_30_genes, cols = -Labels, names_to = "gene", values_to = "expression")
ggplot(df_long, aes(x = Labels, y = expression, fill = Labels)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
  stat_compare_means(method = "wilcox.test", label = "p.format", comparisons = list(c("extinct", "persist")), label.y = max(df_long$expression) + 2) +
  labs(title = "Combined Gene Expression by Clinical Outcome", x = "Clinical Outcome", y = "Gene Expression") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#Bulk RNA-seq analysis - all the processed files are already supplied in the repository under the "data" directory
GSE25055_processed = read.delim("GSE25055/GSE25055_processed.txt")
GSE25065_processed = read.delim("GSE25065/GSE25065_processed.txt")
GSE20271_processed = read.delim("GSE20271/GSE20271_processed.txt")
GSE20194_processed = read.delim("GSE20194/GSE20194_processed.txt")
#GSE32646_processed = read.delim("GSE32646/GSE32646_processed.txt")
#GSE41998_processed_ac_response = read.delim("GSE41998/GSE41998_processed_ac_response.txt")
#GSE41998_processed_pcr = read.delim("GSE41998/GSE41998_processed_pcr=0or1.txt")
#GSE50948_processed = read.delim("GSE50948/GSE50948_processed.txt")
#METABRIC_processed = read.delim("METABRIC/METABRIC_processed_TNBC_OS_5yrs_latest_Dec2023.txt")
all(setdiff(variables, colnames(GSE25055_processed)) == setdiff(variables, colnames(GSE25065_processed)))

all(setdiff(variables, colnames(GSE25055_processed)) == setdiff(variables, colnames(METABRIC_processed)))

length(variables)

variables_new_b_cells = variables[-c(which(variables %in% c( "MDH1B","POTEG","RSPO1","SMYD1","THRSP","ZFP57")))]

colnames(GSE25055_processed)[grep("HLA.DRB4", colnames(GSE25055_processed))] = "HLA-DRB4"

newrma_aggregrated=aggregate(GSE25055_processed[, -c(13286)],
                             by = list(Gene = rownames(GSE25055_processed)),
                             FUN = mean,
                             na.rm = TRUE)
newrma_aggregrated$Labels = GSE25055_processed$Labels
GSE25055_processed = newrma_aggregrated
setdiff(variables, colnames(GSE25055_processed))
write.table(GSE25055_processed, "GSE25055/GSE25055_processed.txt", col.names = TRUE, row.names = TRUE, sep="\t", quote = F)

colnames(METABRIC_processed)[grep("HLA.DRB4", colnames(METABRIC_processed))] = "HLA-DRB4"
which(colnames(METABRIC_processed) %in% "Labels")
newrma_aggregrated=aggregate(METABRIC_processed[, -c(20388)],
                             by = list(Gene = rownames(METABRIC_processed)),
                             FUN = mean,
                             na.rm = TRUE)
newrma_aggregrated$Labels = METABRIC_processed$Labels
METABRIC_processed = newrma_aggregrated
setdiff(variables, colnames(METABRIC_processed))
write.table(METABRIC_processed, "METABRIC/METABRIC_processed_TNBC_OS_5yrs_latest_Dec2023.txt", col.names = TRUE, row.names = TRUE, sep="\t", quote = F)

variables_new_b_cells = variables[-c(which(variables %in% c("MDH1B","POTEG","RSPO1","SMYD1","THRSP","ZFP57")))]
write.table(variables_new_b_cells, "/home/shayantan/TNBC/top_feats_B_cells.txt", col.names = T, row.names = F, sep="\t", quote = F)


setdiff(variables_new_b_cells, colnames(GSE20194_processed))
