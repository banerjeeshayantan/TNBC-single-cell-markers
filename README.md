# TNBC-single-cell-markers
Reproducible code for identifying TIL-based markers from single-cell RNA-seq data of TNBC patients

### Abstract
Triple-negative breast cancer (TNBC) is an aggressive subtype often marked by resistance to neoadjuvant chemotherapy (NAC), making treatment particularly challenging. Tumor-infiltrating lymphocytes (TILs), crucial players in the immune landscape of tumors, have been associated with treatment outcomes, but the prognostic potential of TIL-derived gene markers in pre-NAC samples from TNBC patients remains understudied. In this research, we analyzed the single-cell transcriptional profiles of approximately 5,000 cells from four chemosensitive and four chemoresistant TNBC patients using publicly available datasets. Leveraging standard single-cell analysis, we identified differentially expressed gene signatures within the TIL subpopulation, highlighting significant immune activation pathways differentiating chemoresistant from chemosensitive tumors. By employing robust feature selection and repeated cross-validation across microarray and RNA-seq datasets, we developed a stable set of 30 TIL-based gene markers with notable prognostic relevance for NAC response in TNBC. These markers achieved an AUROC of 0.78 in the training set and validated with AUROCs of 0.8, 0.658, and 0.736 across five independent test datasets, demonstrating consistency across diverse platforms and sequencing technologies. Furthermore, increased expression of these gene signatures correlated with improved recurrence-free survival (RFS) in a cohort of 220 TNBC patients. This study enhances our understanding of the TIL transcriptional landscape in NAC response, identifying potential biomarkers and therapeutic targets for improving treatment outcomes in TNBC.


### Datasets required to run the codes
All associated datasets required to run the R codes and the python notebooks have been deposited in Zenodo: https://doi.org/10.5281/zenodo.14789164


