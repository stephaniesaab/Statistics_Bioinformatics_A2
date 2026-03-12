# Statistics_Bioinformatics_A2
Assignment 2 for Statistics for Bioinformatics

### Question 1.
Principal component analysis (PCA) was performed on the gene expression dataset containing measurements for the 156 most differentially expressed genes across T cell samples. PCA
reduces the dimensionality of the dataset while preserving the major sources of variation in gene expression.

## Visualisation
### Heatmap

The heat map shows naïve (NAI) samples clustering together and show high expression of naïve‑associated genes such as CCR7 and SELL, while effector (EFFE) samples form a distinct cluster characterized by elevated expression of effector‑related genes like granzyme B and HLA‑DR. Memory (MEM) samples consistently appear between these two groups in the dendrogram, reflecting their intermediate transcriptional profile and aligning with the paper’s conclusion that MEM cells occupy a transitional state between NAI and EFFE. These results visually confirm the strong separation of NAI and EFFE cells, the intermediate nature of MEM cells, and the presence of gene modules that drive these distinctions, directly mirroring the lineage relationships described in the original study.

<img width="1200" height="700" alt="q1 p2" src="https://github.com/user-attachments/assets/981b613c-4a79-4feb-848d-4fdd8146e808" />

### Scree plot

The scree plot indicated that the first two principal components explain a substantial proportion of the total variance and captures significant biological differences between the 3 cell types. This corresponds with the observation made by Holmes et al. (2005) about 80% of the variance falling within the first two PCs.  This suggests that the high-dimensional dataset can be effectively summarized using a small number of components and justifies dimensionality reduction as the lower variance PCAs will contribute to noise.  

<img width="1000" height="304" alt="q1 p4" src="https://github.com/user-attachments/assets/b173a583-b182-4371-83a9-5cacd275304a" />


## Principal Component Analysis
### PCA

Visualization of the samples in the PC1–PC2 space revealed clear clustering according to T cell type. Naïve (NAI), effector (EFFE), and memory (MEM) T cells form distinct groups, indicating systematic differences in their gene expression profiles. Memory T cells appear positioned between naïve and effector cells along the principal component axes, suggesting that their gene expression patterns are intermediate between these two states. This observation is consistent with findings reported by Holmes et al. (2005), which proposed that memory T cells exhibit transcriptional profiles between naïve and effector T cells, favoring the parallel differentiation model over the linear differentiation model

<img width="1500" height="700" alt="q1 p5" src="https://github.com/user-attachments/assets/0a6a7c2a-e251-401b-ad83-452cf50bc597" />

### PCA genes loading

The PCA loading plot of genes highlights those genes contributing most strongly to the principal components. Genes located farther from the origin have larger loadings, which plays a
greater role in driving the observed variation among samples. Based on the previous PCA plot, genes at the far left are associated with upregulation in Effector cells, genes on the far right are associated with upregulation in Naive T cells and genes found near the origin indicate a weak differential expression between groups. 

<img width="1000" height="500" alt="q1 p6" src="https://github.com/user-attachments/assets/224804bd-24b8-4160-bf66-3c0dd5ccdfbb" />

The results of the statistical analysis replicate and support the paper’s conclusion that memory T cells represent a transitional transcriptional state, providing strong evidence for a **parallel differentiation model** in which naïve cells give rise to both effector and memory populations through a shared intermediate state.

