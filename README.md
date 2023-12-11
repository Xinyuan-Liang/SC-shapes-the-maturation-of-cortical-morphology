# SC-shapes-the-maturation-of-cortical-morphology
This repository provides core code and relevant toolboxes for data analysis in the article entitled "Liang X, Sun L, Liao X, Lei T, Xia M, Duan D, Zeng Z, Xu Z, Men W, Wang Y, Tan S, Gao J, Qin S, Tao S, Dong Q, Zhao T, He Y, Structural connectome architecture shapes the maturation of cortical morphology from childhood to adolescence. bioRxiv, 2022, https://www.biorxiv.org/content/10.1101/2022.12.15.520527v2"

## Overview:

Content includes standalone software, source code, and demo data. All data required for reproducing our findings have been made publicly available, including regional cortical thickness matrix, structural connectivity matrix in three node parcellation resolutions for each participant, intermediate results during analysis, and data for visualizing the main figures. Raw data is available from the corresponding authors upon reasonable request. Of note, due to limitation of data size, some data is placed in another public cloud storage: https://pan.bnu.edu.cn/l/a1nOfo

## Installation guide:

Please use the “add path” method in MATLAB, "pip install" method in Python and "install.packages()" method in R to add toolboxes and scripts in the code folder. These procedures are not time-consuming.

## Code:

Following five analyses were performed using open-source packages:
1. **Spin test**: The spin test was conducted using an open Matlab code package (https://github.com/frantisekvasa/rotate_parcellation) (1, 2).

2. **Rewired test**: Surrogate networks that preserve the nodal degree and approximate edge length distribution of the empirical network were generated using a Matlab code package (https://www.brainnetworkslab.com/coderesources) (3). 

3. **Group SC generation**: Group-level connectome backbone was generated by using a consensus approach that preserves the connection length distributions of individual networks (https://www.brainnetworkslab.com/coderesources) (4).

4. **SVR model**: The cortical thickness (CT) maturation degree was predicted with an SVR model using an open Matlab code package (https://github.com/ZaixuCui/Pattern_Regression_Clean/tree/master/SVR) (5) and LIBSVM matlab toolbox (https://www.csie.ntu.edu.tw/~cjlin/libsvm/).

5. **Gene Ontology Enrichment Analysis**: Gene Ontology enrichment analysis was conducted using the ToppGene Suite (https://toppgene.cchmc.org/) (6).
6. **Visualization**: BrainNet Viewer (www.nitrc.org/projects/bnv) and netneurotools (https://github.com/netneurolab/netneurotools) toolbox.

The main analysis process was performed step-by-step using the following code.
1. **Mixed linear model**: s1_mixed_model_group_diff.m

2. **White matter network constraints on the spatial maturation of cortical thickness**: s2_Liang_corrNeighborCT.m s3_Liang_corrNeighborCT_system.m

3. **Generalized additive models were used to fit the maturation curves of nodal CT with age** : s4_GAM.R 

4. **Longitudinal analysis**: s5_lxy_sc_corr_longitudinal.m

5. **Random walk model**: s6_Liang_ramdom_walk.m

6. **Computing transcription level differences for 4 neurodevelopmental gene sets (7) between dominant and non-dominant regions**: s7_BrainSpan_dominant.m

7. **Preprocessing of AHBA gene expression data**: The regional microarray expression data were preprocessed using a recommended pipeline with the abagen toolbox (0.1.3) (https://github.com/rmarkello/abagen) (8, 9).
    s8_AHBA_gene_expression.py

    Custom code and toolboxes were tested on a 64-bit Windows 10 PC (Intel Core i7-9700, 56GB RAM). Data was analyzed using Matlab 2018b, Python 3.7.12, and R 3.6.3. We thank the authors and developers for providing these wonderful tools for data analysis.

## Data：

1.	Multiscale Desikan-Kiliany parcellation files (10) were downloaded using netneurotools toolbox (https://github.com/netneurolab/netneurotools).
2.	Basic information of participants: child_info.mat
3.	Cortical thickness  for each participant: cortical_thickness.mat
4.	Maturation extents (*t*-value between child and adolescent groups) and significances of cortical thickness at each parcellation resolution: 
   TVector.txt, PVector.txt
5.	White matter network matrix for each child: sc_child_indivi.mat
6.	Group-level connectome backbone at each parcellation resolution:
   Group_sc.mat, Group_sc_nonei.mat (For each node, we excluded all its spatially adjoining neighbors.)
7.	Matrix of the Euclidean distance between every two nodes: dist.mat
8.	CT maturation rates and white matter network for each individual with longitudinal scans: CT_info.mat, SC_individual.mat
9.	Transition probabilities matrix of the WM network: TP.mat
10.	The tissue samples in Brainspan dataset belong to dominant regions:
    dominant_regions.mat (‘0’ and ‘1’ denote non-dominant and dominant regions.)
11.	Neurodevelopmental processes related genes (7): neurodev_process.mat
12.	The relevant data for visualizing the figures are  provided in  https://pan.bnu.edu.cn/l/a1nOfo.

## References
1.	F. Vasa et al., Adolescent Tuning of Association Cortex in Human Structural Brain Networks. Cereb Cortex 28, 281-294 (2018).
2.	A. F. Alexander-Bloch et al., On testing for spatial correspondence between maps of human brain structure and function. Neuroimage 178, 540-551 (2018).
3.	R. F. Betzel, D. S. Bassett, Specificity and robustness of long-distance connections in weighted, interareal connectomes. Proceedings of the National Academy of Sciences 115, E4880-E4889 (2018).
4.	R. F. Betzel, A. Griffa, P. Hagmann, B. Mišić, Distance-dependent consensus thresholds for generating group-representative structural brain networks. Network neuroscience 3, 475-496 (2019).
5.	Z. Cui, G. Gong, The effect of machine learning regression algorithms and sample size on individualized behavioral prediction with functional connectivity features. Neuroimage 178, 622-637 (2018).
6.	J. Chen, E. E. Bardes, B. J. Aronow, A. G. Jegga, ToppGene Suite for gene list enrichment analysis and candidate gene prioritization. Nucleic acids research 37, W305-W311 (2009).
7. H. J. Kang et al., Spatio-temporal transcriptome of the human brain. Nature 478, 483-489 (2011).
8.	R. D. Markello et al., Standardizing workflows in imaging transcriptomics with the abagen toolbox. Elife 10, e72129 (2021).
9.	A. Arnatkevic̆iūtė, B. D. Fulcher, A. Fornito, A practical guide to linking brain-wide gene expression and neuroimaging data. Neuroimage 189, 353-367 (2019).	
10.	L. Cammoun et al., Mapping the human connectome at multiple scales with diffusion spectrum MRI. Journal of neuroscience methods 203, 386-397 (2012).



