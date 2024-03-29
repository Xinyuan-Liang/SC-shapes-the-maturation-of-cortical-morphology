# SC-shapes-the-maturation-of-cortical-morphology

   This repository provides core code and relevant toolboxes for data analysis in the article "Liang X, Sun L, Liao X, Lei T, Xia M, Duan D, Zeng Z, Li Q, Xu Z, Men W, Wang Y, Tan S, Gao J, Qin S, Tao S, Dong Q, Zhao T, He Y (2024) Structural connectome architecture shapes the maturation of cortical morphology from childhood to adolescence. Nat Commun. 15,784. https://doi.org/10.1038/s41467-024-44863-6"

   ## Overview:

   Contents include standalone software, source code, and demo data. All data necessary to reproduce our results have been made publicly available, including regional cortical thickness (CT) matrices, structural connectivity (SC) matrices, intermediate results during the analysis, and data for visualizing the figures. All data are available in Releases (https://github.com/Xinyuan-Liang/SC-shapes-the-maturation-of-cortical-morphology/releases/tag/v1.0.0).

   ## Installation guide:

   All the scripts (*.m files in the code folder) and the toolboxes can be executed by adding the appropriate environment paths according to the corresponding code language. Use the "add path" function for Matlab, the "pip install" function for Python, and the "install.packages()" function for R. These procedures are not time-consuming.

   ## Code:

   The following analyses were carried out using open source packages:

   1. **Spin test**: The spin test was conducted using an open Matlab code package (https://github.com/frantisekvasa/rotate_parcellation) [1, 2].
   2. **Rewired test**: Surrogate networks that preserve the nodal degree and approximate edge length distribution of the empirical network were generated using a Matlab code package (https://www.brainnetworkslab.com/coderesources) [3]. 
   3. **Group SC generation**: The group-level connectome backbone was generated by using a consensus approach that preserves the connection length distributions of individual networks (https://www.brainnetworkslab.com/coderesources) [4].
   4. **Support vector regression (SVR) model**: The CT maturation degree was predicted using an SVR model with nodal diffusive profiles at each neighboring scale separately as input features. This analysis was performed using an open Matlab code package (https://github.com/ZaixuCui/Pattern_Regression_Clean/tree/master/SVR) [5] and LIBSVM Matlab toolbox (https://www.csie.ntu.edu.tw/~cjlin/libsvm/).
   5. **Generalized additive models (GAM)**: The generalized additive models were performed using the R package mgcv (https://cran.r-project.org/web/packages/mgcv/index.html).
   6. **Allen Human Brain Atlas (AHBA) datasets**: The AHBA datasets were preprocessed using the Python toolbox abagen (https://github.com/rmarkello/abagen).
   7. **Gene Ontology Enrichment Analysis**: Gene ontology enrichment analysis was conducted using online tools including ToppGene Suite (https://toppgene.cchmc.org/) [6] and REViGO (http://revigo.irb.hr).
   8. **Visualization**: BrainNet Viewer (www.nitrc.org/projects/bnv) and netneurotools toolbox(https://github.com/netneurolab/netneurotools).

   The main analyses were carried out step-by-step using the source code as follows.

   1. **Mixed linear model**: s1_mixed_model_group_diff.m

   2. **Calculate the spatial correlation between the empirical CT maturation extent and the mean of its directly connected neighbor nodes**: s2_Liang_corrNeighborCT.m 

   3. **Calculate the spatial correlation between the empirical CT maturation extent and the mean of its directly connected neighbor nodes** (excluding all spatially adjoining neighbors): s3_Liang_corrNeighborCT_system.m

   4. **Generalized additive models (GAM) were used to fit the maturation curves of nodal CT with age**: s4_GAM.R 

   5. **Longitudinal analysis**: s5_lxy_sc_corr_longitudinal.m

   6. **Random walk model**: s6_Liang_ramdom_walk.m

   7. **Calculation of transcription level differences for 4 neurodevelopmental gene sets [7] between dominant and non-dominant regions**: s7_BrainSpan_dominant.m

   8. **Preprocessing of gene expression data from the Allen Human Brain Atlas (AHBA)**: Regional microarray expression data were preprocessed using a recommended pipeline with the abagen toolbox (0.1.3) (https://github.com/rmarkello/abagen) [8, 9].
      s8_AHBA_gene_expression.py

Custom code and toolboxes were tested on a 64-bit Windows 10 PC (Intel Core i7-9700, 56GB RAM). Data analysis was performed using Matlab 2018b, Python 3.7.12, and R 3.6.3. We thank the authors and developers for providing these data analysis tools.

   ## Data：
   All data are available in Releases (https://github.com/Xinyuan-Liang/SC-shapes-the-maturation-of-cortical-morphology/releases/tag/v1.0.0).
   
   For items 1-7, we provided data at three resolutions (219, 448, and 1000 nodes) for the CBDP dataset and data at the highest resolution (1000 nodes) for the HCPD dataset.

   1.	Basic information of participants: child_info.mat
   2.	Regional cortical thickness for each participant: cortical_thickness.mat
   3.	Maturation extents (*t*-value between child and adolescent groups) and significances of cortical thickness: TVector.txt, PVector.txt
   4.	White matter network matrix for each participant in the child group: sc_child_indivi.mat
   5.	Connectome backbone at group level: Group_sc.mat (raw SC), Group_sc_nonei.mat (SC excluding all spatially adjoining neighbors).
   6.	Matrix of the Euclidean distance between every two nodes: dist.mat
   7.	Transition probability matrices of the WM network: TP.mat

   For the GAM and longitudinal analyses (items 8-9), we provide the following data in the CBDP dataset at the highest resolution (1000 nodes).

   8.	CT maturation rates derived from the GAM analysis: CTmaturation_rate_GAM.mat
   9.	CT maturation rates and white matter network for each individual with longitudinal scans: CTmaturation_rate_longitudinal.mat, SC_individual_longitudinal.mat

   For the gene-based analysis (items 10-12), we provide the following data.

   10.	The tissue samples in the Brainspan dataset belong to dominant regions: tissue_dominant_regions.mat (‘0’denotes non-dominant regions and‘1’denotes dominant regions.)
   11.	Genes associated with neurodevelopmental processes [7]: neurodev_process.mat
   12.	Gene expression data of the neocortex from the Brainspan dataset: Brainspan_neocortex.mat
   13.	The detailed enrichment analysis results for all genes: AHBA_results.xlsx

   Finally, some supporting data (items 14-15) can be found as follows.

   14.	Multiscale Desikan-Kiliany parcellation files [10] were downloaded using the netneurotools toolbox (https://github.com/netneurolab/netneurotools).
   15.	The data used to visualize the figures are available at https://github.com/Xinyuan-Liang/SC-shapes-the-maturation-of-cortical-morphology/releases/tag/v1.0.0

   ## References

   1.	F. Vasa et al., Adolescent Tuning of Association Cortex in Human Structural Brain Networks. Cereb Cortex 28, 281-294 (2018).
   2.	A. F. Alexander-Bloch et al., On testing for spatial correspondence between maps of human brain structure and function. Neuroimage 178, 540-551 (2018).
   3.	R. F. Betzel, D. S. Bassett, Specificity and robustness of long-distance connections in weighted, interareal connectomes. Proceedings of the National Academy of Sciences 115, E4880-E4889 (2018).
   4.	R. F. Betzel, A. Griffa, P. Hagmann, B. Mišić, Distance-dependent consensus thresholds for generating group-representative structural brain networks. Network neuroscience 3, 475-496 (2019).
   5.	Z. Cui, G. Gong, The effect of machine learning regression algorithms and sample size on individualized behavioral prediction with functional connectivity features. Neuroimage 178, 622-637 (2018).
   6.	J. Chen, E. E. Bardes, B. J. Aronow, A. G. Jegga, ToppGene Suite for gene list enrichment analysis and candidate gene prioritization. Nucleic acids research 37, W305-W311 (2009).
   7.	H. J. Kang et al., Spatio-temporal transcriptome of the human brain. Nature 478, 483-489 (2011).
   8.	R. D. Markello et al., Standardizing workflows in imaging transcriptomics with the abagen toolbox. Elife 10, e72129 (2021).
   9.	A. Arnatkevic̆iūtė, B. D. Fulcher, A. Fornito, A practical guide to linking brain-wide gene expression and neuroimaging data. Neuroimage 189, 353-367 (2019).	
   10.	L. Cammoun et al., Mapping the human connectome at multiple scales with diffusion spectrum MRI. Journal of neuroscience methods 203, 386-397 (2012).


