# SC-shapes-the-maturation-of-cortical-morphology

   This repository provides core code and relevant toolboxes for data analysis in the article "Liang X, Sun L, Liao X, Lei T, Xia M, Duan D, Zeng Z, Xu Z, Men W, Wang Y, Tan S, Gao J, Qin S, Tao S, Dong Q, Zhao T, He Y, Structural connectome architecture shapes the maturation of cortical morphology from childhood to adolescence. bioRxiv, 2022, https://www.biorxiv.org/content/10.1101/2022.12.15.520527v2"

   ## Overview:

   Contents include standalone software, source code, and demo data. All data required for reproducing our findings have been made publicly available, including regional cortical thickness (CT) matrices, structural connectivity (SC) matrices, intermediate results during analysis, and data for visualizing the main figures.  All data are placed in Releases (https://github.com/Xinyuan-Liang/SC-shapes-the-maturation-of-cortical-morphology/releases/tag/v1.0.0).

   ## Installation guide:

   Please use the “add path” method in Matlab, "pip install" method in Python and "install.packages()" method in R to add toolboxes and scripts in the code folder. These procedures are not time-consuming.

   ## Code:

   Te following analyses were performed using open-source packages:

   1. **Spin test**: The spin test was conducted using an open Matlab code package (https://github.com/frantisekvasa/rotate_parcellation) [1, 2].

   2. **Rewired test**: Surrogate networks that preserve the nodal degree and approximate edge length distribution of the empirical network were generated using a Matlab code package (https://www.brainnetworkslab.com/coderesources) [3]. 

   3. **Group SC generation**: Group-level connectome backbone was generated by using a consensus approach that preserves the connection length distributions of individual networks (https://www.brainnetworkslab.com/coderesources) [4].

   4. **Support vector regression (SVR) model**: The CT maturation degree was predicted with an SVR model with nodal diffusive profiles at each neighboring scale separately as input features. This analysis was performed using an open Matlab code package (https://github.com/ZaixuCui/Pattern_Regression_Clean/tree/master/SVR) [5] and LIBSVM Matlab toolbox (https://www.csie.ntu.edu.tw/~cjlin/libsvm/).

   5. **Gene Ontology Enrichment Analysis**: Gene Ontology enrichment analysis was conducted using online tools including ToppGene Suite (https://toppgene.cchmc.org/) [6] and REViGO (http://revigo.irb.hr).
   6. **Visualization**: BrainNet Viewer (www.nitrc.org/projects/bnv) and netneurotools (https://github.com/netneurolab/netneurotools) toolbox.

   The main analyses were performed step-by-step using the source code as follows.

   1. **Mixed linear model**: s1_mixed_model_group_diff.m

   2. **Calculated the spatial correlation between the empirical CT maturation extent and the mean of its directly connected neighbor nodes**: s2_Liang_corrNeighborCT.m 

   3. **Calculated the spatial correlation between the empirical CT maturation extent and the mean of its directly connected neighbor nodes** (excluded all spatially adjoining neighbors): s3_Liang_corrNeighborCT_system.m

   4. **Generalized additive models (GAM) were used to fit the maturation curves of nodal CT with age**: s4_GAM.R 

   5. **Longitudinal analysis**: s5_lxy_sc_corr_longitudinal.m

   6. **Random walk model**: s6_Liang_ramdom_walk.m

   7. **Computing transcription level differences for 4 neurodevelopmental gene sets [7] between dominant and non-dominant regions**: s7_BrainSpan_dominant.m

   8. **Preprocessing of Allen Human Brain Atlas (AHBA) gene expression data**: The regional microarray expression data were preprocessed using a recommended pipeline with the abagen toolbox (0.1.3) (https://github.com/rmarkello/abagen) [8, 9].
      s8_AHBA_gene_expression.py

      Custom code and toolboxes were tested on a 64-bit Windows 10 PC (Intel Core i7-9700, 56GB RAM). Data was analyzed using Matlab 2018b, Python 3.7.12, and R 3.6.3. We thank the authors and developers for providing these tools for data analysis.

   ## Data：

   For the 1-7 items, we provide data in the CBDP dataset at three resolutions (219, 448, and 1000 nodes), and in the HCPD dataset at the highest resolution (1000 nodes).

   1.	Basic information of participants: child_info.mat
   2.	Regional cortical thickness for each participant: cortical_thickness.mat
   3.	Maturation extents (*t*-value between child and adolescent groups) and significances of cortical thickness: TVector.txt, PVector.txt
   4.	White matter network matrix for each participant in child group: sc_child_indivi.mat
   5.	Group-level connectome backbone: Group_sc.mat (raw SC), Group_sc_nonei.mat (SC excluded all spatially adjoining neighbors).
   6.	Matrix of the Euclidean distance between every two nodes: dist.mat
   7.	Transition probabilities matrix of the WM network: TP.mat

   For the GAM and longitudinal analyses (items 8-9), we provide following data in the CBDP dataset at the highest resolution (1000 nodes).

   8.	CT maturation rates derived from GAM analysis: CTmaturation_rate_GAM.mat
   9.	CT maturation rates and white matter network for each individual with longitudinal scans (1000 nodes): CTmaturation_rate_longitudinal.mat, SC_individual_longitudinal.mat

   For the gene analysis (items 10-12), we provide following data.

   10.	The tissue samples in Brainspan dataset belong to dominant regions: tissue_dominant_regions.mat (‘0’ and ‘1’ denote non-dominant and dominant regions.)
   11.	Neurodevelopmental processes related genes [7]: neurodev_process.mat
   12.	Gene expression information of the neocortex from the Brainspan dataset: Brainspan_neocortex.mat

   Finally, some supporting data (items 13-14) can be found as follows.

   13.	Multiscale Desikan-Kiliany parcellation files [10] were downloaded using netneurotools toolbox (https://github.com/netneurolab/netneurotools).
   14.	The data for visualizing the figures are  provided in releases https://github.com/Xinyuan-Liang/SC-shapes-the-maturation-of-cortical-morphology/releases/tag/v1.0.0

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


