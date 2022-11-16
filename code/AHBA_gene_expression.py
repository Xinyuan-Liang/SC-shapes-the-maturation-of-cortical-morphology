# Preprocessed the microarray expression data from AHBA datasets with the abagen toolbox (https://github.com/rmarkello/abagen)
# Written by Xinyuan Liang, xyliang@mail.bnu.edu.cn
# State Key Laboratory of Cognitive Neuroscience and Learning &
# IDG/McGovern Institute of Brain Research, 
# Beijing Normal University,
# Beijing, PR China.

import abagen
import pandas as pd
import numpy as np

#probe selection differential stability, differential stability 0.1
roinum = '125'
expression_out = abagen.get_expression_data('/Data/xyliang/gene_expression/Camm'+roinum+'_CBDPatlas_mni.nii', \
atlas_info = '/Data/xyliang/gene_expression/atl-Cammoun2012_space-MNI152NLin2009aSym_info_gene_'+roinum+'.csv', \
return_donors = True,ibf_threshold=0.5,data_dir='/Data/xyliang/AHBA')
expression,stability = abagen.keep_stable_genes(list(expression_out.values()),threshold=0.1,percentile=False,return_stability=True)
expression = pd.concat(expression).groupby('label').mean() 
np.savetxt('/Data/xyliang/gene_expression/C'+roinum+'_stability.csv', stability, delimiter=',')
expression.to_csv('/Data/xyliang/gene_expression/C'+roinum+'_DS01_expression.csv')



