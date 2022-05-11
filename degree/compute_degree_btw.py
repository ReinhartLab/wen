import numpy as np
import networkx as nx
from scipy import stats as st

cases = '/home/senthilp/caesar/camcan/cc700/freesurfer_output/12.txt'
SUBDIR = '/home/senthilp/caesar/camcan/cc700/freesurfer_output'
with open(cases) as f:
     case_list = f.read().splitlines()

for subject in case_list:

    print(f'Processing subject... {subject}')

    connectivity = np.load(f'{SUBDIR}/{subject}/mri/shen_corr/shen_corr_5Hz_lcmv.npy')

    connectivity = np.array(connectivity)
    result = np.zeros((connectivity.shape))
    total_connection = len(connectivity) - 1

    for i, row in enumerate(connectivity):
        
        for j, value in enumerate(row):
                
            '''
            * 'less': the mean of the underlying distribution of the sample is
                  significantly less than the given population mean (`value`)
            '''
            stat1, pvalue1 = st.ttest_1samp(row, value, alternative='less')
            stat2, pvalue2 = st.ttest_1samp(connectivity[j], value, alternative='less')
        
            #print(pvalue1, row.mean(), value)
            #print(pvalue2, connectivity[j].mean(), value)
        
            if pvalue1 < 0.01 and pvalue2 < 0.01:
            
                result[i][j] = 1
    
    degree_count = np.sum(result, axis=1)
    degree_avg = degree_count / total_connection
    np.save(f'{subject}_degree.npy', degree_avg)

    G = nx.from_numpy_array(result, create_using=nx.DiGraph)
    btw = nx.betweenness_centrality(G, normalized=True)
    btw_npy = np.array(list(btw.values()))
    np.save(f'{subject}_betweeness.npy', btw_npy)

    del G, btw, btw_npy, connectivity, result, total_connection
