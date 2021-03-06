import numpy as np
import networkx as nx
from scipy import stats as st
import os

cases = '/home/senthilp/caesar/camcan/cc700/freesurfer_output/full.txt'
SUBDIR = '/home/senthilp/caesar/camcan/cc700/freesurfer_output'
n = 268
#freqList = [2,4,6,8,10,12,16,20,24,32,64]
freqList = [4]

with open(cases) as f:
     case_list = f.read().splitlines()

for subject in case_list:

    print(f'Processing subject... {subject}')
    for freq in freqList:
        if os.path.exists(os.path.join(SUBDIR,subject,'mri','shen_corr',f'shen_corr_{freq}Hz_lcmv.npy')):        
            if not os.path.exists(os.path.join(SUBDIR,subject,'mri','shen_corr',f'{subject}_cluster_{freq}Hz.npy')):

                print(f'subject... {subject}, freq {freq}')
                connectivity = np.load(f'{SUBDIR}/{subject}/mri/shen_corr/shen_corr_{freq}Hz_lcmv.npy')

                connectivity = np.array(connectivity)
                result = np.zeros((connectivity.shape))
                total_connection = len(connectivity) - 1

                for i, row in enumerate(connectivity):
        
                    for j, value in enumerate(row):
                
           # '''
           # * 'less': the mean of the underlying distribution of the sample is
            #      significantly less than the given population mean (`value`)
           # '''
                        stat1, pvalue1 = st.ttest_1samp(row, value, alternative='less',nan_policy='omit')
                        stat2, pvalue2 = st.ttest_1samp(connectivity[j], value, alternative='less',nan_policy='omit')
        
            #print(pvalue1, row.mean(), value)
            #print(pvalue2, connectivity[j].mean(), value)
            
                        if pvalue1 < 0.01 and pvalue2 < 0.01:
            
                            result[i][j] = 1
    
                degree_count = np.sum(result, axis=1)
                degree_avg = degree_count / total_connection
                np.save(f'{SUBDIR}/{subject}/mri/shen_corr/{subject}_degree_{freq}Hz.npy', degree_avg)

                G = nx.from_numpy_array(result, create_using=nx.DiGraph)
                btw = nx.betweenness_centrality(G, normalized=True)
                btw_npy = np.array(list(btw.values()))
                np.save(f'{SUBDIR}/{subject}/mri/shen_corr/{subject}_betweeness_{freq}Hz.npy', btw_npy)

    # Shortest path average
                all_pairs = nx.floyd_warshall(G)
                s = [sum(t.values()) for t in all_pairs.values()]
                short_path_avg = np.array(s)/n
                np.save(f'{SUBDIR}/{subject}/mri/shen_corr/{subject}_shortPath_{freq}Hz.npy', short_path_avg)

    # Cluster
                cluster_coeff = list(nx.clustering(G).values())
                cluster_coeff = np.array(cluster_coeff)
                np.save(f'{SUBDIR}/{subject}/mri/shen_corr/{subject}_cluster_{freq}Hz.npy', cluster_coeff)

                del G, btw, btw_npy, connectivity, result, total_connection, degree_avg, short_path_avg, cluster_coeff

