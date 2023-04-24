import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
import sys
sys.setrecursionlimit(100000)

import argparse

parser = argparse.ArgumentParser(description='A program that visualises the SoFA datasets as heatmaps splits the DeepMicro entries into train and test subsets')

parser.add_argument("-d", "--dataset_name", help="Name of the dataset used for this treatment")

parser.add_argument("-p", "--pipeline_path", help="Path to the whole pipeline folder")

parser.add_argument("-r", "--data_ref", help="Full reference to the dataset's folder")

args = parser.parse_args()

dataset_no_iter = args.dataset_name.split('_iteration')[0]

label_table = pd.read_csv(args.pipeline_path+'/Inputs/Label_'+dataset_no_iter+'.csv', header = None)

annots_table = pd.read_csv(args.pipeline_path+'/Outputs/'+ args.data_ref +'/SoFAs/annots_table_temp_'+args.dataset_name+'.csv', index_col=0, header=0)
otu_table = pd.read_csv(args.pipeline_path+'/SoFA_calculation/outputs/'+ dataset_no_iter +'/'+args.dataset_name+'_stripped.tsv', sep ='\t', index_col=0, header=0)
df_colors = pd.read_csv(args.pipeline_path+'/Outputs/'+ args.data_ref +'/SoFAs/df_colors_temp_'+args.dataset_name+'_SoFA.csv', index_col=0)
df_colors_otu = pd.read_csv(args.pipeline_path+'/Outputs/'+ args.data_ref +'/SoFAs/df_colors_temp_'+args.dataset_name+'_OTU.csv', index_col=0)


scaler = MinMaxScaler()

otu_table_scaled = pd.DataFrame(scaler.fit_transform(otu_table.T.values)).T
otu_table_scaled.columns = otu_table.columns
otu_table_scaled.index = otu_table.index


colmap = {'0': (0, 0, 0, 0.7), '1': (0, 0, 0, 0.7)}

clustergrid = sns.clustermap(annots_table, col_colors=df_colors, tree_kws={'colors':[colmap[s] for s in label_table[0].astype('str')]}, row_cluster=True, cmap = 'flare', method = 'ward', figsize=(20,20))
plt.savefig(args.pipeline_path+'/Outputs/'+ args.data_ref +'/SoFAs/Heatmap_SoFAs_'+args.dataset_name+'.png')

clustergrid2 = sns.clustermap(otu_table_scaled, col_colors=df_colors_otu, tree_kws={'colors':[colmap[s] for s in label_table[0].astype('str')]}, row_cluster=True, cmap = 'flare', method = 'ward', figsize=(20,20))
plt.savefig(args.pipeline_path+'/Outputs/'+ args.data_ref +'/SoFAs/Heatmap_OTUs_'+args.dataset_name+'.png')