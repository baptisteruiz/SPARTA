import pandas as pd
import csv
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler

import argparse

parser = argparse.ArgumentParser(description='A program that visualises the SoFA datasets as heatmaps splits the DeepMicro entries into train and test subsets')

parser.add_argument("-d", "--dataset_name", help="Name of the dataset used for this treatment")

parser.add_argument("-p", "--pipeline_path", help="Path to the whole pipeline folder")

parser.add_argument("-r", "--data_ref", help="Full reference to the dataset's folder")


parser.add_argument("--tfigm", action="store_true", help="Indicates whether the pipeline focuses on TF-IGM recalculated data")
parser.add_argument("--scaling", action="store_true", help="Indicates whether the pipeline focuses on scaled data")

args = parser.parse_args()

dataset_no_iter = args.dataset_name.split('_iteration')[0]


##FUNCTIONS
def create_test_train_subsets(full_dataset, indices):
    dataset_train = full_dataset.drop(indices, axis=0)
    dataset_test = full_dataset.iloc[indices,:]

    return(dataset_train, dataset_test)


##GET DATA (RANDOMLY CHOSEN INDICES, SoFAs AND LABELS)
with open(args.pipeline_path+'/SoFA_calculation/outputs/'+dataset_no_iter+'/test_indices.csv', newline='') as f:
    reader = csv.reader(f)
    indices = list(reader)[0]

indices = [int(i) for i in indices]

if args.tfigm:
    entree_dm_annots = pd.read_csv(args.pipeline_path+'/DeepMicro/data/entree_DeepMicro_'+args.dataset_name+'_tfigm.csv', header = None)
    

elif args.scaling:
    entree_dm_annots = pd.read_csv(args.pipeline_path+'/DeepMicro/data/entree_DeepMicro_'+args.dataset_name+'_scaled.csv', header = None)

else:
    entree_dm_annots = pd.read_csv(args.pipeline_path+'/SoFA_calculation/outputs/'+dataset_no_iter+'/entree_DeepMicro_'+args.dataset_name+'.csv', header = None)


entree_dm_otus = pd.read_csv(args.pipeline_path+'/SoFA_calculation/outputs/'+dataset_no_iter+'/entree_DeepMicro_'+args.dataset_name+'_OTU.csv', header = None)
labels = pd.read_csv(args.pipeline_path+'/Inputs/Label_'+dataset_no_iter+'.csv', header = None)


##CREATE HEATMAP REPRESENTATION OF THE SoFA PROFILE

scaler = MinMaxScaler()

annots_table = pd.DataFrame(scaler.fit_transform(entree_dm_annots)).T.astype('float64')


sample_name_labels = pd.read_csv(args.pipeline_path+'/SoFA_calculation/outputs/'+dataset_no_iter+'/score_annot_selon_otus_'+args.dataset_name+'.tsv', sep='\t', header = 0, index_col=0).columns
annot_name_labels = pd.read_csv(args.pipeline_path+'/SoFA_calculation/outputs/'+dataset_no_iter+'/score_annot_selon_otus_'+args.dataset_name+'.tsv', sep='\t', header = 0, index_col=0).index

label_table = pd.read_csv(args.pipeline_path+'/Inputs/Label_'+dataset_no_iter+'.csv', header = None)
label_table.index = sample_name_labels

annots_table.columns = sample_name_labels
annots_table.index = annot_name_labels

lut = dict(zip(label_table[0].unique(), "gr"))
col_colors = label_table[0].map(lut)

black_indices = [None for i in range(len(label_table.values))]
for k in indices:
    black_indices[k] = 'black'

df_colors = pd.DataFrame(data={'Control': col_colors[col_colors == 'g'], 'Sick': col_colors[col_colors == 'r'], 'Test_dataset': black_indices})


annots_table.to_csv(args.pipeline_path+'/Outputs/'+ args.data_ref +'/SoFAs/annots_table_temp_'+args.dataset_name+'.csv')
df_colors.to_csv(args.pipeline_path+'/Outputs/'+ args.data_ref +'/SoFAs/df_colors_temp_'+args.dataset_name+'.csv')

# clustergrid = sns.clustermap(annots_table, col_colors=df_colors, tree_kws={'colors':[colmap[s] for s in label_table[0].astype('str')]}, row_cluster=True, cmap = 'flare', method = 'ward', figsize=(20,20))
# plt.savefig(args.pipeline_path+'/Outputs/'+ args.data_ref +'/SoFAs/Heatmap_'+args.dataset_name+'.png')

##SEPARATE INTO TRAIN AND TEST


annots_train , annots_test = create_test_train_subsets(entree_dm_annots, indices)
otu_train , otu_test = create_test_train_subsets(entree_dm_otus, indices)
labels_train, labels_test = create_test_train_subsets(labels, indices)



if args.tfigm:
    annots_train.to_csv(args.pipeline_path+'/DeepMicro/data/entree_DeepMicro_'+args.dataset_name+'_tfigm.csv', header = None, index = None)
    annots_test.to_csv(args.pipeline_path+'/DeepMicro/data/entree_DeepMicro_'+args.dataset_name+'_test_tfigm.csv', header = None, index = None)
elif args.scaling:
    annots_train.to_csv(args.pipeline_path+'/DeepMicro/data/entree_DeepMicro_'+args.dataset_name+'_scaled.csv', header = None, index = None)
    annots_test.to_csv(args.pipeline_path+'/DeepMicro/data/entree_DeepMicro_'+args.dataset_name+'_test_scaled.csv', header = None, index = None)
else:
    annots_train.to_csv(args.pipeline_path+'/DeepMicro/data/entree_DeepMicro_'+args.dataset_name+'.csv', header = None, index = None)
    annots_test.to_csv(args.pipeline_path+'/DeepMicro/data/entree_DeepMicro_'+args.dataset_name+'_test.csv', header = None, index = None)


otu_train.to_csv(args.pipeline_path+'/DeepMicro/data/entree_DeepMicro_'+args.dataset_name+'_OTU.csv', header = None, index = None)
otu_test.to_csv(args.pipeline_path+'/DeepMicro/data/entree_DeepMicro_'+args.dataset_name+'_OTU_test.csv', header = None, index = None)

labels_train.to_csv(args.pipeline_path+'/DeepMicro/data/Label_'+dataset_no_iter+'.csv', header = None, index = None)
labels_test.to_csv(args.pipeline_path+'/DeepMicro/data/Label_'+dataset_no_iter+'_test.csv', header = None, index = None)
