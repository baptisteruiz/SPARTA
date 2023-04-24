import pandas as pd
import csv
import joblib
from sklearn.metrics import roc_auc_score
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='A program that raandomly extracts a portion of the datasets to be used as benchmarking tests')


parser.add_argument("-d","--dataset_name", help="Name of the dataset used for this treatment")

parser.add_argument("-p", "--pipeline_path", help="Path to the whole pipeline folder")

parser.add_argument("-r", "--data_ref", help="Full reference to the dataset's folder")

parser.add_argument("-i", "--iteration", help="Reference for the method's iteration")

parser.add_argument("--tfigm", action="store_true", help="Indicates whether the pipeline focuses on TF-IGM recalculated data")
parser.add_argument("--scaling", action="store_true", help="Indicates whether the pipeline focuses on scaled data")

args = parser.parse_args()

print('Dataset name :', args.dataset_name)
dataset_name_no_otu = args.dataset_name.split('_OTU')[0]
print('Dataset name no otu:', dataset_name_no_otu)
dataset_name_no_otu_no_iteration = dataset_name_no_otu.split('_iteration')[0]
print('Dataset name no otu no iteration:', dataset_name_no_otu)
if args.dataset_name[-4:] == '_OTU':
    otu = '_OTU'
    otu_bool = True
else:
    otu = ''
    otu_bool = False
dataset_name_no_iteration = dataset_name_no_otu_no_iteration + otu
print('OTU bool', otu_bool)

##IMPORT LAST TRAINED MODEL
#clf = joblib.load(args.pipeline_path+'/DeepMicro/saved_classifier.joblib')

if args.tfigm:
    raw = pd.read_csv(args.pipeline_path+'/DeepMicro/data/entree_DeepMicro_'+args.dataset_name+'_test_tfigm.csv', header=None)
    clf = joblib.load(args.pipeline_path+'/DeepMicro/entree_DeepMicro_'+args.dataset_name+'_tfigm_saved_classifier.joblib')
    with open(args.pipeline_path+'/DeepMicro/results/entree_DeepMicro_'+args.dataset_name+'_tfigm_threshold_opt.txt', 'r') as f:
        threshold = float(f.readlines()[0])
elif args.scaling:
    raw = pd.read_csv(args.pipeline_path+'/DeepMicro/data/entree_DeepMicro_'+args.dataset_name+'_test_scaled.csv', header=None)
    clf = joblib.load(args.pipeline_path+'/DeepMicro/entree_DeepMicro_'+args.dataset_name+'_scaled_saved_classifier.joblib')
    with open(args.pipeline_path+'/DeepMicro/results/entree_DeepMicro_'+args.dataset_name+'_scaled_threshold_opt.txt', 'r') as f:
        threshold = float(f.readlines()[0])
else:
    raw = pd.read_csv(args.pipeline_path+'/DeepMicro/data/entree_DeepMicro_'+args.dataset_name+'_test.csv', header=None)
    clf = joblib.load(args.pipeline_path+'/DeepMicro/entree_DeepMicro_'+args.dataset_name+'_saved_classifier.joblib')
    with open(args.pipeline_path+'/DeepMicro/results/entree_DeepMicro_'+args.dataset_name+'_threshold_opt.txt', 'r') as f:
        threshold = float(f.readlines()[0])

Xtest = raw.values.astype(np.float64)

label = pd.read_csv(args.pipeline_path+'/DeepMicro/data/Label_'+dataset_name_no_otu_no_iteration+'_test.csv', header = None)
Ytest = label.values.astype(int)

predictions = clf.predict_proba(Xtest)
# auc_score = round(roc_auc_score(Ytest, predictions[:, 1]), 4)
# print("auc_score", str(auc_score))


# with open(args.pipeline_path+'/Outputs/'+ args.data_ref +'/Classification_performances/Data_leak_free_'+dataset_name_no_iteration+'_performances.txt', 'a') as f:

#     f.writelines('Run '+ args.iteration +' (AUC): '+str(auc_score) +'\n')

pred_categ = (clf.predict_proba(Xtest)[:,1] >= threshold).astype('int')
print("Ytest : ", Ytest)
print("pred : ", pred_categ)

with open(args.pipeline_path+'/SoFA_calculation/outputs/'+dataset_name_no_otu_no_iteration+'/test_indices.csv', newline='') as f:
    reader = csv.reader(f)
    indices = list(reader)[0]

indices = [int(i) for i in indices]
print("Indices : ", indices)

df_colors = pd.read_csv(args.pipeline_path+'/Outputs/'+ args.data_ref +'/SoFAs/df_colors_temp_'+dataset_name_no_otu+'.csv', index_col=0)

for i in range(len(indices)):
    if pred_categ[i] != Ytest [i]:
        df_colors.iloc[indices[i],:]['Test_dataset'] = 'purple'

if otu_bool:
    df_colors.to_csv(args.pipeline_path+'/Outputs/'+ args.data_ref +'/SoFAs/df_colors_temp_'+dataset_name_no_otu+'_OTU.csv')
    print('Out path: ', args.pipeline_path+'/Outputs/'+ args.data_ref +'/SoFAs/df_colors_temp_'+dataset_name_no_otu+'_OTU.csv')
else:
    df_colors.to_csv(args.pipeline_path+'/Outputs/'+ args.data_ref +'/SoFAs/df_colors_temp_'+dataset_name_no_otu+'_SoFA.csv')
    print('Out path: ', args.pipeline_path+'/Outputs/'+ args.data_ref +'/SoFAs/df_colors_temp_'+dataset_name_no_otu+'_SoFA.csv')