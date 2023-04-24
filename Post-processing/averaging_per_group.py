import pandas as pd
import numpy as np

import argparse

parser = argparse.ArgumentParser(description='A program that averages the RF importance scores of the functional annotations, and associates them to OTUs')

parser.add_argument("-d","--dataset_name", help="Name of the dataset used for this treatment")

parser.add_argument("-p", "--pipeline_path", help="Path to the whole pipeline folder")

parser.add_argument("--tfigm", action="store_true", help="Indicates whether the pipeline focuses on TF-IGM recalculated data")
parser.add_argument("--scaling", action="store_true", help="Indicates whether the pipeline focuses on scaled data")

args = parser.parse_args()

dataset_name = args.dataset_name
pipeline_path = args.pipeline_path

###STEP 1: calculate the average score of each annotation for the designated group
def average_per_group(score_dataframe, label_refs, label_value = None):

    if label_value != None:
        table_group = score_dataframe.loc[:, [i[0] == label_value for i in label_refs.values]]
    else:
        table_group = score_dataframe

    index_filters = []
    for i in list(score_dataframe.index):
        if ((i[0].isnumeric()) or ("GO" in i)) and (i != '16s_rrna'):
            index_filters.append(i)
    

    filtered_table_group = table_group.loc[index_filters].astype(float)

    filtered_table_group["average"] = filtered_table_group.mean(axis=1)


    return(filtered_table_group)


###STEP 2: calculate the average prevalence of each OTU in each group

def average_count_per_group(count_dataframe,  label_refs, label_value = None):

    if label_value != None:
        count_group = count_dataframe.loc[:, [i[0] == label_value for i in label_refs.values]].astype(float)
    else:
        count_group = count_dataframe.astype(float)

    count_group["average"] = count_group.mean(axis=1)


    return(count_group)


##STEP 3: separate the average scores into significant/non significant categories

def separate_significant(dataframe_work, dataframe_ref_significance):


    cutoff_index = dataframe_ref_significance.index[dataframe_ref_significance["ID"] == 'CUTOFF']
    list_significant = list(dataframe_ref_significance["ID"][:cutoff_index[0]] )

    dataframe_significant = dataframe_work.reindex(list_significant)
    list_significant.append('CUTOFF')
    dataframe_non_significant = dataframe_work[~dataframe_work.index.isin(list_significant)]

    return(dataframe_significant, dataframe_non_significant)



label_refs = pd.read_csv(pipeline_path + '/Inputs/Label_'+dataset_name+'.csv', header = None)

score_db = pd.read_csv(pipeline_path + '/SoFA_calculation/outputs/'+dataset_name+'/score_annot_selon_otus_'+dataset_name+'.tsv', sep ="\t", header = 0, index_col=0)

avg_score_total = average_per_group(score_db, label_refs)
avg_score_total.index = [ind.replace("_",".") for ind in list(avg_score_total.index)]
dict_averaged_score_dfs = {}


count_db = pd.read_csv(pipeline_path + '/SoFA_calculation/outputs/'+dataset_name+'/'+dataset_name+'_stripped.tsv', sep ="\t", header =0, index_col=0)

avg_count_total = average_count_per_group(count_db, label_refs)
dict_averaged_count_dfs = {}

df_max_scores = pd.DataFrame()
df_max_counts = pd.DataFrame()

for label in np.unique(label_refs.values):
    avg_score_group = average_per_group(score_db, label_refs, label)
    avg_count_group = average_count_per_group(count_db, label_refs, label)

    avg_score_group.index = [ind.replace("_",".") for ind in list(avg_score_group.index)]

    dict_averaged_score_dfs[str(label)] = avg_score_group
    dict_averaged_count_dfs[str(label)] = avg_count_group

    df_max_scores[label] = avg_score_group['average']
    df_max_counts[label] = avg_count_group['average']


# max_score_groups = []

# for index in list(avg_score_total.index):

#     index_max = list(dict_averaged_score_dfs.keys())[0]
#     for key in list(dict_averaged_score_dfs.keys()):

#         if dict_averaged_score_dfs[key]['average'][index] > dict_averaged_score_dfs[index_max]['average'][index]:
#             index_max = index
#     max_score_groups.append(index_max)
avg_score_total['Representative_group'] = df_max_scores.idxmax(axis=1).values

# max_count_groups = []
# for index in list(avg_count_total.index):
#     index_max = list(dict_averaged_count_dfs.keys())[0]
#     for key in list(dict_averaged_count_dfs.keys()):
#         if dict_averaged_count_dfs[key]['average'][index] > dict_averaged_count_dfs[index_max]['average'][index]:
#             index_max = index
#     max_count_groups.append(index_max)
avg_count_total['Representative_group'] = df_max_counts.idxmax(axis=1).values

dict_averaged_score_dfs['total'] = avg_score_total
dict_averaged_count_dfs['total'] = avg_count_total

# avg_difference_score = avg_score_ctrl["average"] - avg_score_sick["average"]
# avg_difference_count = avg_count_ctrl["average"] - avg_count_sick["average"]


df_ref_significance_OTU = pd.read_csv(pipeline_path + '/Post-processing/outputs/'+dataset_name+'_OTU_rf.csv')

if args.tfigm:
    dataset_name = dataset_name + '_tfigm'
elif args.scaling:
    dataset_name = dataset_name + '_scaled'

for key in list(dict_averaged_score_dfs.keys()):
    dict_averaged_score_dfs[key].to_csv(pipeline_path + '/Post-processing/outputs/'+dataset_name+'_avg_difference_score_'+key+'.csv', sep = ',', index = True)
for key in list(dict_averaged_count_dfs.keys()):
    dict_averaged_count_dfs[key].to_csv(pipeline_path + '/Post-processing/outputs/'+dataset_name+'_avg_difference_count_'+key+'.csv', sep = ',', index = True)
