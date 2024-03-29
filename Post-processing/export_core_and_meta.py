import pandas as pd
import numpy as np
from collections import defaultdict
from Bio.ExPASy import Enzyme
from goatools import obo_parser
import requests
import os
import ast

import argparse

parser = argparse.ArgumentParser(description='A program that averages the RF importance scores of the functional annotations, and associates them to OTUs')

parser.add_argument("-d","--dataset_name", help="Name of the dataset used for this treatment")

parser.add_argument("-p", "--pipeline_path", help="Path to the whole pipeline folder")

parser.add_argument("-r", "--data_ref", help="Full reference to the dataset's folder")

parser.add_argument("-n", "--nb_repeats", help="Number of times the pipeline has been repeated")

parser.add_argument("-i", "--nb_iterations", help="Number of iterations per run of the pipeline")

parser.add_argument("--tfigm", action="store_true", help="Indicates whether the pipeline focuses on TF-IGM recalculated data")
parser.add_argument("--scaling", action="store_true", help="Indicates whether the pipeline focuses on scaled data")

args = parser.parse_args()

dataset_name = args.dataset_name

pipeline_path = args.pipeline_path

###FUNCTION

def extract_core_associates(dataframe, core_list):
    '''
    This function automatically extracts the core-significant from the lists of associated variables
    '''

    if 'Linked_Reactions' in dataframe.columns:
        col_ref = 'Linked_Reactions'
        new_col = 'Significant_linked_Reactions'
    else:
        col_ref = 'Linked_OTUs'
        new_col = 'Significant_linked_OTUs'

    signif_vars = []
    for vars_assoc in dataframe[col_ref].values:
        vars_list = []
        for var in ast.literal_eval(vars_assoc):
            if var in core_list:
                vars_list.append(var)
        signif_vars.append(vars_list)
    
    dataframe[new_col] = signif_vars
    return dataframe

    


###GATHER ALL SIGNIFICANT LISTS

link_to_test_folder = pipeline_path+'/Meta_Outputs/'+args.data_ref
iterations=int(args.nb_iterations)
dict_lists = {}
median_perfs_otus = defaultdict(list)
median_perfs_sofas = defaultdict(list)

for i in range(iterations):
    iter_folder = args.data_ref+'_'+str(i+1)
    for run_nb in os.listdir(link_to_test_folder+'/'+iter_folder+'/Selection_outputs'):
        otu_db = pd.read_csv(link_to_test_folder+'/'+iter_folder+'/Selection_outputs/'+run_nb+'/OTU_'+dataset_name+'.csv')
        score_db = pd.read_csv(link_to_test_folder+'/'+iter_folder+'/Selection_outputs/'+run_nb+'/scores_'+dataset_name+'.csv')
        
        signif_list_sofa = score_db.truncate(after = score_db.loc[score_db['ID']=='CUTOFF'].index.values[0]-1)['ID'].values
        signif_list_otu = otu_db.truncate(after = otu_db.loc[otu_db['ID']=='CUTOFF'].index.values[0]-1)['ID'].values

        if 'OTU' not in dict_lists.keys():
            dict_lists['OTU'] = {}
        
        if 'SoFA' not in dict_lists.keys():
            dict_lists['SoFA'] = {}
        
        if run_nb not in dict_lists['OTU'].keys():
            dict_lists['OTU'][run_nb]=defaultdict(list)
        if run_nb not in dict_lists['SoFA'].keys():
            dict_lists['SoFA'][run_nb]=defaultdict(list)
        
        dict_lists['OTU'][run_nb][i+1] = signif_list_otu
        dict_lists['SoFA'][run_nb][i+1] = signif_list_sofa

    for run_nb in os.listdir(link_to_test_folder+'/'+iter_folder+'/Classification_performances'):
        otu_classif_perfs = pd.read_csv(link_to_test_folder+'/'+iter_folder+'/Classification_performances/'+run_nb+'/OTU_classif_perfs.csv')['Test performance'].values
        sofa_classif_perfs = pd.read_csv(link_to_test_folder+'/'+iter_folder+'/Classification_performances/'+run_nb+'/SoFA_classif_perfs.csv')['Test performance'].values

        median_perfs_otus[run_nb].append(np.median(otu_classif_perfs))
        median_perfs_sofas[run_nb].append(np.median(sofa_classif_perfs))


for i in median_perfs_sofas.keys():
    median_perfs_sofas[i] = np.mean(median_perfs_sofas[i])

for i in median_perfs_otus.keys():
    median_perfs_otus[i] = np.mean(median_perfs_otus[i])

max_ind_sofas = max(median_perfs_sofas, key=median_perfs_sofas.get)
max_ind_otus = max(median_perfs_otus, key=median_perfs_otus.get)

###CREATE CORE AND META LISTS

if not os.path.exists(link_to_test_folder+'/core_and_meta_all_runs'):
    os.makedirs(link_to_test_folder+'/core_and_meta_all_runs')
if not os.path.exists(link_to_test_folder+'/core_and_meta_all_runs/All_iterations'):
    os.makedirs(link_to_test_folder+'/core_and_meta_all_runs/All_iterations')
if not os.path.exists(link_to_test_folder+'/core_and_meta_all_runs/Best_iteration'):
    os.makedirs(link_to_test_folder+'/core_and_meta_all_runs/Best_iteration')


core_meta_dfs_archive = {'OTU':{'core':defaultdict(defaultdict), 'meta':defaultdict(defaultdict)}, 'SoFA': {'core':defaultdict(defaultdict), 'meta':defaultdict(defaultdict)}}

for profile in ['OTU','SoFA']:
    for run in dict_lists[profile].keys():


        full_list = sum([list(dict_lists[profile][run][i]) for i in dict_lists[profile][run].keys()],[])
        core_meta = {'core':[], 'meta':{'ID':[],'Count':[]}}

        for annot in np.unique(full_list):
            if full_list.count(annot) == iterations:
                core_meta['core'].append(annot)
            else:
                core_meta['meta']['ID'].append(annot)
                core_meta['meta']['Count'].append(full_list.count(annot))

        
        core = pd.DataFrame(core_meta['core'], columns=['ID'])
        meta = pd.DataFrame(core_meta['meta']).sort_values(by='Count', ascending=False)

        # if profile =='SoFA':
        #     core = add_reaction_names(core)
        #     meta = add_reaction_names(meta)

        
        # core.to_csv(link_to_test_folder +'/core_and_meta_all_runs/core_'+profile+'_'+dataset_name+'_iteration_'+run.split('_')[-1]+'.csv', index=0)
        # meta.to_csv(link_to_test_folder +'/core_and_meta_all_runs/meta_'+profile+'_'+dataset_name+'_iteration_'+run.split('_')[-1]+'.csv', index=0)


        core_meta_dfs_archive[profile]['core'][run] = core
        core_meta_dfs_archive[profile]['meta'][run] = meta
        

##Getting the reference info for OTUs and annotations
otu_db = pd.read_csv(link_to_test_folder+'/'+args.data_ref+'_1/Selection_outputs/run_1/OTU_'+dataset_name+'.csv')[['Average_importance','ID','Linked_Reactions','Family']]
sofa_db = pd.read_csv(link_to_test_folder+'/'+args.data_ref+'_1/Selection_outputs/run_1/scores_'+dataset_name+'.csv')[['Average_importance','ID','Name','Linked_OTUs','Family']]

otu_db = otu_db[otu_db['ID']!='CUTOFF']
sofa_db = sofa_db[sofa_db['ID']!='CUTOFF']

##Adding OTU taxonomies
otu_names_ref = pd.read_csv(link_to_test_folder+'/'+args.data_ref+'_1/SoFAs/OTU_name_table_ref.tsv', sep='\t')
otu_taxonomies = []

for otu in otu_db['ID'].values:
    otu_taxo = otu_names_ref[otu_names_ref['observation_name']==otu]['taxonomic_affiliation'].values
    otu_taxonomies.append(otu_taxo)

otu_db['Taxonomy'] = otu_taxonomies

##Adding each variable's presence in each profile

labels = np.unique(pd.read_csv(pipeline_path+'/Inputs/Label_'+dataset_name+'.csv', header=None).values)
label_avg_counts ={'OTU':{}, 'SoFA':{}}

if 'tf_igm' in args.data_ref:
    data_call = dataset_name + '_tfigm'
elif 'scaled' in args.data_ref:
    data_call =args.data_ref
else:
    data_call = dataset_name

for lab in labels:
    otu_avg_label = pd.read_csv(pipeline_path+'/Post-processing/outputs/'+data_call+'_avg_difference_count_'+str(lab)+'.csv', index_col=0)['average']
    sofa_avg_label = pd.read_csv(pipeline_path+'/Post-processing/outputs/'+data_call+'_avg_difference_score_'+str(lab)+'.csv', index_col=0)['average']

    otu_db['average score in '+str(lab)] = otu_avg_label.loc[list(otu_db['ID'].values)].values
    sofa_db['average score in '+str(lab)] = sofa_avg_label.loc[list(sofa_db['ID'].values)].values
    


##Formatting the Core and Meta output files
for run in dict_lists[profile].keys():
    core_otu_df = core_meta_dfs_archive['OTU']['core'][run]
    meta_otu_df = core_meta_dfs_archive['OTU']['meta'][run]
    core_sofa_df = core_meta_dfs_archive['SoFA']['core'][run]
    meta_sofa_df = core_meta_dfs_archive['SoFA']['meta'][run]

    core_otu_df_plus_info = otu_db[otu_db['ID'].isin(core_otu_df['ID'])]
    core_sofa_df_plus_info = sofa_db[sofa_db['ID'].isin(core_sofa_df['ID'])]

    meta_otu_df_plus_info = otu_db[otu_db['ID'].isin(meta_otu_df['ID'])]
    meta_sofa_df_plus_info = sofa_db[sofa_db['ID'].isin(meta_sofa_df['ID'])]

    meta_otu_counts = []
    for otu_to_count in meta_otu_df_plus_info['ID'].values:
        otu_count = meta_otu_df[meta_otu_df['ID']==otu_to_count]['Count'].values[0]
        meta_otu_counts.append(otu_count)
    
    meta_sofa_counts = []
    for sofa_to_count in meta_sofa_df_plus_info['ID'].values:
        sofa_count = meta_sofa_df[meta_sofa_df['ID']==sofa_to_count]['Count'].values[0]
        meta_sofa_counts.append(sofa_count)

    meta_otu_df_plus_info['Count'] = meta_otu_counts
    meta_sofa_df_plus_info['Count'] = meta_sofa_counts


    core_otu_df_plus_info = extract_core_associates(core_otu_df_plus_info, core_sofa_df['ID'].values)
    meta_otu_df_plus_info = extract_core_associates(meta_otu_df_plus_info, core_sofa_df['ID'].values)
    
    core_sofa_df_plus_info = extract_core_associates(core_sofa_df_plus_info, core_otu_df['ID'].values)
    meta_sofa_df_plus_info = extract_core_associates(meta_sofa_df_plus_info, core_otu_df['ID'].values)

    core_sofa_df_plus_info = core_sofa_df_plus_info[['Average_importance','ID','Name','Linked_OTUs','Significant_linked_OTUs','Family']+['average score in '+str(lab) for lab in labels]]
    meta_sofa_df_plus_info = meta_sofa_df_plus_info[['ID','Count','Average_importance','Name','Linked_OTUs','Significant_linked_OTUs','Family']+['average score in '+str(lab) for lab in labels]].sort_values(by='Count', ascending=False)
    core_otu_df_plus_info = core_otu_df_plus_info[['Average_importance','ID','Taxonomy','Linked_Reactions','Significant_linked_Reactions','Family']+['average score in '+str(lab) for lab in labels]]
    meta_otu_df_plus_info = meta_otu_df_plus_info[['ID','Count','Average_importance','Taxonomy','Linked_Reactions','Significant_linked_Reactions','Family']+['average score in '+str(lab) for lab in labels]].sort_values(by='Count', ascending=False)

    core_sofa_df_plus_info.to_csv(link_to_test_folder +'/core_and_meta_all_runs/All_iterations/core_annots_'+dataset_name+'_iteration_'+run.split('_')[-1]+'.csv', index=0)
    meta_sofa_df_plus_info.to_csv(link_to_test_folder +'/core_and_meta_all_runs/All_iterations/meta_annots_'+dataset_name+'_iteration_'+run.split('_')[-1]+'.csv', index=0)
    core_otu_df_plus_info.to_csv(link_to_test_folder +'/core_and_meta_all_runs/All_iterations/core_OTU_'+dataset_name+'_iteration_'+run.split('_')[-1]+'.csv', index=0)
    meta_otu_df_plus_info.to_csv(link_to_test_folder +'/core_and_meta_all_runs/All_iterations/meta_OTU_'+dataset_name+'_iteration_'+run.split('_')[-1]+'.csv', index=0)

    if run == max_ind_sofas:
        core_sofa_df_plus_info.to_csv(link_to_test_folder +'/core_and_meta_all_runs/Best_iteration/core_annots_'+dataset_name+'_iteration_'+run.split('_')[-1]+'.csv', index=0)
        meta_sofa_df_plus_info.to_csv(link_to_test_folder +'/core_and_meta_all_runs/Best_iteration/meta_annots_'+dataset_name+'_iteration_'+run.split('_')[-1]+'.csv', index=0)

    if run == max_ind_otus:
        core_otu_df_plus_info.to_csv(link_to_test_folder +'/core_and_meta_all_runs/Best_iteration/core_OTU_'+dataset_name+'_iteration_'+run.split('_')[-1]+'.csv', index=0)
        meta_otu_df_plus_info.to_csv(link_to_test_folder +'/core_and_meta_all_runs/Best_iteration/meta_OTU_'+dataset_name+'_iteration_'+run.split('_')[-1]+'.csv', index=0)

    core_meta_dfs_archive['OTU']['core'][run] = core_otu_df_plus_info
    core_meta_dfs_archive['OTU']['meta'][run] = meta_otu_df_plus_info
    core_meta_dfs_archive['SoFA']['core'][run] = core_sofa_df_plus_info
    core_meta_dfs_archive['SoFA']['meta'][run] = meta_sofa_df_plus_info

    




