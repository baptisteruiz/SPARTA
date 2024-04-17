import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
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


###COMPILING INFO ABOUT ALL RUNS AND GATHERING ALL SIGNIFICANT LISTS

path_save = pipeline_path + '/Meta_Outputs/'+args.data_ref

median_perf_values_dict = {}

median_otu_test_perfs = defaultdict(list)
median_sofa_test_perfs = defaultdict(list)


dict_lists = {}

for run_nb in range(1, int(args.nb_repeats) + 1):
     
    path_ref = path_save +'/'+args.data_ref+'_'+str(run_nb)

    for i in range(int(args.nb_iterations)):
    #Gather test classification performances

        local_path = path_ref + '/Classification_performances/run_'+str(i+1)+'/SoFA_classif_perfs.csv'   
        print(local_path)

        dataframe_sofa = pd.read_csv(local_path, sep = ',')
        median_sofa_test_perfs[i].append(np.median(list(dataframe_sofa['Test performance'].values)))

        local_path_otu = path_ref + '/Classification_performances/run_'+str(i+1)+'/OTU_classif_perfs.csv'
        dataframe_otu = pd.read_csv(local_path_otu, sep = ',')
        median_otu_test_perfs[i].append(np.median(list(dataframe_otu['Test performance'].values)))
    
    median_perf_values_dict[dataset_name] = [median_sofa_test_perfs, median_otu_test_perfs]

    #for run_folder in os.listdir(path_ref+'/Selection_outputs'):
    for i in range(int(args.nb_iterations)):
        run_folder = 'run_'+str(i+1)
    #Gather the lists of selected variables per run and iteration
        otu_db = pd.read_csv(path_ref+'/Selection_outputs/'+run_folder+'/OTU_'+dataset_name+'.csv')
        score_db = pd.read_csv(path_ref+'/Selection_outputs/'+run_folder+'/scores_'+dataset_name+'.csv')
        
        signif_list_sofa = score_db.truncate(after = score_db.loc[score_db['ID']=='CUTOFF'].index.values[0]-1)['ID'].values
        signif_list_otu = otu_db.truncate(after = otu_db.loc[otu_db['ID']=='CUTOFF'].index.values[0]-1)['ID'].values

        if 'OTU' not in dict_lists.keys():
            dict_lists['OTU'] = {}
        
        if 'SoFA' not in dict_lists.keys():
            dict_lists['SoFA'] = {}
        
        if run_folder not in dict_lists['OTU'].keys():
            dict_lists['OTU'][run_folder]=defaultdict(list)
        if run_folder not in dict_lists['SoFA'].keys():
            dict_lists['SoFA'][run_folder]=defaultdict(list)
        
        dict_lists['OTU'][run_folder][run_nb] = signif_list_otu
        dict_lists['SoFA'][run_folder][run_nb] = signif_list_sofa
        
    



###PLOTTING MEDIAN PERFORMANCES

ticks = ['Original Data', 'First Selection', 'Second Selection', 'Third Selection', 'Fourth Selection']

sns.set_theme()
sns.set_context('paper')
sns.set_style('whitegrid')

fig,ax = plt.subplots(1,constrained_layout=True, figsize=(10,10))

best_runs = []
max_vals = []
for i, disease_name in enumerate(median_perf_values_dict.keys()):
    
    
    sofa_test_perfs = median_perf_values_dict[disease_name][0]

    otu_test_perfs = median_perf_values_dict[disease_name][1]
    
    median_sofas_perfs_tfigm = [np.mean(sofa_test_perfs[i]) for i in sofa_test_perfs.keys()]
    median_otu_perfs = [np.mean(otu_test_perfs[i]) for i in otu_test_perfs.keys()]
    

    max_ind_sofas = median_sofas_perfs_tfigm.index(max(median_sofas_perfs_tfigm))
    max_ind_otus = median_otu_perfs.index(max(median_otu_perfs))


    best_runs+=[max_ind_sofas] + [max_ind_otus]

    max_vals.append(max(sofa_test_perfs[max_ind_sofas]))
    max_vals.append(max(otu_test_perfs[max_ind_otus]))

    dataframe = pd.DataFrame.from_dict({'AUC score': sofa_test_perfs[max_ind_sofas]+otu_test_perfs[max_ind_otus]})
    dataframe['Profile']=['SoFAs' for i in range(len(sofa_test_perfs[max_ind_sofas]))]+['OTUs' for i in range(len(otu_test_perfs[max_ind_otus]))]
    dataframe['Disease']=[disease_name.split('_')[1] for i in range(len(sofa_test_perfs[max_ind_sofas])+len(otu_test_perfs[max_ind_otus]))]
    dataframe['Best run']=[max_ind_sofas for i in range(len(sofa_test_perfs[max_ind_sofas]))] + [max_ind_otus for i in range(len(otu_test_perfs[max_ind_otus]))]
    dataframe['Median']=[median_sofas_perfs_tfigm[max_ind_sofas] for i in range(len(sofa_test_perfs[max_ind_sofas]))] + [median_otu_perfs[max_ind_otus] for i in range(len(otu_test_perfs[max_ind_otus]))]
    
    if i == 0:
        df_cumul = dataframe
    else:
        df_cumul = pd.concat([df_cumul, dataframe])

medians = sns.boxplot(
            y='Disease',
            x='AUC score',
            hue = 'Profile',
            data=df_cumul,
            ax=ax,
            palette = 'pastel',
            fliersize = 0)
sns.stripplot(y='Disease', x='AUC score', hue = 'Profile', data=df_cumul, jitter = True, dodge=True, ax=ax, palette = 'dark')


max_val = max(max_vals)
plt.text(max_val + 0.01, -0.5, 'Optimal number of selections', size=15, fontweight = 'bold')
for xi,yi,tx, col in zip([max_val + 0.03 for i in max_vals], (np.arange(len(dataframe['Median'].values))/2)-0.2, best_runs, ['blue','darkorange']*6):
     plt.text( xi, yi,tx, size=18, fontweight = 'bold', c=col)
# plt.text(1.01, -0.8, 'Optimal selection', size=10, fontweight = 'bold')

handles, labels = ax.get_legend_handles_labels()
l = plt.legend(handles[2:], labels[2:], loc='upper left', borderaxespad=0., fontsize = 15)

plt.yticks(size=15)
ax.set_ylabel('Disease',size=20)

plt.xticks(size=15)
ax.set_xlabel('Median AUC',size=20)

plt.savefig(path_save +'/median_OTU_vs_SoFA_(best_vs_best).png', dpi=300)


###CREATE CORE AND META LISTS

if not os.path.exists(path_save+'/core_and_meta_all_runs'):
    os.makedirs(path_save+'/core_and_meta_all_runs')
if not os.path.exists(path_save+'/core_and_meta_all_runs/All_iterations'):
    os.makedirs(path_save+'/core_and_meta_all_runs/All_iterations')
if not os.path.exists(path_save+'/core_and_meta_all_runs/Best_iteration'):
    os.makedirs(path_save+'/core_and_meta_all_runs/Best_iteration')


core_meta_dfs_archive = {'OTU':{'core':defaultdict(defaultdict), 'meta':defaultdict(defaultdict)}, 'SoFA': {'core':defaultdict(defaultdict), 'meta':defaultdict(defaultdict)}}

for profile in ['OTU','SoFA']:
    for run in dict_lists[profile].keys():

        full_list = sum([list(dict_lists[profile][run][i]) for i in dict_lists[profile][run].keys()],[])
        core_meta = {'core':[], 'meta':{'ID':[],'Count':[]}}

        for annot in np.unique(full_list):
            if full_list.count(annot) == int(args.nb_repeats):
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
otu_db = pd.read_csv(path_save+'/'+args.data_ref+'_1/Selection_outputs/run_1/OTU_'+dataset_name+'.csv')[['ID','Linked_Reactions','Family']]
sofa_db = pd.read_csv(path_save+'/'+args.data_ref+'_1/Selection_outputs/run_1/scores_'+dataset_name+'.csv')[['ID','Name','Linked_OTUs','Family']]

otu_db = otu_db[otu_db['ID']!='CUTOFF']
sofa_db = sofa_db[sofa_db['ID']!='CUTOFF']

##Adding OTU taxonomies
otu_names_ref = pd.read_csv(path_save+'/'+args.data_ref+'_1/SoFAs/OTU_name_table_ref.tsv', sep='\t')
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

    core_sofa_df_plus_info = core_sofa_df_plus_info[['ID','Name','Linked_OTUs','Significant_linked_OTUs','Family']+['average score in '+str(lab) for lab in labels]]
    meta_sofa_df_plus_info = meta_sofa_df_plus_info[['ID','Count','Name','Linked_OTUs','Significant_linked_OTUs','Family']+['average score in '+str(lab) for lab in labels]].sort_values(by='Count', ascending=False)
    core_otu_df_plus_info = core_otu_df_plus_info[['ID','Taxonomy','Linked_Reactions','Significant_linked_Reactions','Family']+['average score in '+str(lab) for lab in labels]]
    meta_otu_df_plus_info = meta_otu_df_plus_info[['ID','Count','Taxonomy','Linked_Reactions','Significant_linked_Reactions','Family']+['average score in '+str(lab) for lab in labels]].sort_values(by='Count', ascending=False)

    core_sofa_df_plus_info.to_csv(path_save +'/core_and_meta_all_runs/All_iterations/core_annots_'+dataset_name+'_iteration_'+run.split('_')[-1]+'.csv', index=0)
    meta_sofa_df_plus_info.to_csv(path_save +'/core_and_meta_all_runs/All_iterations/meta_annots_'+dataset_name+'_iteration_'+run.split('_')[-1]+'.csv', index=0)
    core_otu_df_plus_info.to_csv(path_save +'/core_and_meta_all_runs/All_iterations/core_OTU_'+dataset_name+'_iteration_'+run.split('_')[-1]+'.csv', index=0)
    meta_otu_df_plus_info.to_csv(path_save +'/core_and_meta_all_runs/All_iterations/meta_OTU_'+dataset_name+'_iteration_'+run.split('_')[-1]+'.csv', index=0)
    
    if run == 'run_'+str(max_ind_sofas):
        core_sofa_opti_df = core_sofa_df_plus_info
        meta_sofa_opti_df= meta_sofa_df_plus_info

    if run == 'run_'+str(max_ind_otus):
        core_otu_opti_df = core_otu_df_plus_info
        meta_otu_opti_df = meta_otu_df_plus_info

if max_ind_sofas == 0:
    core_sofa_opti_df = sofa_db

if max_ind_otus == 0:
    core_otu_opti_df = otu_db

opti_core_annots = core_sofa_opti_df['ID'].values
opti_core_otus = core_otu_opti_df['ID'].values

core_otu_opti_df = extract_core_associates(core_otu_opti_df, opti_core_annots)
core_sofa_opti_df = extract_core_associates(core_sofa_opti_df, opti_core_otus)

if max_ind_otus > 0:
    meta_otu_opti_df = extract_core_associates(meta_otu_opti_df, opti_core_annots)
    meta_otu_opti_df.to_csv(path_save +'/core_and_meta_all_runs/Best_iteration/meta_OTU_'+dataset_name+'_iteration_'+str(max_ind_otus)+'.csv', index=0)

if max_ind_sofas > 0:
    meta_sofa_opti_df = extract_core_associates(meta_sofa_opti_df, opti_core_otus)
    meta_sofa_opti_df.to_csv(path_save +'/core_and_meta_all_runs/Best_iteration/meta_annots_'+dataset_name+'_iteration_'+str(max_ind_sofas)+'.csv', index=0)

    
core_sofa_opti_df.to_csv(path_save +'/core_and_meta_all_runs/Best_iteration/core_annots_'+dataset_name+'_iteration_'+str(max_ind_sofas)+'.csv', index=0)
core_otu_opti_df.to_csv(path_save +'/core_and_meta_all_runs/Best_iteration/core_OTU_'+dataset_name+'_iteration_'+str(max_ind_otus)+'.csv', index=0)
