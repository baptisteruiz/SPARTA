import pandas as pd
import numpy as np
import ast
import os.path
import seaborn as sns
import matplotlib.pyplot as plt
from progress.bar import IncrementalBar
from collections import defaultdict

import argparse

parser = argparse.ArgumentParser(description='A program that averages the RF importance scores of the functional annotations, and associates them to OTUs')

parser.add_argument("-d","--dataset_name", help="Name of the dataset used for this treatment")

parser.add_argument("-p", "--pipeline_path", help="Path to the whole pipeline folder")

parser.add_argument("-r", "--data_ref", help="Full reference to the dataset's folder")

parser.add_argument("-i", "--iteration", help="Reference for the method's iteration")

parser.add_argument("--tfigm", action="store_true", help="Indicates whether the pipeline focuses on TF-IGM recalculated data")
parser.add_argument("--scaling", action="store_true", help="Indicates whether the pipeline focuses on scaled data")

args = parser.parse_args()

dataset_name = args.dataset_name
dataset_no_iter = args.dataset_name.split('_iteration')[0]

pipeline_path = args.pipeline_path

# dataset_name = "abundance_whatever_2_iteration_2"
# dataset_no_iter = dataset_name.split('_iteration')[0]
# pipeline_path = "/home/bruiz/Documents/Code/full_pipeline"
# iteration = 2
# data_ref = "abundance_whatever_2_tf_igm_relative"


def truncate_signif(df):
    random_index = df.index[df["ID"] == 'CUTOFF']
    new_df = df.iloc[:random_index[0]]

    return(new_df)

def extract_associated(df):
    list_relevant_otus_multi = df["Linked_OTUs"].values.tolist()
    
    list_relevant_otus = []
    for otulist in list_relevant_otus_multi:
        #print(type(otulist))
        otulist = ast.literal_eval(otulist)
        list_relevant_otus.extend(otulist)
    
    list_relevant_otus = list(set(list_relevant_otus))

    return list_relevant_otus

def extract_associated_reaction(df):
    list_relevant_reac_multi = df["Linked_Reactions"].values.tolist()
    
    list_relevant_reac = []
    for reaclist in list_relevant_reac_multi:
        #print(reaclist)
        reaclist = ast.literal_eval(reaclist)
        list_relevant_reac.extend(reaclist)
    
    list_relevant_reac = list(set(list_relevant_reac))

    return list_relevant_reac

def highlight_significant_otus(df, ref_list):

    significant_refs = []
    bar = IncrementalBar('Highlighting OTUs', max = len(df["Linked_OTUs"].values.tolist()))

    for otu_list in df["Linked_OTUs"]:
        found_significant = []
        if str(otu_list) not in ['nan', 'CUTOFF']:
            for otu in ast.literal_eval(otu_list):
                if otu in ref_list:
                    found_significant.append(otu)
        significant_refs.append(found_significant)

        bar.next()

    df["Significant_linked_OTUs"] = significant_refs 
    return df 

def highlight_significant_annots(df, ref_list):

    significant_refs = []
    bar = IncrementalBar('Highlighting Reactions', max = len(df["Linked_Reactions"].values.tolist()))

    for reac_list in df["Linked_Reactions"]:
        found_significant = []
        if str(reac_list) not in ['nan', 'CUTOFF']:
            for reac in ast.literal_eval(reac_list):
                if reac in ref_list:
                    found_significant.append(reac)
        significant_refs.append(found_significant)

        bar.next()

    df["Significant_linked_Reactions"] = significant_refs 
    return df 

def extract_id(df):
    list_relevant_id = df["ID"].values.tolist()
    #print(list_relevant_otus)
    return(list_relevant_id)

def compare_extract(list1, list2):
    list_not_found = []
    for item in list1:
        if item not in list2:
            list_not_found.append(item)
    
    percent = len(list_not_found)/len(list1) * 100

    return list_not_found, percent

def scatter_presentation(ax, x, y, palette):
### Takes charge of the aesthetics of the final scatter plot

    sns.set_context("paper")
    sns.set_style("darkgrid")
    z = y

    ax.scatter(x, y, c=z, cmap=palette, s=50)
    sns.despine()




os.chdir(os.path.dirname(__file__))

### Random Forest importance ranking for the OTU counts:
otu_table_full = pd.read_csv(pipeline_path + '/Post-processing/outputs/'+dataset_name+'_OTU_rf.csv')


if args.tfigm:
    dataset_name = dataset_name + '_tfigm'
    dataset_no_iter_ref = dataset_no_iter + '_tfigm'
elif args.scaling:
    dataset_name = dataset_name + '_scaled'
    dataset_no_iter_ref = dataset_no_iter + '_scaled'
else:
    dataset_no_iter_ref = dataset_no_iter

### Random Forest importance ranking for the annotation scores:
annotation_table_full = pd.read_csv(pipeline_path + '/Post-processing/outputs/'+dataset_name+'_annots_rf.csv')


#removing random
# annotation_table_full = annotation_table_full[annotation_table_full['ID'] != 'random']
# otu_table_full = otu_table_full[otu_table_full['ID'] != 'random']


# ### Only keep annotations that are associated with sick
# scores_difference_control_sick = pd.read_csv("/home/bruiz/Documents/Code/average_analysis_per_group/output/avg_difference_score_total.csv", sep =",", header = 0, index_col=0)
# annots_retain = []

# for ind_diff in range(len(scores_difference_control_sick["average"])):
#     if scores_difference_control_sick["average"][ind_diff] > 0: #<0: retain sick; >0: retain control
#         retain_annot = list(scores_difference_control_sick.index.values)[ind_diff]
#         annots_retain.append(retain_annot)

# annotation_table_full_sick = annotation_table_full[annotation_table_full["ID"].isin(annots_retain)]
# print("annots:", annotation_table_full_sick)
#Extract significant (above random) OTUs from OTU table, and annotations from the annotations table
annotation_table_significant = truncate_signif(annotation_table_full)
otu_table_significant = truncate_signif(otu_table_full)

#annotation_table_significant_sick = annotation_table_significant[annotation_table_significant["ID"].isin(annots_retain)]
#Re-format significant OTUs/annotations and OTUs/annotations associated with significant counterparts as lists
associated_OTU_list = extract_associated(annotation_table_significant)
significant_otu_list = extract_id(otu_table_significant)

associated_reaction_list = extract_associated_reaction(otu_table_significant)
significant_reaction_list = extract_id(annotation_table_significant)

#Highlight within the original tables which associated OTUs/annotations are significant

annotation_table_full_highlighted = highlight_significant_otus(annotation_table_full, significant_otu_list)
otu_table_full_highlighted = highlight_significant_annots(otu_table_full, significant_reaction_list)

### Original count/score files, after averaging steps

db_difference_sick_ctrl_scores = pd.read_csv(pipeline_path + '/Post-processing/outputs/'+dataset_no_iter_ref+'_avg_difference_score_total.csv', sep =",", header = 0,index_col=0)
db_difference_sick_ctrl_counts = pd.read_csv(pipeline_path + '/Post-processing/outputs/'+dataset_no_iter_ref+'_avg_difference_count_total.csv', sep =",", header = 0,index_col=0)

represented_groups = np.unique(np.concatenate([db_difference_sick_ctrl_scores['Representative_group'].values, db_difference_sick_ctrl_counts['Representative_group'].values]))
score_db_repres = {}
count_db_repres = {}

for group in represented_groups:

    score_db_group = pd.read_csv(pipeline_path + '/Post-processing/outputs/'+dataset_no_iter_ref+'_avg_difference_score_'+str(group)+'.csv', sep =",", header = 0,index_col=0)
    count_db_group = pd.read_csv(pipeline_path + '/Post-processing/outputs/'+dataset_no_iter_ref+'_avg_difference_count_'+str(group)+'.csv', sep =",", header = 0,index_col=0)

    score_db_repres[group] = score_db_group
    count_db_repres[group] = count_db_group


### Create lists of annotations that are more representative of one group or of another, for future reference

dict_rep_annots = defaultdict(list)


for annot in db_difference_sick_ctrl_scores.index.values.tolist():
    if annot !='CUTOFF':
        group = db_difference_sick_ctrl_scores['Representative_group'][annot]
        dict_rep_annots[group].append(annot)


# print("sick_ref_annot: ", list_sick_rep_annots)
# print("crtl_ref_annot: ", list_ctrl_rep_annots)

linked_otus_length = []
linked_reacs_length = []

linked_significant_otus_length = []
linked_significant_reacs_length = []

linked_reacs_opposite_length = []
linked_significant_reacs_opposite_length = []

fam_annot = []

### Re-calculate scores, and influence of significant OTUs in scores
bar = IncrementalBar('Looking up reaction scores', max = annotation_table_full_highlighted.shape[0] -1)

for line in range(annotation_table_full_highlighted.shape[0]):

    i = annotation_table_full_highlighted.iloc[line]['Linked_OTUs']
    i_signif = annotation_table_full_highlighted.iloc[line]['Significant_linked_OTUs']
    id_annot = annotation_table_full_highlighted.iloc[line]['ID']

    if id_annot == 'CUTOFF':
        fam_annot.append('NA')
        
    else:
        
        family = db_difference_sick_ctrl_scores['Representative_group'][id_annot]
        count_db = count_db_repres[family]
        fam_annot.append(family)

        if str(i) == 'nan':
            linked_otus_length.append(0)
            linked_significant_otus_length.append(0)
        else:
            score_otu = 0
            score_otu_signif = 0
            for otu in ast.literal_eval(i):
                path_to_OTU = pipeline_path + '/SoFA_calculation/outputs/'+dataset_no_iter+'/esmecata_outputs_annot/annotation_reference/' + otu + ".tsv"
                OTU_esmecata_annot = pd.read_csv(path_to_OTU, sep ="\t", header = 0)

                if "GO" in id_annot:
                    factor = OTU_esmecata_annot["GO"].astype(str).str.count(id_annot).sum()
                
                else:
                    factor = OTU_esmecata_annot["EC"].astype(str).str.count(id_annot).sum()

                score_otu += factor * count_db["average"][otu].astype(float)

                if otu in i_signif:
                    score_otu_signif += factor * count_db["average"][otu].astype(float)

            linked_otus_length.append(score_otu)
            linked_significant_otus_length.append(score_otu_signif)
    bar.next()

### Percentage of total linked reactions to OTU that are significant (takes account of group representation)

fam = []
list_same_fam = []
list_same_fam_signif = []
list_opposite_fam = []
list_opposite_fam_signif = []

bar = IncrementalBar('Looking up OTU scores', max = otu_table_full_highlighted.shape[0] - 1)

for i in list(otu_table_full_highlighted[otu_table_full_highlighted["ID"] != 'CUTOFF']["ID"].index):
    
    otu = otu_table_full_highlighted["ID"][i]

    if otu == 'CUTOFF':
        linked_reacs_length.append(0)
        linked_significant_reacs_length.append(0)
        fam.append('NA')
        list_same_fam.append('NA')
        list_same_fam_signif.append('NA')
        list_opposite_fam.append('NA')
        list_opposite_fam_signif.append('NA')


    else:
        
        score_otu = 0
        score_otu_signif = 0

        score_otu_opposite = 0
        score_otu_signif_opposite = 0
        
        list_same_fam_i = []
        list_same_fam_signif_i = []
        list_opposite_fam_i = []
        list_opposite_fam_signif_i = []

        path_to_OTU = pipeline_path + '/SoFA_calculation/outputs/'+dataset_no_iter+'/esmecata_outputs_annot/annotation_reference/' + otu + ".tsv"

        if os.path.exists(path_to_OTU):

            OTU_esmecata_annot = pd.read_csv(path_to_OTU, sep ="\t", header = 0)

            family = db_difference_sick_ctrl_counts['Representative_group'][otu]
            score_ref_local = dict_rep_annots[family]
            fam.append(family)

            for reac in ast.literal_eval(otu_table_full_highlighted['Linked_Reactions'][i]):



                if "GO" in reac:
                    factor = OTU_esmecata_annot["GO"].str.count(reac).sum()
                
                else:
                    factor = OTU_esmecata_annot["EC"].str.count(reac).sum()
            
                if reac in score_ref_local:
                    list_same_fam_i.append(reac)
                    score_otu += factor
                    if reac in otu_table_full_highlighted['Significant_linked_Reactions'][i]:
                        score_otu_signif += factor
                        list_same_fam_signif_i.append(reac)


                else:
                    list_opposite_fam_i.append(reac)
                    score_otu_opposite += factor
                    if reac in otu_table_full_highlighted['Significant_linked_Reactions'][i]:
                        score_otu_signif_opposite += factor
                        list_opposite_fam_signif_i.append(reac)

        else:
            fam.append('NA')

            
        linked_reacs_length.append(score_otu)
        linked_significant_reacs_length.append(score_otu_signif)

        linked_reacs_opposite_length.append(score_otu_opposite)
        linked_significant_reacs_opposite_length.append(score_otu_signif_opposite)

        list_same_fam.append(list_same_fam_i)
        list_same_fam_signif.append(list_same_fam_signif_i)
        list_opposite_fam.append(list_opposite_fam_i)
        list_opposite_fam_signif.append(list_opposite_fam_signif_i)

    
    bar.next()

# for i in otu_table_full_highlighted['Significant_linked_Reactions']:
#     if str(i) == 'nan':
#         linked_significant_reacs_length.append(0)
#     else:
#         linked_significant_reacs_length.append(len(i))

signif_otu_percentage = []
signif_reac_percentage = []

for k in range(len(linked_otus_length)):
    if linked_otus_length[k] != 0:
        signif_otu_percentage.append(linked_significant_otus_length[k] / linked_otus_length[k])
    else:
        signif_otu_percentage.append(0)

# aberrant_index = []
# for j in range(len(signif_otu_percentage)):
#     if signif_otu_percentage[j] > 1:
#         signif_otu_percentage[j] = -1
#         aberrant_index.append(j)


for k in range(len(linked_reacs_length)):
    if linked_reacs_length[k] != 0:
        signif_reac_percentage.append(linked_significant_reacs_length[k] / linked_reacs_length[k])
    else:
        signif_reac_percentage.append(0)


### Plotting results

cutoff_annots_index = annotation_table_full_highlighted.index[annotation_table_full_highlighted["ID"] == 'CUTOFF']
cutoff_otus_index = otu_table_full_highlighted.index[otu_table_full_highlighted["ID"] == 'CUTOFF']


fig, axs = plt.subplots(2, 1, constrained_layout=True)
fig.set_size_inches(20,14)

# axs[0].scatter(x=[i for i in range(len(annotation_table_full_highlighted['Significant_linked_OTUs']) -1)], y=signif_otu_percentage, c='r')
# axs[0].set_ylabel('Impact_of_Significant_Linked_OTUs_on_score')

# axs[1].scatter(x=[i for i in range(len(otu_table_full_highlighted['Significant_linked_Reactions']) -1)], y=signif_reac_percentage, c='b')
# axs[1].set_ylabel('Percentage_of_Significant_Linked_Reactions')

scatter_presentation(axs[0], x=[i for i in range(len(signif_otu_percentage))], y=signif_otu_percentage, palette = sns.dark_palette("red", as_cmap=True))
axs[0].axvline(x=cutoff_annots_index[0] - 0.5, color='r', ls= '--', label='Cutoff')
axs[0].text(cutoff_annots_index[0] - 0.5, max(signif_otu_percentage) + 0.05, "Cutoff (x = "+str(cutoff_annots_index[0])+")", color='r', weight = 'bold', size = 10)
axs[0].set_ylabel('Impact of Significant Linked OTUs on score', size = 12, weight = 'bold')
axs[0].set_xlabel('Functional annotations ranked by feature importance', size = 12, weight = 'bold')

scatter_presentation(axs[1], x=[i for i in range(len(signif_reac_percentage))], y=signif_reac_percentage, palette = sns.dark_palette("#69d", as_cmap=True))
axs[1].axvline(x=cutoff_otus_index[0] - 0.5, color='b', ls= '--', label='Cutoff')
axs[1].text(cutoff_otus_index[0] - 0.5, max(signif_reac_percentage) + 0.05, "Cutoff (x = "+str(cutoff_otus_index[0])+ ")", color='b', weight = 'bold', size = 10)
axs[1].set_ylabel('Percentage of Significant Linked Reactions', size = 12, weight = 'bold')
axs[1].set_xlabel('OTUs ranked by feature importance', size = 12, weight = 'bold')


plt.savefig(pipeline_path+'/Outputs/'+ args.data_ref + '/Selection_outputs/run_'+str(args.iteration)+'/' +dataset_name + '_correlation_plots.png', dpi = 100)


signif_otu_percentage.insert(len(annotation_table_significant.index), 'CUTOFF')
signif_reac_percentage.insert(len(otu_table_significant.index), 'CUTOFF')

fam.insert(len(otu_table_significant.index), 'NA')
list_same_fam.insert(len(otu_table_significant.index), 'NA')
list_same_fam_signif.insert(len(otu_table_significant.index), 'NA')
list_opposite_fam.insert(len(otu_table_significant.index), 'NA')
list_opposite_fam_signif.insert(len(otu_table_significant.index), 'NA')


linked_reacs_length.insert(len(otu_table_significant.index), 'NA')
linked_significant_reacs_length.insert(len(otu_table_significant.index), 'NA')
linked_reacs_opposite_length.insert(len(otu_table_significant.index), 'NA')
linked_significant_reacs_opposite_length.insert(len(otu_table_significant.index), 'NA')

annotation_table_full_highlighted['signif_otu_percentage'] = signif_otu_percentage
otu_table_full_highlighted['signif_reac_percentage'] = signif_reac_percentage

print(fam)
otu_table_full_highlighted['Family'] = fam

otu_table_full_highlighted['Annotations_from_same_Family'] = list_same_fam
otu_table_full_highlighted['Score_Annotations_from_same_Family'] = linked_reacs_length

otu_table_full_highlighted['Significant_Annotations_from_same_Family'] = list_same_fam_signif
otu_table_full_highlighted['Score_Significant_Annotations_from_same_Family'] = linked_significant_reacs_length

otu_table_full_highlighted['Annotations_from_opposite_Family'] = list_opposite_fam
otu_table_full_highlighted['Score_Annotations_from_opposite_Family'] = linked_reacs_opposite_length


otu_table_full_highlighted['Significant_Annotations_from_opposite_Family'] = list_opposite_fam_signif
otu_table_full_highlighted['Score_Significant_Annotations_from_opposite_Family'] = linked_significant_reacs_opposite_length


annotation_table_full_highlighted['Family'] = fam_annot


annotation_table_full_highlighted.to_csv(pipeline_path + "/Outputs/"+ args.data_ref + "/Selection_outputs/run_" +str(args.iteration)+ "/scores_"+dataset_no_iter+".csv", index = True)
otu_table_full_highlighted.to_csv(pipeline_path + "/Outputs/"+ args.data_ref + "/Selection_outputs/run_"+str(args.iteration)+"/OTU_"+dataset_no_iter+".csv", index = True)

#Compare the lists of OTUs/annotations obtained with each method
associated_not_in_otu_significant, non_assoc_percent_associated_otu = compare_extract(associated_OTU_list,significant_otu_list)
otu_significant_not_in_associated, non_assoc_percent_otu_significant = compare_extract(significant_otu_list,associated_OTU_list)

associated_not_in_reac_significant, non_assoc_percent_associated_reac = compare_extract(associated_reaction_list,significant_reaction_list)
reac_significant_not_in_associated, non_assoc_percent_reac_significant = compare_extract(significant_reaction_list,associated_reaction_list)


with open(pipeline_path + '/Outputs/'+ args.data_ref + '/Selection_outputs/run_'+str(args.iteration)+'/'+args.dataset_name+'_output_comparison_signif_rf.txt', 'w') as f:
    f.writelines("--------------------- OTU comparison analysis ------------------- \n")
    f.writelines("Percentage of non asociated OTUs in the significant list: " + str(non_assoc_percent_otu_significant) + "% \n")
    f.writelines("Detail [" + str(len(otu_significant_not_in_associated)) + "]:\n")
    f.writelines(str(otu_significant_not_in_associated))
    f.writelines('\n \n')
    f.writelines("Percentage of non significant OTUs in the associated list: " + str(non_assoc_percent_associated_otu) + "% \n")
    f.writelines("Detail [" + str(len(associated_not_in_otu_significant)) + "]:\n")
    f.writelines(str(associated_not_in_otu_significant))
    f.writelines('\n \n')
    f.writelines("--------------------- Functional annotations comparison analysis ------------------- \n")
    f.writelines("Percentage of non asociated annotations in the significant list: " + str(non_assoc_percent_reac_significant) + "% \n")
    f.writelines("Detail [" + str(len(reac_significant_not_in_associated)) + "]:\n")
    f.writelines(str(reac_significant_not_in_associated))
    f.writelines('\n \n')
    f.writelines("Percentage of non significant reactions in the associated list: " + str(non_assoc_percent_associated_reac) + "% \n")
    f.writelines("Detail [" + str(len(associated_not_in_reac_significant)) + "]:\n")
    f.writelines(str(associated_not_in_reac_significant))