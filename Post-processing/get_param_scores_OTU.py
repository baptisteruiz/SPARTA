import pandas as pd
from progress.bar import IncrementalBar
import os.path
from operator import methodcaller

import argparse

parser = argparse.ArgumentParser(description='A program that averages the RF importance scores of the OTUs, and associates them to functional annotations')

parser.add_argument("-d","--dataset_name", help="Name of the dataset used for this treatment")

parser.add_argument("-p", "--pipeline_path", help="Path to the whole pipeline folder")

parser.add_argument("--tfigm", action="store_true", help="Indicates whether the pipeline focuses on TF-IGM recalculated data")

args = parser.parse_args()

dataset_name = args.dataset_name
dataset_no_iter = args.dataset_name.split('_iteration')[0]
pipeline_path = args.pipeline_path

# dataset_name = "abundance_whatever_2_iteration_2"
# dataset_no_iter = dataset_name.split('_iteration')[0]
# pipeline_path = "/home/bruiz/Documents/Code/full_pipeline"

def df_average(dataframe):
    nb_of_parameters = max(dataframe.iloc[:,0])

    bar = IncrementalBar('Averaging', max = nb_of_parameters + 1)

    averaged_data = []

    for i in range(0,nb_of_parameters + 1):
        subset = dataframe[dataframe.iloc[:,0] == i]
        if i == nb_of_parameters:
            avg_subset = subset.quantile(q=0.95, axis=0)
        else:
            avg_subset = subset.mean(axis=0)
        avg_importance_i = avg_subset[1]

        averaged_data.append([i,avg_importance_i])

        bar.next()
    
    averaged_data = pd.DataFrame(averaged_data)
    return(averaged_data)

def get_otu_codes(dataframe, id_dataframe):
    list_of_id = id_dataframe.iloc[:,0]
    #list_of_id.loc[list_of_id.shape[0]] = 'random'
    dataframe[len(dataframe.columns)] = list_of_id
    dataframe = dataframe.drop(0 , axis = 1)
    

    dataframe = dataframe.sort_values(1, axis = 0, ascending = False)

    dataframe.columns = ["Average_importance","ID"]

    return dataframe

def find_relevant_reactions(dataframe, path):

    bar = IncrementalBar('Checking out reactions', max = len([f for f in os.listdir(path) if not f.startswith('.')]))

    found_reac = {}
    for filename in [f for f in os.listdir(path) if not f.startswith('.')]:
        #print(filename)
        otu_data = pd.read_csv(path + "/" + filename, sep ='\t', keep_default_na=False)
        otu_name = filename.replace(".tsv","")
        go_list_otu = otu_data["GO"].values.tolist()
        ec_list_otu = otu_data["EC"].values.tolist()
        annot_list_otu = go_list_otu + ec_list_otu
        
        annot_list_otu = list(filter(('').__ne__, annot_list_otu))

        annot_list_otu = list(map(methodcaller("split", ","), annot_list_otu))

        annot_list_otu_flattened = []

        for annotlist in annot_list_otu:
            annot_list_otu_flattened.extend(annotlist)

        annot_list_otu_flattened = list(set(annot_list_otu_flattened))

        found_reac[otu_name] = annot_list_otu_flattened
        bar.next()

    dataframe["Linked_Reactions"] = dataframe["ID"].map(found_reac)
    return dataframe

# def add_reaction_names(dataframe):

#     if not(os.path.isfile('filename.txt')):
#         url = "https://ftp.expasy.org/databases/enzyme/enzyme.dat"
#         r = requests.get(url, allow_redirects=True)
#         open('enzyme.dat', 'wb').write(r.content)

#     bar = IncrementalBar('Recovering Names', max = len(dataframe.iloc[:,1]))

#     reaction_names = []
#     for id in dataframe.iloc[:,1]:
#         #print("id=", id)
#         handle = open("/home/bruiz/Documents/Code/Sorties_classement_scores_rf/enzyme.dat")
#         records = Enzyme.parse(handle)
#         reaction_names.append("Not Found")
#         for record in records:
#             if record["ID"] == id:
#                 reaction_names[-1] = record["DE"]
#                 #print("record_id = ", record["DE"])
#         bar.next()
#     #print(len(reaction_names))
#     dataframe[len(dataframe.columns)+1] = reaction_names
#     return dataframe


data_score_rf = pd.read_csv(pipeline_path + '/DeepMicro/results/entree_DeepMicro_'+dataset_name+'_OTU_best_features_random_rf.txt', header = None)
data_labels = pd.read_csv(pipeline_path+'/SoFA_calculation/outputs/'+dataset_no_iter+'/'+dataset_name+'_stripped.tsv', sep ='\t')
path_to_annotation_ref = pipeline_path + '/SoFA_calculation/outputs/'+dataset_no_iter+'/esmecata_outputs_annot/annotation_reference'

data_score_rf_average = df_average(data_score_rf)
sorted_score_codes = get_otu_codes(data_score_rf_average, data_labels)
print()
data_score_codes_plus_reactions = find_relevant_reactions(sorted_score_codes, path_to_annotation_ref)

# data_score_codes_names = add_reaction_names(sorted_score_codes)

data_score_codes_plus_reactions.to_csv(pipeline_path + '/Post-processing/outputs/'+dataset_name+'_OTU_rf.csv', index = None)

    
