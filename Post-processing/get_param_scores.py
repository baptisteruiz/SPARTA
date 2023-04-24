import pandas as pd
from progress.bar import IncrementalBar
import os.path
import requests
from Bio.ExPASy import Enzyme
from goatools import obo_parser

import argparse

parser = argparse.ArgumentParser(description='A program that averages the RF importance scores of the functional annotations, and associates them to OTUs')

parser.add_argument("-d","--dataset_name", help="Name of the dataset used for this treatment")

parser.add_argument("-p", "--pipeline_path", help="Path to the whole pipeline folder")

parser.add_argument("--tfigm", action="store_true", help="Indicates whether the pipeline focuses on TF-IGM recalculated data")
parser.add_argument("--scaling", action="store_true", help="Indicates whether the pipeline focuses on scaled data")


args = parser.parse_args()

dataset_name = args.dataset_name
dataset_no_iter = args.dataset_name.split('_iteration')[0]
pipeline_path = args.pipeline_path

# dataset_name = "abundance_whatever_2_iteration_2"
# dataset_no_iter = dataset_name.split('_iteration')[0]
# pipeline_path = "/home/bruiz/Documents/Code/full_pipeline"


def df_average(dataframe):
    '''
    This function calculates the average importance of each parameter in the RF classifications
    '''
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

def get_metabolic_codes(dataframe, id_dataframe):
    '''
    This function associates the parameter IDs (GO terms or EC numbers) to their scores and orders them by decreasing values
    '''
    list_of_id = id_dataframe.iloc[:,0]
    list_of_id = [id.replace('_', '.') for id in list_of_id]
    #list_of_id.append("random")
    dataframe[len(dataframe.columns)] = list_of_id
    dataframe = dataframe.drop(0 , axis = 1)
    

    dataframe = dataframe.sort_values(1, axis = 0, ascending = False)
    return dataframe

def add_reaction_names(dataframe):
    '''
    This function queries the OBO and ExPASy data banks to get the names of the IDs
    '''
    
    bar = IncrementalBar('Recovering Names', max = len(dataframe.iloc[:,1]))

    data_folder = pipeline_path + '/Post-processing/data'

    # Check if we have the ./data directory already
    if(not os.path.isfile(data_folder)):
        # Emulate mkdir -p (no error if folder exists)
        try:
            os.mkdir(data_folder)
        except OSError as e:
            if(e.errno != 17):
                raise e
    else:
        raise Exception('Data path (' + data_folder + ') exists as a file. '
                    'Please rename, remove or change the desired location of the data path.')
    
    #Check if we have enzyme file, if not download it
    if not(os.path.isfile(data_folder + '/enzyme.dat')):
        url = "https://ftp.expasy.org/databases/enzyme/enzyme.dat"
        r = requests.get(url, allow_redirects=True)
        open(data_folder + '/enzyme.dat', 'wb').write(r.content)

    #Same with GO file
    if(not os.path.isfile(data_folder+'/go-basic.obo')):
        go_obo_url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
        r = requests.get(go_obo_url, allow_redirects=True)
        open(data_folder + '/go-basic.obo', 'wb').write(r.content)

    
    go = obo_parser.GODag(data_folder+'/go-basic.obo')
    
    reaction_names = []

    for id in dataframe.iloc[:,1]:
        #print("id=", id)
        reaction_names.append("Not Found")

        if "GO" in id:
            if id in go:
                go_term = go[id]
                reaction_names[-1] = go_term.name
        else:
            handle = open(data_folder + '/enzyme.dat')
            records = Enzyme.parse(handle)
            de_found = next((item["DE"] for item in records if item["ID"] == id), "Not Found")
            reaction_names[-1] = de_found
            #print("record_id = ", record["DE"])
        
        bar.next()
    #print(len(reaction_names))
    dataframe[len(dataframe.columns)+1] = reaction_names
    dataframe.columns = ["Average_importance","ID","Name"]
    return dataframe

def find_relevant_otus(dataframe, path):

    bar = IncrementalBar('Checking out OTUs', max = len([f for f in os.listdir(path) if not f.startswith('.')]))

    found_otu = {}
    for filename in [f for f in os.listdir(path) if not f.startswith('.')]:
        #print(filename)
        otu_data = pd.read_csv(path + "/" + filename, sep ='\t', keep_default_na=False)
        otu_name = filename.replace(".tsv","")
        for id in dataframe["ID"]:
            
            go = False
            ec = False

            for golist in otu_data["GO"].values:
                
                if (id in golist) and (not pd.isna(golist)):
                    go = True
            
            for eclist in otu_data["EC"].values:
            
                if (id in eclist) and (not pd.isna(eclist)):
                    ec = True
            
            if (go) | (ec):
                if id not in found_otu.keys():
                    found_otu[id] = []
                found_otu[id].append(otu_name)
        bar.next()

    dataframe["Linked_OTUs"] = dataframe["ID"].map(found_otu)
    return dataframe

#os.chdir(os.path.dirname(__file__))

if args.tfigm:
    data_score_rf = pd.read_csv(pipeline_path + '/DeepMicro/results/entree_DeepMicro_'+dataset_name+'_tfigm_best_features_random_rf.txt', header = None)
elif args.scaling:
    data_score_rf = pd.read_csv(pipeline_path + '/DeepMicro/results/entree_DeepMicro_'+dataset_name+'_scaled_best_features_random_rf.txt', header = None)
else:
    data_score_rf = pd.read_csv(pipeline_path + '/DeepMicro/results/entree_DeepMicro_'+dataset_name+'_best_features_random_rf.txt', header = None)

data_labels = pd.read_csv(pipeline_path + '/SoFA_calculation/outputs/'+dataset_no_iter+'/score_annot_selon_otus_'+dataset_name+'.tsv', header = 0, sep ='\t')
path_to_annotation_ref = pipeline_path + '/SoFA_calculation/outputs/'+dataset_no_iter+'/esmecata_outputs_annot/annotation_reference'

data_score_rf_average = df_average(data_score_rf)
print()
sorted_score_codes = get_metabolic_codes(data_score_rf_average, data_labels)

data_score_codes_names = add_reaction_names(sorted_score_codes)
print()
data_score_codes_names_plus_otus = find_relevant_otus(data_score_codes_names, path_to_annotation_ref)

#print(data_score_codes_names_plus_otus)
if(not os.path.isfile(pipeline_path + '/Post-processing/outputs')):
        # Emulate mkdir -p (no error if folder exists)
        try:
            os.mkdir(pipeline_path + '/Post-processing/outputs')
        except OSError as e:
            if(e.errno != 17):
                raise e

if args.tfigm:
    data_score_codes_names_plus_otus.to_csv(pipeline_path + '/Post-processing/outputs/'+dataset_name+'_tfigm_annots_rf.csv', index = None)
elif args.scaling:
    data_score_codes_names_plus_otus.to_csv(pipeline_path + '/Post-processing/outputs/'+dataset_name+'_scaled_annots_rf.csv', index = None)
else:
    data_score_codes_names_plus_otus.to_csv(pipeline_path + '/Post-processing/outputs/'+dataset_name+'_annots_rf.csv', index = None)
