import pandas as pd
import numpy as np
from collections import defaultdict
from Bio.ExPASy import Enzyme
from goatools import obo_parser
import requests
import os

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

def add_reaction_names(dataframe):
    '''
    This function queries the OBO and ExPASy data banks to get the names of the IDs
    '''
    
    #bar = IncrementalBar('Recovering Names', max = len(dataframe.iloc[:,1]))

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

    for id in dataframe.iloc[:,0]:
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
        
        #bar.next()
    #print(len(reaction_names))
    dataframe["Names"] = reaction_names
    #dataframe.columns = ["Average_importance","ID","Name"]
    return dataframe


###GATHER ALL SIGNIFICANT LISTS

link_to_test_folder = pipeline_path+'/Meta_Outputs/'+args.data_ref
iterations=int(args.nb_iterations)
dict_lists = {}
for i in range(iterations):
    iter_folder = args.data_ref+'_'+str(i+1)
    for run_nb in os.listdir(link_to_test_folder+'/'+iter_folder+'/Selection_outputs'):
        otu_db = pd.read_csv(link_to_test_folder+'/'+iter_folder+'/Selection_outputs/'+run_nb+'/OTU_'+dataset_name+'.csv')
        score_db = pd.read_csv(link_to_test_folder+'/'+iter_folder+'/Selection_outputs/'+run_nb+'/scores_'+dataset_name+'.csv')
        
        signif_list_sofa = score_db.truncate(after = score_db.loc[score_db['ID']=='CUTOFF'].index.values[0]-1)['ID'].values
        signif_list_otu = otu_db.truncate(after = otu_db.loc[otu_db['ID']=='CUTOFF'].index.values[0]-1)['ID'].values
        print(signif_list_otu)

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


###CREATE CORE AND META LISTS

dict_core_meta = defaultdict(defaultdict)
if not os.path.exists(link_to_test_folder+'/core_and_meta_all_runs'):
    os.makedirs(link_to_test_folder+'/core_and_meta_all_runs')

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

        if profile =='SoFA':
            core = add_reaction_names(core)
            meta = add_reaction_names(meta)

        core.to_csv(link_to_test_folder +'/core_and_meta_all_runs/core_'+profile+'_'+dataset_name+'_iteration_'+run+'.csv', index=0)
        meta.to_csv(link_to_test_folder +'/core_and_meta_all_runs/meta_'+profile+'_'+dataset_name+'_iteration_'+run+'.csv', index=0)
