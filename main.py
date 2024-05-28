from DMmodif_export_test_ver import *

from esmecata.proteomes import retrieve_proteomes
from esmecata.clustering import make_clustering
from esmecata.annotation import annotate_proteins
from esmecata.eggnog import annotate_with_eggnog

import os
import csv
import shutil
import psutil
import sys
import traceback
import requests
from operator import methodcaller
from Bio.ExPASy import Enzyme
from goatools import obo_parser
from Bio import SeqIO
from Bio import __version__ as biopython_version
import multiprocessing

import pandas as pd
import numpy as np
import subprocess
import math
from scipy.stats import mannwhitneyu
from sklearn.model_selection import train_test_split
from kneebow.rotor import Rotor

from collections import Counter, defaultdict

import seaborn as sns
import matplotlib.pyplot as plt

import logging
from datetime import datetime

from tqdm import tqdm, trange
from tqdm.contrib.logging import logging_redirect_tqdm

import argparse

parser = argparse.ArgumentParser(description='A program that averages the RF importance scores of the functional annotations, and associates them to OTUs')

parser.add_argument("-d","--dataset_name", help="Name of the dataset", required=True)
parser.add_argument("-t", "--treatment", default=None, help="Data treatment for the functional table (can be: 'tf_igm', default: no treatment)")
parser.add_argument("-s", "--scaling", default=None, help="Scaling method to apply to the taxonomic table (can be: 'relative', default: no scaling)")
parser.add_argument("-i", "--iterations", default=5, help="Number of iterations of the method (default: 5 iterations)")
parser.add_argument("-f", "--forests", default=20, help="Amount of trained classifiers per iteration of the command (default: 20 forests)")
parser.add_argument("-r", "--runs", default=10, help="Amount of pipeline runs (default: 10 runs)")
parser.add_argument("--eggnog", default=False, help="Path to the eggnog database for the EsMeCaTa pipeline. If not given, the pipeline will be launhed with the 'UniProt' workflow by default.")
parser.add_argument("--annotations_only", default=False, action='store_true', help="This is a flag that signals that the input is a functional table. If True, all steps involving taxonomic tables will be skipped, and SPARTA will iteratively classify and select on the given functional table alone.")
parser.add_argument("--reference_test_sets", default=False, action='store_true', help="This option allows the user to give their own test sets to be used during classification.")
parser.add_argument("--esmecata_relaunch", default=False, action='store_true', help="This option allows the user to force a re-run of the EsMeCaTa pipeline over an already existing output. This is particularly useful if a previous run of the pipeline was botched at this step.")
parser.add_argument("--keep_temp", default=False, action='store_true', help="This option allows the user to keep the contents of the 'Outputs_temp' folder at the end of the run.")

pd.options.mode.chained_assignment = None

args_passed = parser.parse_args()


logger = logging.getLogger(__name__)

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    with logging_redirect_tqdm():
        for i in trange(9):
            if i == 4:
                logger.info("console logging redirected to `tqdm.write()`")
    # logging restored


### FORMATTING STEP
def pre_formatting(dataset_full, dataset_name, data_ref_output_name, pipeline_path):
    '''
    This function pre-formats the taxonomic data by:
        - Creating an input for the EsMeCaTa pipeline (OTU name + taxonomic description)
        - Removing the metadata from the table
    '''

    #Import the original file's description of the OTU taxonomies
    
    dataset_compos = dataset_full.loc[dataset_full[0].str.contains('k__')]

    #Replace the separators to create compatibility with EsMeCaTa
    new_cols = list(dataset_compos[0].values)

    new_cols = [line.replace('k__','') for line in new_cols]
    new_cols = [line.replace('|p__',';') for line in new_cols]
    new_cols = [line.replace("|c__",";") for line in new_cols]
    new_cols = [line.replace("|o__",";") for line in new_cols]
    new_cols = [line.replace("|f__",";") for line in new_cols]
    new_cols = [line.replace("|g__",";") for line in new_cols]
    new_cols = [line.replace("|s__",";") for line in new_cols]
    new_cols = [line.replace("_"," ") for line in new_cols]

    #Formatting the dataset for DeepMicro and EsMeCaTa, and giving the OTUs standardised aliases (OTU_x)
    otu_names = ['OTU_' + str(i+1) for i in range(len(new_cols))]
    formatted_data = pd.DataFrame({'observation_name':otu_names, 'taxonomic_affiliation':new_cols}, index = None)

    
    dataset_compos.index = otu_names
    dataset_compos = dataset_compos.iloc[:,1:]
    dataset_compos.columns=dataset_full[dataset_full[0] == 'sampleID'].values[0][1:]


    cols_to_change = dataset_compos.columns
    dataset_compos[cols_to_change] = dataset_compos[cols_to_change].astype(float)

    #Writing EsMeCaTa input
    formatted_data.to_csv(pipeline_path+"/Outputs_temp/"+data_ref_output_name+"/SoFA_calculation/"+dataset_name+".tsv", sep = "\t", index = None)
    

    return formatted_data, dataset_compos

def data_to_deepmicro(dataset_compos):
    '''
    This function transforms taxonomic and functional profiles to a format compatible with DeepMicro
    '''

    dataset_compos_deepmicro = dataset_compos.transpose()
    #dataset_compos_deepmicro.columns = dataset_compos_deepmicro.iloc[0]
    #dataset_compos_deepmicro.to_csv(pipeline_path+"/Outputs_temp/"+data_ref_output_name+"/DeepMicro_data/entree_DeepMicro_"+dataset_name+"_OTU.csv", sep = ",", header = False, index=False)
    return dataset_compos_deepmicro


def absolute_to_relative(db):
    '''
    This function transforms a dataframe of absolute abundances into relative abundances
    '''
    for col_num in db:
        col_sum = db[col_num].sum(axis = 0)
        db[col_num] = [x/col_sum for x in db[col_num].values]
    return db





###ESMECATA LAUNCHER (TEMPORARY?)

def personalized_make_clustering(proteome_folder, output_folder, nb_cpu=1, clust_threshold=1, mmseqs_options=None, linclust=None, remove_tmp=None):
    """From the proteomes found by esmecata proteomes, create protein cluster for each taxonomic affiliations. 
        Personnalization: Certain features were removed, notably for graphical representations of the outputs.
    Args:
        proteome_folder (str): pathname to folder from esmecata folder
        output_folder (str): pathname to the output folder
        nb_cpu (int): number of CPUs to be used by mmseqs
        clust_threshold (float): threshold to select protein cluster according to the representation of protein proteome in the cluster
        mmseqs_options (str): use alternative mmseqs option
        linclust (bool): use linclust
        remove_tmp (bool): remove the tmp files
    """
    starttime = time.time()
    logger.info('|EsMeCaTa|clustering| Begin clustering.')

    # Check if mmseqs is in path.
    mmseqs_path = which('mmseqs')
    if not mmseqs_path:
        logger.critical('|EsMeCaTa|clustering| mmseqs not available in path, esmecata will not be able to cluster the proteomes.')
        sys.exit(1)

    if not is_valid_dir(proteome_folder):
        logger.critical('|EsMeCaTa|clustering| Input must be a folder %s.', proteome_folder)
        sys.exit(1)

    # Use the proteomes folder created by retrieve_proteome.py.
    proteome_tax_id_pathname = os.path.join(proteome_folder, 'proteome_tax_id.tsv')

    if not is_valid_path(proteome_tax_id_pathname):
        logger.critical(f"|EsMeCaTa|clustering| Missing output from esmecata proteomes in {proteome_tax_id_pathname}.")
        sys.exit(1)

    reference_proteins_path = os.path.join(output_folder, 'reference_proteins')
    is_valid_dir(reference_proteins_path)

    cluster_founds_path = os.path.join(output_folder, 'cluster_founds')
    is_valid_dir(cluster_founds_path)

    computed_threshold_path = os.path.join(output_folder, 'computed_threshold')
    is_valid_dir(computed_threshold_path)

    already_performed_clustering = [reference_protein_file.replace('.tsv', '') for reference_protein_file in os.listdir(reference_proteins_path)]

    # Create a dictionary with observation_name as key and the pathname to the proteomes associated to this observation_name as value.
    observation_name_fasta_files = {}
    with open(proteome_tax_id_pathname, 'r') as proteome_tax_file:
        csvreader = csv.DictReader(proteome_tax_file, delimiter='\t')
        for line in csvreader:
            proteomes = line['proteome'].split(',')
            tax_name = line['name'].replace(' ','_')
            proteomes_path = [os.path.join(proteome_folder, 'proteomes', proteome+'.faa.gz') for proteome in proteomes]
            if tax_name not in already_performed_clustering:
                observation_name_fasta_files[tax_name] = proteomes_path
            else:
                logger.info('|EsMeCaTa|clustering| Already performed clustering for %s.', tax_name)
    
    

    is_valid_dir(output_folder)

    # Create metadata file.
    clustering_metadata = {}
    clustering_metadata['tool_options'] = {'proteome_folder': proteome_folder, 'output_folder': output_folder, 'nb_cpu':nb_cpu,
                                        'clust_threshold':clust_threshold, 'mmseqs_options': mmseqs_options, 'linclust':linclust,
                                        'remove_tmp': remove_tmp}

    clustering_metadata['tool_dependencies'] = {}
    subprocess_output = subprocess.check_output(['mmseqs', 'version'])
    mmseqs_version = subprocess_output.decode('utf-8')
    clustering_metadata['tool_dependencies']['mmseqs_version'] = mmseqs_version
    clustering_metadata['tool_dependencies']['mmseqs_path'] = mmseqs_path
    clustering_metadata['tool_dependencies']['python_package'] = {}
    clustering_metadata['tool_dependencies']['python_package']['Python_version'] = sys.version
    clustering_metadata['tool_dependencies']['python_package']['biopython'] = biopython_version
    clustering_metadata['tool_dependencies']['python_package']['esmecata'] = esmecata_version

    # Create tmp folder for mmseqs analysis.
    mmseqs_tmp_path = os.path.join(output_folder, 'mmseqs_tmp')
    is_valid_dir(mmseqs_tmp_path)

    # Create output folder containing shared representative proteins.
    reference_proteins_representative_fasta_path = os.path.join(output_folder, 'reference_proteins_representative_fasta')
    is_valid_dir(reference_proteins_representative_fasta_path)

    reference_proteins_consensus_fasta_path = os.path.join(output_folder, 'reference_proteins_consensus_fasta')
    is_valid_dir(reference_proteins_consensus_fasta_path)

    proteome_taxon_id_file = os.path.join(proteome_folder, 'proteome_tax_id.tsv')
    clustering_taxon_id_file = os.path.join(output_folder, 'proteome_tax_id.tsv')

    if os.path.exists(clustering_taxon_id_file):
        if not os.path.samefile(proteome_taxon_id_file, clustering_taxon_id_file):
            os.remove(clustering_taxon_id_file)
            shutil.copyfile(proteome_taxon_id_file, clustering_taxon_id_file)
    else:
        shutil.copyfile(proteome_taxon_id_file, clustering_taxon_id_file)

    proteomes_taxa_names = get_proteomes_tax_name(proteome_taxon_id_file)

    all_tax_names = set(list(proteomes_taxa_names.values()))

    # For each OTU run mmseqs easy-cluster on them to found the clusters that have a protein in each proteome of the OTU.
    # We take the representative protein of a cluster if the cluster contains a protein from all the proteomes of the OTU.
    # If this condition is not satisfied the cluster will be ignored.
    # Then a fasta file containing all the representative proteins for each OTU is written in representative_fasta folder.

    for proteomes_tax_name in all_tax_names:
        # Get proteomes associated with taxon name.
        observation_name_proteomes = observation_name_fasta_files[proteomes_tax_name]

        # Change space with '_' to avoid issue.
        proteomes_tax_name = proteomes_tax_name.replace(' ', '_')
        # If the computed threshold file exists, mmseqs has already been run.
        mmseqs_tmp_cluster = os.path.join(mmseqs_tmp_path, proteomes_tax_name)
        # Run mmseqs on organism.
        # Delete previous mmseqs2 run if it exists to avoid overwritting issues.
        if os.path.exists(mmseqs_tmp_cluster):
            shutil.rmtree(mmseqs_tmp_cluster)
        mmseqs_tmp_clustered_tabulated, mmseqs_tmp_representative_fasta, mmseqs_consensus_fasta = run_mmseqs(proteomes_tax_name, observation_name_proteomes, mmseqs_tmp_path, nb_cpu, mmseqs_options, linclust)

        # Extract protein clusters from mmseqs results.
        cluster_proteomes_output_file = os.path.join(cluster_founds_path, proteomes_tax_name+'.tsv')
        protein_clusters = extrat_protein_cluster_from_mmseqs(mmseqs_tmp_clustered_tabulated, cluster_proteomes_output_file)

        # Compute proteome representativeness ratio.
        computed_threshold_file = os.path.join(computed_threshold_path, proteomes_tax_name+'.tsv')
        number_proteomes, rep_prot_organims, computed_threshold_cluster = compute_proteome_representativeness_ratio(protein_clusters,
                                                                                                                    observation_name_proteomes, computed_threshold_file)

        # Filter protein cluster for each protein cluster.
        cluster_proteomes_filtered_output_file = os.path.join(reference_proteins_path, proteomes_tax_name+'.tsv')
        protein_cluster_to_keeps = filter_protein_cluster(protein_clusters, number_proteomes, rep_prot_organims, computed_threshold_cluster,
                                                        clust_threshold, cluster_proteomes_filtered_output_file)

        logger.info('|EsMeCaTa|clustering| %d protein clusters kept for %s.', len(protein_cluster_to_keeps), proteomes_tax_name)

        # Create BioPython records with the representative proteins kept.
        new_records = [record for record in SeqIO.parse(mmseqs_tmp_representative_fasta, 'fasta') if record.id.split('|')[1] in protein_cluster_to_keeps]

        # Do not create fasta file when 0 sequences were kept.
        if len(new_records) > 0:
            # Create output proteome file for OTU.
            representative_fasta_file = os.path.join(reference_proteins_representative_fasta_path, proteomes_tax_name+'.faa')
            SeqIO.write(new_records, representative_fasta_file, 'fasta')
        else:
            logger.info('|EsMeCaTa|clustering| 0 protein clusters %s, no fasta created.', proteomes_tax_name)
        del new_records

        # Create BioPython records with the consensus proteins kept.
        consensus_new_records = [record for record in SeqIO.parse(mmseqs_consensus_fasta, 'fasta') if record.id.split('|')[1] in protein_cluster_to_keeps]

        # Do not create fasta file when 0 sequences were kept.
        if len(consensus_new_records) > 0:
            # Create output proteome file for OTU.
            consensus_fasta_file = os.path.join(reference_proteins_consensus_fasta_path, proteomes_tax_name+'.faa')
            SeqIO.write(consensus_new_records, consensus_fasta_file, 'fasta')
        else:
            logger.info('|EsMeCaTa|clustering| 0 protein clusters %s, no fasta created.', proteomes_tax_name)
        del consensus_new_records

        if remove_tmp:
            shutil.rmtree(mmseqs_tmp_cluster)

    # Compute number of protein clusters kept.
    stat_file = os.path.join(output_folder, 'stat_number_clustering.tsv')
    compute_stat_clustering(output_folder, stat_file)
    output_figure_file = os.path.join(output_folder, 'representativeness_clustering_ratio.svg')
    #create_proteome_representativeness_lineplot(clustering_taxon_id_file, computed_threshold_path, output_figure_file)

    proteome_ratio_lineplots_path = os.path.join(output_folder, 'proteome_ratio_lineplots')
    is_valid_dir(proteome_ratio_lineplots_path)
    #create_proteome_representativeness_lineplot_per_taxon_rank(clustering_taxon_id_file, computed_threshold_path, proteome_ratio_lineplots_path)

    endtime = time.time()
    duration = endtime - starttime
    clustering_metadata['esmecata_clustering_duration'] = duration
    clustering_metadata_file = os.path.join(output_folder, 'esmecata_metadata_clustering.json')
    if os.path.exists(clustering_metadata_file):
        metadata_files = [metadata_file for metadata_file in os.listdir(output_folder) if 'esmecata_metadata_clustering' in metadata_file]
        clustering_metadata_file = os.path.join(output_folder, 'esmecata_metadata_clustering_{0}.json'.format(len(metadata_files)))
        with open(clustering_metadata_file, 'w') as ouput_file:
            json.dump(clustering_metadata, ouput_file, indent=4)
    else:
        with open(clustering_metadata_file, 'w') as ouput_file:
            json.dump(clustering_metadata, ouput_file, indent=4)

    logger.info('|EsMeCaTa|clustering| Clustering complete in {0}s.'.format(duration))



def proteome_check(esmecata_prot_out):
    '''
    This function checks if a proteome has been downloaded for each OTU given as input to EsMeCaTa
    '''

    path_to_proteome_count = esmecata_prot_out+'/proteome_tax_id.tsv'

    proteome_count_db = pd.read_csv(path_to_proteome_count, sep='\t')

    list_of_proteomes = []

    for proteome_comma in proteome_count_db['proteome'].values:
        for proteome in proteome_comma.split(','):
            if proteome not in list_of_proteomes:
                list_of_proteomes.append(proteome)

    proteome_count = len(list_of_proteomes)

    dir_path = esmecata_prot_out+'/proteomes'

    count = 0

    for path in os.listdir(dir_path):
        # check if current path is a file
        if os.path.isfile(os.path.join(dir_path, path)):
            count += 1

    all_proteomes_checked = (count == proteome_count)
    logger.info("All proteomes downloaded: "+ str(all_proteomes_checked))
    return all_proteomes_checked


def get_proteomes_tax_name(proteomes_tax_id):
    """ Extract tax_id associated with observation name.

    Args:
        proteomes_tax_id (str): pathname to the proteomes_tax_id file.

    Returns:
        proteomes_taxa_names (dict): dict containing observation names (as key) associated with tax name used for proteomes (as value)
    """
    proteomes_taxa_names = {}
    with open(proteomes_tax_id, 'r') as proteome_tax_file:
        csvreader = csv.DictReader(proteome_tax_file, delimiter='\t')
        for line in csvreader:
            observation_name = line['observation_name']
            tax_name = line['name']
            proteomes_taxa_names[observation_name] = tax_name.replace(' ', '_')

    return proteomes_taxa_names


def check_annotation(reference_proteins_consensus_fasta_folder, annotation_reference_folder, proteomes_tax_id):
    proteomes_taxa_names = get_proteomes_tax_name(proteomes_tax_id)
    annotated_observation_names = [observation_name.replace('.tsv', '') for observation_name in os.listdir(annotation_reference_folder)]
    expected_annotations = [observation_name for observation_name in proteomes_taxa_names
                                                if os.path.exists(os.path.join(reference_proteins_consensus_fasta_folder, proteomes_taxa_names[observation_name] + '.faa'))]
    
    logger.info('Annotated observations: '+str(annotated_observation_names))
    logger.info('Expected observations: '+str(expected_annotations))

    empty_output_otus = []
    for otu_name in proteomes_taxa_names.keys():    
        if otu_name not in expected_annotations:
            empty_output_otus.append(otu_name)

    return len(annotated_observation_names) == len(expected_annotations), empty_output_otus





def esmecata_plus_check(esmecata_input, pipeline_path, data_ref_output_name, dataset_name, eggnog_path=False):
    '''
    This function runs EsMeCaTa on the input files. Checks are made at the end of the 'Proteomes' and 'Annotation' processes 
    '''
    nb_cpu_available = multiprocessing.cpu_count()
    if nb_cpu_available <= 3:
        nb_cpu_available = 1
    if nb_cpu_available > 3:
        nb_cpu_available -= 2

    input_location = pipeline_path+'/Outputs_temp/'+data_ref_output_name+'/SoFA_calculation/'+dataset_name+'.tsv'
    esmecata_prot_out = pipeline_path+'/EsMeCaTa_outputs/'+dataset_name+'/esmecata_outputs_proteomes'
    esmecata_cluster_out = pipeline_path+'/EsMeCaTa_outputs/'+dataset_name+'/esmecata_outputs_clustering'
    esmecata_annots_out = pipeline_path+'/EsMeCaTa_outputs/'+dataset_name+'/esmecata_outputs_annots'
    

    ## EsMeCaTa
    count_check = False
    retries = 0
    while (not count_check) and (retries<20):
        try:
            retrieve_proteomes(input_location, esmecata_prot_out, option_bioservices=True)
        except:
            pass
            
        count_check = proteome_check(esmecata_prot_out)
        retries+=1
    
    if retries >= 20:
        raise Exception("EsMeCaTa has failed 20 times in a row. A connexion error is likely. Aborting...")

    
    if not os.path.exists(esmecata_cluster_out+'/stat_number_clustering.tsv'):
        if os.path.exists(esmecata_cluster_out):
            logger.info('Previous incomplete iteration of the clustering step found: deleting and starting over')
            shutil.rmtree(esmecata_cluster_out, ignore_errors=True)
        make_clustering(esmecata_prot_out, esmecata_cluster_out, nb_cpu = nb_cpu_available, mmseqs_options=None, clust_threshold=0.5, linclust=None, remove_tmp=True)
    else:
        logger.info('Clustering step already done, moving to annotation.')

    try:
        if eggnog_path == False:
            
            annotate_proteins(esmecata_cluster_out, esmecata_annots_out, uniprot_sparql_endpoint=None,
                            propagate_annotation=1, uniref_annotation=None, expression_annotation=None)
        else:
            annotate_with_eggnog(esmecata_cluster_out, esmecata_annots_out, eggnog_path, nb_cpu = nb_cpu_available)
    except:
        pass
    
    


    ## Check
    retries = 0
    annotation_reference_folder = esmecata_annots_out+'/annotation_reference/'
    reference_proteins_consensus_fasta_folder = esmecata_cluster_out+'/reference_proteins_consensus_fasta/'
    proteomes_tax_id = esmecata_cluster_out+'/proteome_tax_id.tsv'
    check_outputs, empty_ids = check_annotation(reference_proteins_consensus_fasta_folder, annotation_reference_folder, proteomes_tax_id)
    logger.info('All annotations found: '+str(check_outputs))
    while (not check_outputs) and (retries<20):
        try:
            if eggnog_path == False:
                
                annotate_proteins(esmecata_cluster_out, esmecata_annots_out, uniprot_sparql_endpoint=None,
                                propagate_annotation=1, uniref_annotation=None, expression_annotation=None)
            else:
                annotate_with_eggnog(esmecata_cluster_out, esmecata_annots_out, eggnog_path, nb_cpu = 10)
        except:
            pass
        retries += 1
        check_outputs, empty_ids = check_annotation(reference_proteins_consensus_fasta_folder, annotation_reference_folder, proteomes_tax_id)
        logger.info('All annotations found: '+str(check_outputs))

    if retries >= 20:
        raise Exception("EsMeCaTa has failed 20 times in a row. A connexion error is likely. Aborting...")

    empty_output_df = pd.DataFrame(columns=['protein_cluster','cluster_members','gene_name','GO','EC','KEGG_reaction'])
    for otu_id in empty_ids:
        empty_output_df.to_csv(annotation_reference_folder+'/'+otu_id+'.tsv', sep='\t', index=False)
    
    ## Clean

    shutil.rmtree(esmecata_prot_out, ignore_errors=True)
    shutil.rmtree(esmecata_cluster_out, ignore_errors=True)


### CALCULATING AND TRANSFORMING SoFAs


def tf_igm(dataline, lbd):
    '''
    This function applies the TF-IGM transformation to a row of a dataset
    '''
    line_sorted = np.array(sorted(dataline, reverse = True))
    fk1 = line_sorted[0]
    multipliers_array = np.array([i+1 for i in range(len(line_sorted))])
    fk_sum_list = line_sorted * multipliers_array
    fk_sum = np.sum(fk_sum_list)
    
    for data_index in range(len(dataline)):
        if fk_sum != 0:
            dataline[data_index] = math.sqrt(dataline[data_index]) * (1 + lbd * fk1 / fk_sum)
        else:
            dataline[data_index] = 0
    
    return(dataline)


def tf_igm_apply(deepmicro_sofa):
    '''
    This function runs the TF-IGM transformation function on all lines of a dataset
    '''
    
    lbd = 7

    sum_vector = deepmicro_sofa.sum(axis=1)

    freq_matrix = deepmicro_sofa.div(sum_vector, axis=0).fillna(0)

    for i in deepmicro_sofa.columns.values:
        #logger.info(i)
        score_test = tf_igm(freq_matrix.loc[:,i].values, lbd)
        deepmicro_sofa[i] = score_test

    return deepmicro_sofa
    



def sofa_calculation(pipeline_path, dataset_name, data_ref_output_name, otu_table_stripped, args_passed):
    '''
    This function calculates the scores of the functional annotations from EsMeCaTa's outputs and the original OTU abundances.

    INPUTS: 
        - original args
        - OTU table
    
    OUTPUTS:
        - sofa_table: a table containing all of the scores of functional annotations
        - deepmicro_sofa: the same table, but formatted as a DeepMicro input. If a transformation is given as an argument, 
    '''

    esmecata_output_path = pipeline_path+'/EsMeCaTa_outputs/'+dataset_name+'/esmecata_outputs_annots/annotation_reference'
    otus_GOs = {}
    otus_ECs = {}
    all_gos = []
    all_ecs = []

    dir_list = os.listdir(esmecata_output_path)

    #Parse the EsMeCaTa outputs and list all annotations expressed per OTU, along with the number of proteins that express them
    for annot_file in tqdm(dir_list, desc="Parsing EsMeCaTa's outputs"):
        base_file = os.path.basename(annot_file)
        if base_file[0] != '.':
            base_filename = os.path.splitext(base_file)[0]
            annot_file_path = os.path.join(esmecata_output_path, annot_file)
            df = pd.read_csv(annot_file_path, sep='\t')
            df = df.replace(np.nan, '')
            go_series = df.GO.map(lambda x: [i.strip() for i in x.split(',')]).apply(pd.Series)
            if go_series.empty is False:
                otu_go_terms = list(go_series.stack())
                otu_go_terms = [go for go in otu_go_terms if go != '']
                otu_go_terms_counter = Counter(otu_go_terms)
                otus_GOs[base_filename] = otu_go_terms_counter
                all_gos.extend(otu_go_terms_counter.keys())

            ec_series = df.EC.map(lambda x: [i.strip() for i in x.split(',')]).apply(pd.Series)
            if ec_series.empty is False:
                otu_ec_numbers = list(ec_series.stack())
                otu_ec_numbers = [ec for ec in otu_ec_numbers if ec != '']
                otu_ec_numbers_counter = Counter(otu_ec_numbers)
                otus_ECs[base_filename] = otu_ec_numbers_counter
                all_ecs.extend(otu_ec_numbers_counter.keys())



    all_gos = list(set(all_gos))
    all_ecs = list(set(all_ecs))
    all_annots = all_gos + all_ecs
    otu_annots = {}


    for otu in tqdm(otu_table_stripped.index, desc='Checking out OTU abundances'):
        otu_annots[otu] = {}
        if otu in otus_ECs:
            otu_annots[otu].update(otus_ECs[otu])
        if otu in otus_GOs:
            otu_annots[otu].update(otus_GOs[otu])
    

    # Compute for each annotation its abundance
    sofa_table = pd.DataFrame(all_annots)
    sofa_table.set_index(0, inplace=True)


    for sample in tqdm(otu_table_stripped.columns, desc='Calculating abundances'):
        if sample != 'OTU':
            otu_annots_dataframe = pd.DataFrame(otu_annots)
            
            for col in otu_annots_dataframe.columns: 
                otu_annots_dataframe[col] = otu_annots_dataframe[col] * otu_table_stripped[sample].loc[col]
            sofa_table[sample] = otu_annots_dataframe.sum(axis=1)


    sofa_table.to_csv(pipeline_path+'/Meta-Outputs/'+data_ref_output_name+'/SoFA_table_'+dataset_name+'.csv')
    deepmicro_sofa = sofa_table.transpose()
    

    if args_passed.treatment == 'tf_igm':
        deepmicro_sofa = tf_igm_apply(deepmicro_sofa)

    #deepmicro_sofa.to_csv(pipeline_path+'/Outputs_temp/'+data_ref_output_name+'/DeepMicro_data/entree_DeepMicro_'+dataset_name+'.csv', sep=',', header = None, index = None)
    
    return sofa_table, deepmicro_sofa

### AVERAGING PER GROUP AND GETTING INFORMATION ABOUT FUNCTIONAL ANNOTATIONS

def add_reaction_names(list_of_annots, pipeline_path):
    '''
    This function queries the OBO and ExPASy data banks to get the names of the IDs
    '''
    
    data_folder = pipeline_path + '/data'

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


    for id in tqdm(list_of_annots, desc='Checking out annotation names...'):
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
            

    dataframe = pd.DataFrame(list(zip(list_of_annots, reaction_names)), columns=["ID","Name"])
                             
    return dataframe


def add_otu_names(otu_list, otu_name_df):
    '''
    This function formats a dataframe linking each OTU to its full name
    '''

    otu_names = []
    for otu in tqdm(otu_list, desc="Checking out OTU names..."):
        otu_name_translated = otu_name_df[otu_name_df['observation_name'] == otu]['taxonomic_affiliation'].values[0]
        otu_names.append(otu_name_translated)
    

    dataframe = pd.DataFrame(list(zip(otu_list, otu_names)), columns=["ID","Name"])

    return dataframe


def find_relevant_otus(dataframe, path, otu_name_df):
    '''
    This function finds and names all OTUs associated with each annotation
    '''

    found_otu = {}
    found_otu_named = {}
    for filename in tqdm([f for f in os.listdir(path) if not f.startswith('.')], desc="Linking OTUs to annotations..."):
        otu_data = pd.read_csv(path + "/" + filename, sep ='\t', keep_default_na=False)
        otu_name = filename.replace(".tsv","")
        otu_name_translated = otu_name_df[otu_name_df['observation_name'] == otu_name]['taxonomic_affiliation'].values[0]
        otu_name_translated_species = otu_name_translated.split(';')[-1]

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
                    found_otu_named[id] = []
                found_otu[id].append(otu_name)
                found_otu_named[id].append(otu_name_translated_species)


    dataframe["Linked_OTUs"] = dataframe["ID"].map(found_otu)
    dataframe["Named_linked_OTUs"] = dataframe["ID"].map(found_otu_named)
    return dataframe

def find_relevant_reactions(dataframe, path):
    '''
    This function finds all annotations associated with each OTU
    '''

    found_reac = {}
    for filename in tqdm([f for f in os.listdir(path) if not f.startswith('.')], desc="Linking annotations to OTUs..."):
        
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

    dataframe["Linked_annotations"] = dataframe["ID"].map(found_reac)
    return dataframe



def get_info_annots(score_db, esmecata_input, pipeline_path, dataset_name, args_passed):
    

    list_of_annots = list(score_db.index)
    annots_with_names = add_reaction_names(list_of_annots, pipeline_path)

    path = pipeline_path+'/EsMeCaTa_outputs/'+dataset_name+'/esmecata_outputs_annots/annotation_reference'

    if not args_passed.annotations_only:
        annots_with_names_and_associated_otus = find_relevant_otus(annots_with_names, path, esmecata_input)
    else:
        annots_with_names_and_associated_otus = annots_with_names

    return(annots_with_names_and_associated_otus)


def get_info_taxons(otu_db, esmecata_input, pipeline_path, dataset_name):
    
    list_of_otus = list(otu_db.index)

    annots_with_names = add_otu_names(list_of_otus, esmecata_input)

    path = pipeline_path+'/EsMeCaTa_outputs/'+dataset_name+'/esmecata_outputs_annots/annotation_reference'

    annots_with_names_and_associated_otus = find_relevant_reactions(annots_with_names, path)

    return(annots_with_names_and_associated_otus)






def average_per_group(score_dataframe, label_refs, label_value = None):
    '''
    This function averages the presence of a functional annotation within a given label group
    '''
    
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


def average_count_per_group(count_dataframe,  label_refs, label_value=None):
    '''
    This function averages the presence of a taxon within a given label group
    '''


    if label_value != None:
        count_group = count_dataframe.loc[:, [i[0] == label_value for i in label_refs.values]].astype(float)
    else:
        count_group = count_dataframe.astype(float)
    
    count_group["average"] = count_group.mean(axis=1)


    return(count_group)


###PROCESSES OF PIPELINE ITERATION

def create_test_train_subsets(full_dataset, indices):
    '''
    This function separates a dataset into a test and a train subset, determined by the input indices
    '''

    dataset_train = full_dataset.drop(indices, axis=0)
    dataset_test = full_dataset.loc[indices]

    return(dataset_train, dataset_test)


def separate_test_train(labels_full):
    if not labels_full.values.shape[1] > 1:
        label_flatten = labels_full.values.reshape((labels_full.values.shape[0]))
    else:
        raise Exception("Please submit the label file in the form of an ordered vector")

    indices = np.arange(len(label_flatten))

    y_train, y_test, sample_train, sample = train_test_split(label_flatten, indices, test_size=0.2, stratify=label_flatten)

    return y_train, y_test, sample



def run_deep_micro(set_test, set_train, label_test, label_train, dataset_name, iteration_nb, run_nb, pipeline_path, args_passed, data_dir, profile):
    
    try:

        ##IMPORT SEPARATE TEST DATASET

        Xtest_ext = set_test.values.astype(np.float64)

        Ytest_ext = label_test.astype(int)

        best_auc = 0
        threshold_opt = 0
        perf_dict = {}
        repeats = int(args_passed.forests)
        best_feature_records = []

        # hyper-parameter grids for classifiers
    

        if repeats > 1:
            for i in tqdm(range(repeats), desc="Training models for the "+profile+" profile (Run: "+str(run_nb)+", Iteration: "+str(iteration_nb)+")"):
                #logger.info('BEST FEATURE RECORDS:', best_feature_records)
                best_auc, threshold_opt, perf_dict, best_feature_records = run_exp(i, best_auc, threshold_opt, perf_dict, Xtest_ext, Ytest_ext, set_train, label_train, dataset_name, best_feature_records, data_dir)
        else:
            best_auc, threshold_opt, perf_dict, best_feature_records = run_exp(42, best_auc, threshold_opt, perf_dict, Xtest_ext, Ytest_ext, set_train, label_train, dataset_name, best_feature_records, data_dir)

        perf_df = pd.DataFrame.from_dict(perf_dict, orient='index', dtype=None, columns=['Best parameters', 'Validation set indices','Threshold', 'Training performance', 'Validation performance', 'Test performance'])
    
        return perf_df, best_feature_records

    except OSError as error:
        raise Exception("DeepMicro encountered a problem.")

        
def inflexion_cutoff(datatable):
    
    datatable_ordered = datatable.sort_values(by="Average", ascending=False)

    data = datatable_ordered["Average"].values

    data_rotor = []
    for i in range(len(data)):
        data_rotor.append([i+1,data[i]])

    rotor = Rotor()
    rotor.fit_rotate(data_rotor)

    elbow_idx = rotor.get_elbow_index()

    datatable_truncated = datatable_ordered["Average"].iloc[:elbow_idx+1]

    return datatable_truncated

### PIPELINE LAUNCHER FUNCTIONS

def formatting_step(dataset_name, data_ref_output_name, pipeline_path, args_passed):
    '''
    This script runs the steps to format the data into:
        - a metadata-less version of the original OTU table (otu_table_stripped)
        - an input to launch the esmecata pipeline afterwards (esmecata_input)
        - a taxonomic input for DeepMicro (deepmicro_otu)
    
        If --annotation only: the functional dataset will be called 'otu_table' in context of this function, and will be imported and converted to deepmicro input
    '''

    label_file = pd.read_csv(pipeline_path+"/Inputs/Label_"+dataset_name+".csv", header = None, sep = "\t")


    if not args_passed.annotations_only:
        dataset_full = pd.read_csv(pipeline_path+"/Inputs/"+dataset_name+".txt", header = None, sep = "\t")
        esmecata_input, otu_table_stripped = pre_formatting(dataset_full, dataset_name, data_ref_output_name, pipeline_path)
        if args_passed.scaling == "relative":
            otu_table_stripped = absolute_to_relative(otu_table_stripped)
        #Writing standardised version of the dataset
        otu_table_stripped.to_csv(pipeline_path+"/Outputs_temp/"+data_ref_output_name+"/SoFA_calculation/"+dataset_name+"_stripped.tsv", sep = "\t", index = None)
    
    else:
        dataset_full = pd.read_csv(pipeline_path+"/Inputs/"+dataset_name+".txt", index_col = 0, sep = "\t")
        esmecata_input = None
        otu_table_stripped = dataset_full

    deepmicro_otu = data_to_deepmicro(otu_table_stripped)
    
    return dataset_full, otu_table_stripped, esmecata_input, deepmicro_otu, label_file

def averaging_and_info_step(count_db, score_db, label_refs, esmecata_input, pipeline_path, dataset_name, args_passed):
    '''
    This script runs the steps to build an info database on the taxons and annotations from the dataset. Said databases will be used to add information to the final outputs.
    '''
    
    if not args_passed.annotations_only:
        info_taxons = get_info_taxons(count_db, esmecata_input, pipeline_path, dataset_name)
        avg_count_total = average_count_per_group(count_db, label_refs)

    info_annots = get_info_annots(score_db, esmecata_input, pipeline_path, dataset_name, args_passed)
    avg_score_total = average_per_group(score_db, label_refs)
    

    df_max_scores = pd.DataFrame()
    df_max_counts = pd.DataFrame()

    for label in np.unique(label_refs.values):
        
        if not args_passed.annotations_only:
            avg_count_group = average_count_per_group(count_db, label_refs, label)
            df_max_counts[label] = avg_count_group['average']
            info_taxons['Average presence in '+str(label)] = avg_count_group['average'].values



        avg_score_group = average_per_group(score_db, label_refs, label)
        df_max_scores[label] = avg_score_group['average']
        info_annots['Average presence in '+str(label)] = avg_score_group['average'].values
        
    
    if not args_passed.annotations_only:
        info_taxons['Representative_group'] = df_max_counts.idxmax(axis=1).values
        info_taxons['Average presence (total)'] = avg_count_total['average'].values

    else:
        info_taxons = None

    info_annots['Representative_group'] = df_max_scores.idxmax(axis=1).values
    info_annots['Average presence (total)'] = avg_score_total['average'].values
    
    return info_annots, info_taxons

def run_iterate(run_nb, dataset_name, data_ref, data_ref_output_name, iterations, pipeline_path, dataset_full, label_file, otu_table_stripped, esmecata_input, deepmicro_otu, sofa_table, deepmicro_sofa, info_annots, info_taxons, test_set_dict, bank_of_selections_annots, bank_of_selections_taxons, bank_of_performance_dfs_annots, bank_of_performance_dfs_taxons, args_passed):
    
    
    if args_passed.reference_test_sets:
        #Get the test set references if they are given
        test_set_refs = pd.read_csv(pipeline_path+'/Inputs/Test_sets_'+dataset_name+'.csv')
        test_labels = test_set_refs['Run_'+str(run_nb)].values

        sample_names = list(sofa_table.columns.values)
        test_indices = [sample_names.index(i) for i in test_labels]

        labels_test = label_file.iloc[test_indices].values.reshape((label_file.iloc[test_indices].values.shape[0]))
        #logger.info('Label_test:', labels_test)
        train_indices = []
        for i in range(len(label_file)):
            if i not in test_indices:
                train_indices.append(i)

        labels_train = label_file.iloc[train_indices].values.reshape((label_file.iloc[train_indices].values.shape[0]))

        

    else:
        #Select a test subset
        labels_train, labels_test, test_indices = separate_test_train(label_file)

        test_labels = [sofa_table.columns[i] for i in test_indices]

    #Keeping track of the selected test sets (and writing them at every run, in case of a crash)
    test_set_dict['Run_'+str(run_nb)] = test_labels
    test_set_df = pd.DataFrame.from_dict(test_set_dict)
    test_set_df.to_csv(pipeline_path+'/Meta-Outputs/'+data_ref_output_name+'/Test_sets.csv')

    if not args_passed.annotations_only:
        deepmicro_otu_iteration = deepmicro_otu
    deepmicro_sofa_iteration = deepmicro_sofa

    

    os.mkdir(pipeline_path+'/Meta-Outputs/'+data_ref_output_name+'/Run_'+str(run_nb)+'/Trained_classifiers')
    os.mkdir(pipeline_path+'/Meta-Outputs/'+data_ref_output_name+'/Run_'+str(run_nb)+'/Classification_performances')
    os.mkdir(pipeline_path+'/Meta-Outputs/'+data_ref_output_name+'/Run_'+str(run_nb)+'/Selected_Variables')

    #ITERATED:
    for iteration_number in range(iterations):
        
    #Separate test and train subsets
        if not args_passed.annotations_only:
            otu_train , otu_test = create_test_train_subsets(deepmicro_otu_iteration, test_labels)
        annots_train , annots_test = create_test_train_subsets(deepmicro_sofa_iteration, test_labels)
        

    # annots_train.to_csv(args.pipeline_path+'/DeepMicro/data/entree_DeepMicro_'+args.dataset_name+'.csv', header = None, index = None)
    # annots_test.to_csv(args.pipeline_path+'/DeepMicro/data/entree_DeepMicro_'+args.dataset_name+'_test.csv', header = None, index = None)


    # otu_train.to_csv(args.pipeline_path+'/DeepMicro/data/entree_DeepMicro_'+args.dataset_name+'_OTU.csv', header = None, index = None)
    # otu_test.to_csv(args.pipeline_path+'/DeepMicro/data/entree_DeepMicro_'+args.dataset_name+'_OTU_test.csv', header = None, index = None)

    # labels_train.to_csv(args.pipeline_path+'/DeepMicro/data/Label_'+dataset_no_iter+'.csv', header = None, index = None)
    # labels_test.to_csv(args.pipeline_path+'/DeepMicro/data/Label_'+dataset_no_iter+'_test.csv', header = None, index = None)

        #DeepMicro
        os.mkdir(pipeline_path+'/Meta-Outputs/'+data_ref_output_name+'/Run_'+str(run_nb)+'/Trained_classifiers/Iteration_'+str(iteration_number))
        data_dir = pipeline_path+'/Meta-Outputs/'+data_ref_output_name+'/Run_'+str(run_nb)+'/Trained_classifiers/Iteration_'+str(iteration_number)
        if not args_passed.annotations_only:
            perf_df_otu, best_feature_records_otu = run_deep_micro(otu_test, otu_train, labels_test, labels_train, dataset_name+'_OTU', iteration_number, run_nb, pipeline_path, args_passed, data_dir, "Taxonomic")
        perf_df_sofa, best_feature_records_sofa = run_deep_micro(annots_test, annots_train, labels_test, labels_train, dataset_name+'_Functions', iteration_number, run_nb, pipeline_path, args_passed, data_dir, "Functional")
        

        os.mkdir(pipeline_path+'/Meta-Outputs/'+data_ref_output_name+'/Run_'+str(run_nb)+'/Classification_performances/Iteration_'+str(iteration_number))
        
        if not args_passed.annotations_only:
            perf_df_otu.to_csv(pipeline_path+'/Meta-Outputs/'+data_ref_output_name+'/Run_'+str(run_nb)+'/Classification_performances/Iteration_'+str(iteration_number)+'/Taxonomic_performances.csv')
        perf_df_sofa.to_csv(pipeline_path+'/Meta-Outputs/'+data_ref_output_name+'/Run_'+str(run_nb)+'/Classification_performances/Iteration_'+str(iteration_number)+'/Annotation_performances.csv')
        
        if not args_passed.annotations_only:
            bank_of_performance_dfs_taxons[iteration_number][run_nb] = perf_df_otu
        bank_of_performance_dfs_annots[iteration_number][run_nb] = perf_df_sofa

        if not args_passed.annotations_only:
            best_feature_records_otu_df = pd.concat(best_feature_records_otu, axis=1)
            best_feature_records_otu_df.columns = ['Model_'+str(i) for i in range(int(args_passed.forests))]
            best_feature_records_otu_df.index = deepmicro_otu_iteration.columns
            #best_feature_records_otu_df.to_csv(pipeline_path+'/Meta-Outputs/'+data_ref_output_name+'/Run_'+str(run_nb)+'/Classification_performances/Iteration_'+str(iteration_number)+'/Best_feature_records_taxons.csv')
            
        best_feature_records_sofa_df = pd.concat(best_feature_records_sofa, axis=1)
        best_feature_records_sofa_df.columns = ['Model_'+str(i) for i in range(int(args_passed.forests))]
        best_feature_records_sofa_df.index = deepmicro_sofa_iteration.columns
        #best_feature_records_sofa_df.to_csv(pipeline_path+'/Meta-Outputs/'+data_ref_output_name+'/Run_'+str(run_nb)+'/Classification_performances/Iteration_'+str(iteration_number)+'/Best_feature_records_annotations.csv')
        
        #Average Gini importances
        if not args_passed.annotations_only:
            best_feature_records_otu_df['Average'] = best_feature_records_otu_df.mean(axis=1)
        best_feature_records_sofa_df['Average'] = best_feature_records_sofa_df.mean(axis=1)
        
        #Rotor cutoff

        os.mkdir(pipeline_path+'/Meta-Outputs/'+data_ref_output_name+'/Run_'+str(run_nb)+'/Selected_Variables/Iteration_'+str(iteration_number))

        if not args_passed.annotations_only:
            retained_otus = inflexion_cutoff(best_feature_records_otu_df)
        retained_annots = inflexion_cutoff(best_feature_records_sofa_df)

        if not args_passed.annotations_only:
            selection_plus_info_taxons = info_taxons[info_taxons['ID'].isin(list(retained_otus.index))]
        selection_plus_info_annots = info_annots[info_annots['ID'].isin(list(retained_annots.index))]
        
        
        # Write the selection files with info
        if not args_passed.annotations_only:
            signif_otus = []
            signif_otu_names = []

            for link_otu_list in selection_plus_info_annots['Linked_OTUs'].values:
                signif_links = []
                signif_links_named = []
                for otu in link_otu_list:
                    if otu in retained_otus:
                        signif_links.append(otu)
                        otu_name_translated = esmecata_input[esmecata_input['observation_name'] == otu]['taxonomic_affiliation'].values[0]
                        otu_name_translated_species = otu_name_translated.split(';')[-1]
                        signif_links_named.append(otu_name_translated_species)
                
                signif_otus.append(signif_links)
                signif_otu_names.append(signif_links_named)

            selection_plus_info_annots['Significant_linked_OTUs'] = signif_otus
            selection_plus_info_annots['Significant_linked_Named_OTUs'] = signif_otu_names


            signif_annots = []

            for link_annot_list in selection_plus_info_taxons['Linked_annotations'].values:
                signif_links = []

                for annot in link_annot_list:
                    if annot in retained_annots:
                        signif_links.append(annot)
                
                signif_annots.append(signif_links)


            selection_plus_info_taxons['Significant_linked_annotations'] = signif_annots
            selection_plus_info_taxons.to_csv(pipeline_path+'/Meta-Outputs/'+data_ref_output_name+'/Run_'+str(run_nb)+'/Selected_Variables/Iteration_'+str(iteration_number)+'/Selected_taxons_run_'+str(run_nb)+'_iter_'+str(iteration_number)+'.csv')

        selection_plus_info_annots.to_csv(pipeline_path+'/Meta-Outputs/'+data_ref_output_name+'/Run_'+str(run_nb)+'/Selected_Variables/Iteration_'+str(iteration_number)+'/Selected_annotations_run_'+str(run_nb)+'_iter_'+str(iteration_number)+'.csv')
        


        # Store the selections for later incorporation to the Core/Meta sets
        if not args_passed.annotations_only:
            bank_of_selections_taxons[iteration_number][run_nb] = list(retained_otus.index)
        bank_of_selections_annots[iteration_number][run_nb] = list(retained_annots.index)
        
       
        if not args_passed.annotations_only:
            deepmicro_otu_iteration = otu_table_stripped.loc[retained_otus.index].transpose()
        deepmicro_sofa_iteration = sofa_table.loc[retained_annots.index].transpose()


    return(test_set_dict, bank_of_selections_annots, bank_of_selections_taxons, bank_of_performance_dfs_annots, bank_of_performance_dfs_taxons)

### PLOTTING RESULTS

def get_median_perfs_and_best_iter(bank_of_performance_dfs, median_classifs_per_iteration):
    '''
    This function calculates the median performances of each run for a given selection iteration, and returns the ones corresponding to the best iteration (i.e: best mean of median AUC performances)
    '''

    best_mean_perf = 0

    # Processing the mean performance of each run per level of iteration
    for iteration_lv in bank_of_performance_dfs.keys():
        iteration_lv_perfs = bank_of_performance_dfs[iteration_lv]

        for run_lv in iteration_lv_perfs.keys():
            run_df = iteration_lv_perfs[run_lv]
            test_perf_values = run_df['Test performance'].values
            median_test_perf_values = np.median(test_perf_values)
            median_classifs_per_iteration[iteration_lv].append(median_test_perf_values)
        
        #Finding and recording the best performing iteration
        mean_perf = np.mean(median_classifs_per_iteration[iteration_lv])
        if mean_perf > best_mean_perf:
            best_mean_perf = mean_perf
            best_selec_iter = iteration_lv
    
    return median_classifs_per_iteration[best_selec_iter], mean_perf, best_selec_iter

def plot_classifs(bank_of_performance_dfs_annots, bank_of_performance_dfs_taxons, dataset_name, pipeline_path, data_ref_output_name, args_passed):
    '''
    This function processes the classification performances obtained through the previous iterative process, and plots the best iteration's performances
    '''

    # Gather the best iteration and the corresponding median performances for each iterative level
    median_classifs_best_iteration_annots = defaultdict(list)
    median_classifs_best_iteration_annots, mean_perf_annots, best_selec_iter_annots = get_median_perfs_and_best_iter(bank_of_performance_dfs_annots, median_classifs_best_iteration_annots)
    best_iteration_indices = [best_selec_iter_annots]


    if not args_passed.annotations_only:
        median_classifs_best_iteration_taxons = defaultdict(list)
        median_classifs_best_iteration_taxons, mean_perf_taxons, best_selec_iter_taxons = get_median_perfs_and_best_iter(bank_of_performance_dfs_taxons, median_classifs_best_iteration_taxons)
        best_iteration_indices.append(best_selec_iter_taxons)
        uval, pval = mannwhitneyu(median_classifs_best_iteration_annots, median_classifs_best_iteration_taxons)  
    else:
        median_classifs_best_iteration_taxons = []
        best_selec_iter_taxons = None
    
    # Create a dataframe to be plotted
    performance_scores = median_classifs_best_iteration_annots + median_classifs_best_iteration_taxons
    profiles = ['Functional' for i in range(len(median_classifs_best_iteration_annots))]
    data_name = [dataset_name  for i in range(len(median_classifs_best_iteration_annots))]
    best_run = [best_selec_iter_annots for i in range(len(median_classifs_best_iteration_annots))]
    mean_perfs= [mean_perf_annots for i in range(len(median_classifs_best_iteration_annots))]


    if not args_passed.annotations_only:
        profiles = profiles + ['Taxonomic' for i in range(len(median_classifs_best_iteration_taxons))]
        data_name = data_name + [dataset_name  for i in range(len(median_classifs_best_iteration_taxons))]
        best_run = best_run + [best_selec_iter_taxons for i in range(len(median_classifs_best_iteration_taxons))]
        mean_perfs= mean_perfs + [mean_perf_taxons for i in range(len(median_classifs_best_iteration_taxons))]

    dataframe_to_plot = pd.DataFrame.from_dict({'AUC score': performance_scores, 'Profile': profiles, 'Disease': data_name, 'Best run': best_run, 'Mean Perfs': mean_perfs})

    # Plot the info
    fig,ax = plt.subplots(1,constrained_layout=True, figsize=(10,10))

    medians = sns.boxplot(
            y='Disease',
            x='AUC score',
            hue = 'Profile',
            data=dataframe_to_plot,
            ax=ax,
            palette = 'pastel',
            fliersize = 0)
    
    sns.stripplot(y='Disease', x='AUC score', hue = 'Profile', data=dataframe_to_plot, jitter = True, dodge=True, ax=ax, palette = 'dark')
    

    # Plot presentation
    max_val = max(performance_scores)

    plt.text(max_val + 0.01, -0.5, 'Optimal number\nof selections', size=15, fontweight = 'bold')
    
    for xi,yi,tx, col in zip([max_val + 0.03 for i in range(len(np.unique(profiles)))], (np.arange(len(dataframe_to_plot['Mean Perfs'].values))/2)-0.2, best_iteration_indices, ['blue','darkorange']):
        plt.text( xi, yi,tx, size=18, fontweight = 'bold', c=col)
    
    handles, labels = ax.get_legend_handles_labels()
    l = plt.legend(loc='upper left', borderaxespad=0., fontsize = 15)

    plt.yticks(size=15, rotation=90)
    ax.set_ylabel('Disease',size=20)

    plt.xticks(size=15)
    ax.set_xlabel('Median AUC',size=20)

    title = 'Classification_performances'
    if not args_passed.annotations_only:
        title += '\n(Mann-Whitney test p-value: '+str(pval)+')'

    ax.set_title(title)
    plt.savefig(pipeline_path+'/Meta-Outputs/'+data_ref_output_name+'/median_OTU_vs_SoFA_(best_vs_best).png', dpi=300)

    return(best_selec_iter_annots, best_selec_iter_taxons)

### GATHER CORE AND META FEATURES

def extract_core_associates(dataframe, core_list, esmecata_input=None):
    '''
    This function automatically extracts the core-significant from the lists of associated variables
    '''

    if 'Linked_annotations' in dataframe.columns:
        col_ref = 'Linked_annotations'
        new_col = 'Core_linked_annotaions'
        otu_links = False
    else:
        col_ref = 'Linked_OTUs'
        new_col = 'Significant_linked_OTUs'
        otu_links = True

    signif_vars = []
    for vars_assoc in dataframe[col_ref].values:
        vars_list = []
        for var in vars_assoc:
            if var in core_list:
                vars_list.append(var)
        signif_vars.append(vars_list)
    
    dataframe[new_col] = signif_vars

    if otu_links:
        signif_links_named = []
        for otu_list in signif_vars:
            named_links = []
            for otu in otu_list:
                otu_name_translated = esmecata_input[esmecata_input['observation_name'] == otu]['taxonomic_affiliation'].values[0]
                otu_name_translated_species = otu_name_translated.split(';')[-1]
                named_links.append(otu_name_translated_species)
            signif_links_named.append(named_links)
        dataframe['Named_significant_linked_OTUs'] = signif_links_named

    return dataframe


def create_core_and_meta_dfs(iteration_selections_per_run, runs):
    '''
    This function creates, from a dictionary containing selected variables per iteration level and run:
        - a 'Core' dataframe containing the features that are present in all of the selections for a given iteration
        - a 'Meta' dataframe containing the features that are present in at least one of the selections for a given iteration, along with the number of runs in which it is selected
    '''
    
    full_list = sum([list(iteration_selections_per_run[i]) for i in iteration_selections_per_run.keys()],[])
    core_meta = {'core':[], 'meta':{'ID':[],'Count':[]}}

    for feature in np.unique(full_list):
        if full_list.count(feature) == runs:
            core_meta['core'].append(feature)
        else:
            core_meta['meta']['ID'].append(feature)
            core_meta['meta']['Count'].append(full_list.count(feature))
    
    core = pd.DataFrame(core_meta['core'], columns=['ID'])
    meta = pd.DataFrame(core_meta['meta']).sort_values(by='Count', ascending=False)

    return(core,meta)



def formatting_core_meta_outputs(info_df, core_df, meta_df, runs, zero_case=False):
    '''
    This function adds information inherited from upstream analyses to the core and meta dataframes
    '''
    
    if meta_df is None:
        meta_skip = True
    else:
        meta_skip = False
    if zero_case:
        core_info = info_df[info_df['ID'].isin(list(core_df.index))]
    else:
        core_info = info_df[info_df['ID'].isin(list(core_df['ID'].values))]

    if meta_skip:
        meta_info = None
        
    
    else:
        significance_category = []
        for count in meta_df['Count'].values:
            if count > runs*0.75:
                significance_category.append('Confident')
            else:
                significance_category.append('Candidate')

        meta_info = info_df[info_df['ID'].isin(list(meta_df['ID'].values))]
        meta_info['Significance_count'] = meta_df['Count'].values
        meta_info['Significance_category'] = significance_category

    return core_info, meta_info

def extract_and_write_core_meta(path_core_meta, bank_of_selections_annots, bank_of_selections_taxons, best_selec_iter_annots, best_selec_iter_taxons, info_annots, info_taxons, runs, esmecata_input, otu_table_stripped, sofa_table, args_passed):
    '''
    This function launches the process to identify and write the Core and Meta taxons and annotations from all of the effectuated selections.
    '''

    for iteration in bank_of_selections_annots.keys():
        
        iteration_selections_per_run_annots = bank_of_selections_annots[iteration]
        core_annots, meta_annots = create_core_and_meta_dfs(iteration_selections_per_run_annots, runs)
        core_annot_info, meta_annot_info = formatting_core_meta_outputs(info_annots, core_annots, meta_annots, runs, zero_case=False)

        if iteration == best_selec_iter_annots - 1:
            core_annots_opti, meta_annots_opti = core_annot_info, meta_annot_info


        if not args_passed.annotations_only:
            iteration_selections_per_run_taxons = bank_of_selections_taxons[iteration]
            core_taxons, meta_taxons = create_core_and_meta_dfs(iteration_selections_per_run_taxons, runs)
            core_taxons_info, meta_taxons_info = formatting_core_meta_outputs(info_taxons, core_taxons, meta_taxons, runs, zero_case=False)

            if iteration == best_selec_iter_taxons -1:
                core_taxons_opti, meta_taxons_opti = core_taxons_info, meta_taxons_info
            
            core_taxons_id = list(core_taxons['ID'].values)
            core_annots_id = list(core_annots['ID'].values)
            
            core_annot_info = extract_core_associates(core_annot_info, core_taxons_id, esmecata_input)
            core_taxons_info = extract_core_associates(core_taxons_info, core_annots_id, esmecata_input)

            meta_annot_info = extract_core_associates(meta_annot_info, core_taxons_id, esmecata_input)
            meta_taxons_info = extract_core_associates(meta_taxons_info, core_annots_id, esmecata_input)

            core_taxons_info.to_csv(path_core_meta+'/All_iterations/Core_taxons_iteration_'+str(iteration)+'.csv')
            meta_taxons_info.to_csv(path_core_meta+'/All_iterations/Meta_taxons_iteration_'+str(iteration)+'.csv')
        
        core_annot_info.to_csv(path_core_meta+'/All_iterations/Core_annots_iteration_'+str(iteration)+'.csv')
        meta_annot_info.to_csv(path_core_meta+'/All_iterations/Meta_annots_iteration_'+str(iteration)+'.csv')

    if best_selec_iter_annots == 0:
        core_annots_opti, meta_annots_opti = formatting_core_meta_outputs(info_annots, sofa_table, None, runs, zero_case=True)

    if not args_passed.annotations_only:
        
        core_taxons_opti_id = list(core_taxons_opti['ID'].values)
        core_annots_opti_id = list(core_annots_opti['ID'].values)

        if best_selec_iter_taxons == 0:
            core_taxons_opti, meta_taxons_info = formatting_core_meta_outputs(info_taxons, otu_table_stripped, None, runs, zero_case=True)

        else:
            meta_taxons_opti = extract_core_associates(meta_taxons_opti, core_annots_opti_id, esmecata_input)
        
        if best_selec_iter_annots != 0:
            meta_annots_opti = extract_core_associates(meta_annots_opti, core_taxons_opti_id, esmecata_input)
        
        
        core_annots_opti = extract_core_associates(core_annots_opti, core_taxons_opti_id, esmecata_input)
        core_taxons_opti = extract_core_associates(core_taxons_opti, core_annots_opti_id, esmecata_input)

        

        core_taxons_opti.to_csv(path_core_meta+'/Best_iteration/Core_taxons_iteration_'+str(best_selec_iter_taxons-1)+'.csv')
        meta_taxons_opti.to_csv(path_core_meta+'/Best_iteration/Meta_taxons_iteration_'+str(best_selec_iter_taxons-1)+'.csv')

    core_annots_opti.to_csv(path_core_meta+'/Best_iteration/Core_annots_iteration_'+str(best_selec_iter_annots-1)+'.csv')
    meta_annots_opti.to_csv(path_core_meta+'/Best_iteration/Meta_annots_iteration_'+str(best_selec_iter_annots-1)+'.csv')
        
        



def main(args_passed):
    '''
    This function formats the inputs, formats the output directories, and launches the iterative seection and classification process r times.
    '''
    
    now_begin = datetime.now()
    date_time = now_begin.strftime("%d%m%Y%H%M")


    dataset_name = args_passed.dataset_name
    pipeline_path = os.getcwd()
    runs = int(args_passed.runs)
    iterations = int(args_passed.iterations)

    ## Checking the validity of the args
    if iterations <1:
        raise ValueError("Please input a valid amount of iterations")
    if runs <1:
        raise ValueError("Please input a valid amount of pipeline runs")

    if not (os.path.isfile(pipeline_path+'/Inputs/'+dataset_name+'.txt') and os.path.isfile(pipeline_path+'/Inputs/Label_'+dataset_name+'.csv')):
        raise ValueError("Please give an input for "+dataset_name+".txt and Label_"+dataset_name+".csv")

    if not args_passed.treatment in ["tf_igm", None]:
        raise ValueError("Please enter a valid value for the '-t' argument: 'scaling', 'tf_igm' or None")

    if not args_passed.scaling in ["relative", None]:
        raise ValueError("Please enter a valid value for the '-s' argument: 'relative' or None")
    
    if args_passed.reference_test_sets:
        if not os.path.isfile(pipeline_path+'/Inputs/Test_sets_'+dataset_name+'.csv'):
            raise ValueError("If using the --reference_test_sets flag, please give an input for Test_sets"+dataset_name+".csv")


    ## Setting up the output directories
    data_ref = dataset_name

    if args_passed.treatment != None:
        data_ref = dataset_name+'_'+args_passed.treatment

    if args_passed.scaling != None:
        data_ref = data_ref+'_'+args_passed.scaling
    
    if args_passed.annotations_only:
        data_ref = data_ref+'_annotations_only'

    data_ref_output_name = data_ref+'_'+date_time
   

    #Temporary output directory
    if not os.path.isdir(pipeline_path+'/Outputs_temp/'):
        os.mkdir(pipeline_path+'/Outputs_temp/')
        logger.info('Created Output directory')
    
    if not os.path.isdir(pipeline_path+'/Outputs_temp/'+data_ref_output_name):
        os.mkdir(pipeline_path+'/Outputs_temp/'+data_ref_output_name)
    else:
        raise Exception("Another instance of SPARTA was launched on the same dataset too soon before this one. The names of our output directories include the current date, to the minute, in order to be unique. Right now, there will be a conflict in directory names. Please try again in a minute!")

    #Permanent output directory
    if not os.path.isdir(pipeline_path+'/Meta-Outputs'):
        os.mkdir(pipeline_path+'/Meta-Outputs')
        logger.info('Created Meta-Output directory')

    os.mkdir(pipeline_path+'/Meta-Outputs/'+data_ref_output_name)
    logger.info('Meta-Outputs ready! ID: '+data_ref_output_name)

    

    #Temporary necessities for calculation of SoFAs
    if not os.path.isdir(pipeline_path+'/Outputs_temp/'+data_ref_output_name+'/SoFA_calculation'):
        os.mkdir(pipeline_path+'/Outputs_temp/'+data_ref_output_name+'/SoFA_calculation')
        logger.info('Created Temporary Sofa calculation directory')
    

    ## Formatting step
    if not args_passed.annotations_only:
        logger.info('Input formatting...')
        dataset_full, otu_table_stripped, esmecata_input, deepmicro_otu, label_file = formatting_step(dataset_name, data_ref_output_name, pipeline_path, args_passed)

        ####Time measurement####
        date_time_format = datetime.now()
        formatting_time = date_time_format - now_begin
        formatting_time_seconds = formatting_time.total_seconds()

        stopwatch_file = pipeline_path+"/stopwatch_"+data_ref_output_name+".txt"
        f = open(stopwatch_file, "w")
        f.write("Formatting step length (s): "+str(formatting_time_seconds)+"\n")
        f.close()
        ########################

        #EsMeCaTa output folder
        if not os.path.isdir(pipeline_path+'/EsMeCaTa_outputs'):
            os.mkdir(pipeline_path+'/EsMeCaTa_outputs')

        if os.path.isdir(pipeline_path+'/EsMeCaTa_outputs/'+dataset_name):
            if not args_passed.esmecata_relaunch:
                logger.info('An EsMeCaTa output has been found for your dataset. This output will be used for the rest of the pipeline. If you wish to re-launch EsMeCaTa, please remove the existing output before launching SPARTA.')
            else:
                logger.info('Re-launching EsMeCaTa over previous results')
                if os.path.isdir(pipeline_path+'/EsMeCaTa_outputs/'+dataset_name+'/esmecata_outputs_annots/annotation_reference'):
                    #Removing previous 'annotations_reference' output to ensure the final step of EsMeCaTa does not malfunction
                    shutil.rmtree(pipeline_path+'/EsMeCaTa_outputs/'+dataset_name+'/esmecata_outputs_annots/annotation_reference', ignore_errors=True)
                esmecata_plus_check(esmecata_input, pipeline_path, data_ref_output_name, dataset_name, args_passed.eggnog)

        else:
            os.mkdir(pipeline_path+'/EsMeCaTa_outputs/'+dataset_name)
            logger.info('Launching EsMeCaTa')
            esmecata_plus_check(esmecata_input, pipeline_path, data_ref_output_name, dataset_name, args_passed.eggnog)

        ####Time measurement####
        date_time_esmecata = datetime.now()
        esmecata_time = date_time_esmecata - date_time_format
        esmecata_time_seconds = esmecata_time.total_seconds()

        f = open(stopwatch_file, "a")
        f.write("EsMeCaTa step length (s): "+str(esmecata_time_seconds)+"\n")
        f.close()
        ########################

        ## Calculating the scores of functional annotations
        sofa_table, deepmicro_sofa = sofa_calculation(pipeline_path, dataset_name, data_ref_output_name, otu_table_stripped, args_passed)

        ####Time measurement####
        date_time_score = datetime.now()
        score_calculation_time = date_time_score - date_time_esmecata
        score_calculation_time_seconds = score_calculation_time.total_seconds()

        f = open(stopwatch_file, "a")
        f.write("Functional score calculation step length (s): "+str(score_calculation_time_seconds)+"\n")
        f.close()
        ########################

        ref_time = date_time_score
    
    else:
        dataset_full, sofa_table, esmecata_input, deepmicro_sofa, label_file = formatting_step(dataset_name, data_ref_output_name, pipeline_path, args_passed)

        ####Time measurement####
        date_time_format = datetime.now()
        formatting_time = date_time_format - now_begin
        formatting_time_seconds = formatting_time.total_seconds()

        stopwatch_file = pipeline_path+"/stopwatch_"+data_ref_output_name+".txt"
        f = open(stopwatch_file, "w")
        f.write("Formatting step length (s): "+str(formatting_time_seconds)+"\n")
        f.close()
        ########################

        otu_table_stripped, deepmicro_otu = (None, None)

        ref_time = date_time_format


    ## Calculating average presence of taxons and annotations per label, and collecting info about them

    info_annots, info_taxons = averaging_and_info_step(otu_table_stripped, sofa_table, label_file, esmecata_input, pipeline_path, dataset_name, args_passed)

    #info_annots.to_csv(pipeline_path+'/Meta-Outputs/'+data_ref_output_name+'/annot_info_check.csv')

    # if not args_passed.annotations_only:
    #     info_taxons.to_csv(pipeline_path+'/Meta-Outputs/'+data_ref_output_name+'/taxon_info_check.csv')

    ##A dictionary to keep track of the test sets and selection/performance outputs of the different runs

    test_set_dict = {}
    bank_of_selections_annots = defaultdict(defaultdict)
    bank_of_selections_taxons = defaultdict(defaultdict)
    bank_of_performance_dfs_annots = defaultdict(defaultdict)
    bank_of_performance_dfs_taxons = defaultdict(defaultdict)

    ## Launching the SPARTA runs
    for run_nb in range(1,runs+1):
        
        os.mkdir(pipeline_path+'/Meta-Outputs/'+data_ref_output_name+'/Run_'+str(run_nb))
        test_set_dict, bank_of_selections_annots, bank_of_selections_taxons, bank_of_performance_dfs_annots, bank_of_performance_dfs_taxons = run_iterate(run_nb, dataset_name, data_ref, data_ref_output_name, iterations, pipeline_path, dataset_full, label_file, otu_table_stripped, esmecata_input, deepmicro_otu, sofa_table, deepmicro_sofa, info_annots, info_taxons, test_set_dict, bank_of_selections_annots, bank_of_selections_taxons, bank_of_performance_dfs_annots, bank_of_performance_dfs_taxons, args_passed)

        ####Time measurement####
        date_time_now = datetime.now()
        run_i_time = date_time_now - ref_time
        run_i_time_seconds = run_i_time.total_seconds()

        f = open(stopwatch_file, "a")
        f.write("SPARTA Run "+str(run_nb)+" length (s): "+str(run_i_time_seconds)+"\n")
        f.close()

        ref_time = date_time_now
        ########################

    best_selec_iter_annots, best_selec_iter_taxons = plot_classifs(bank_of_performance_dfs_annots, bank_of_performance_dfs_taxons, dataset_name, pipeline_path, data_ref_output_name, args_passed)
    
    os.mkdir(pipeline_path+'/Meta-Outputs/'+data_ref_output_name+'/Core_and_Meta_outputs')
    os.mkdir(pipeline_path+'/Meta-Outputs/'+data_ref_output_name+'/Core_and_Meta_outputs/All_iterations')
    os.mkdir(pipeline_path+'/Meta-Outputs/'+data_ref_output_name+'/Core_and_Meta_outputs/Best_iteration')

    path_core_meta = pipeline_path+'/Meta-Outputs/'+data_ref_output_name+'/Core_and_Meta_outputs'
    
    extract_and_write_core_meta(path_core_meta, bank_of_selections_annots, bank_of_selections_taxons, best_selec_iter_annots, best_selec_iter_taxons, info_annots, info_taxons, runs, esmecata_input, otu_table_stripped, sofa_table, args_passed)

    if not args_passed.keep_temp:
        shutil.rmtree(pipeline_path+'/Outputs_temp/'+data_ref_output_name, ignore_errors=True)
    # pd.DataFrame.from_dict(bank_of_selections_annots).to_csv(pipeline_path+'/Meta-Outputs/'+data_ref_output_name+'/bank_of_selections_check.csv')

    ####Time measurement####
    end_time = datetime.now()

    post_process_time = end_time - ref_time
    post_process_time_seconds = post_process_time.total_seconds()

    full_process_time = end_time - now_begin
    full_process_time_seconds = full_process_time.total_seconds()

    f = open(stopwatch_file, "a")
    f.write("Post processing step length (s): "+str(post_process_time_seconds)+"\n")
    f.write("TOTAL: FULL PROCESS LENGTH (s): "+str(full_process_time_seconds) +"\n")

    ref_time = run_i_time
    ########################


try:
    main(args_passed)
except ValueError as error:
    logger.error(error)
    raise
