import pandas as pd
import os
import csv
import logging
import pandas as pd
import numpy as np
import os
import math
from collections import Counter
from tqdm import tqdm
import shutil
from datetime import datetime

import multiprocessing


from esmecata.proteomes import retrieve_proteomes
from esmecata.clustering import make_clustering
from esmecata.annotation import annotate_proteins
from esmecata.eggnog import annotate_with_eggnog
from esmecata.utils import get_rest_uniprot_release, get_sparql_uniprot_release, is_valid_file, is_valid_dir, send_uniprot_sparql_query

from ete3 import NCBITaxa


logger = logging.getLogger(__name__)
logging.getLogger("esmecata").setLevel(logging.DEBUG)


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

def pre_formatting(dataset_full, output_file):
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
    formatted_data.to_csv(output_file, sep = "\t", index = None)

    return formatted_data, dataset_compos

def formatting_step(label_filepath, abundance_file, output_folder, annotations_only=None, scaling='no scaling'):
    '''
    This script runs the steps to format the data into:
        - a metadata-less version of the original OTU table (otu_table_stripped)
        - an input to launch the esmecata pipeline afterwards (esmecata_input)
        - a taxonomic input for DeepMicro (deepmicro_otu)
    
        If --annotation only: the functional dataset will be called 'otu_table' in context of this function, and will be imported and converted to deepmicro input
    '''
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    label_file = pd.read_csv(label_filepath)

    if not annotations_only:
        dataset_full = pd.read_csv(abundance_file, header = None, sep = "\t")
        esmecata_input_path = os.path.join(output_folder, 'sofa_calculation.tsv')
        esmecata_input, otu_table_stripped = pre_formatting(dataset_full, esmecata_input_path)
        if scaling == "relative":
            otu_table_stripped = absolute_to_relative(otu_table_stripped)
        #Writing standardised version of the dataset
        otu_table_stripped_filepath = os.path.join(output_folder, 'otu_table_stripped.tsv')
        otu_table_stripped.to_csv(otu_table_stripped_filepath, sep = "\t")
    
    else:
        dataset_full = pd.read_csv(abundance_file, index_col = 0, sep = "\t")
        esmecata_input = None
        esmecata_input_path = None
        otu_table_stripped = dataset_full

    label_file_df = label_file[otu_table_stripped.columns].transpose()

    deepmicro_otu = data_to_deepmicro(otu_table_stripped)

    return dataset_full, otu_table_stripped, esmecata_input, esmecata_input_path, deepmicro_otu, label_file_df


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

def sofa_calculation(esmecata_annotation_reference, output_sofa_table_filepath, otu_table_stripped, treatment='tf_igm'):
    '''
    This function calculates the scores of the functional annotations from EsMeCaTa's outputs and the original OTU abundances.

    INPUTS: 
        - original args
        - OTU table
    
    OUTPUTS:
        - sofa_table: a table containing all of the scores of functional annotations
        - deepmicro_sofa: the same table, but formatted as a DeepMicro input. If a transformation is given as an argument, 
    '''

    esmecata_output_path = esmecata_annotation_reference
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

    sofa_table.to_csv(output_sofa_table_filepath)
    deepmicro_sofa = sofa_table.transpose()

    if treatment == 'tf_igm':
        deepmicro_sofa = tf_igm_apply(deepmicro_sofa)

    #deepmicro_sofa.to_csv(pipeline_path+'/Outputs_temp/'+data_ref_output_name+'/DeepMicro_data/entree_DeepMicro_'+dataset_name+'.csv', sep=',', header = None, index = None)

    return sofa_table, deepmicro_sofa


def proteome_check(esmecata_prot_out):
    '''
    This function checks if a proteome has been downloaded for each OTU given as input to EsMeCaTa
    '''

    path_to_proteome_count = os.path.join(esmecata_prot_out, 'proteome_tax_id.tsv')

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
    return all_proteomes_checked, count

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


def esmecata_plus_check(esmecata_input, esmecata_output_folder, eggnog_path=None, update_ncbi=None):
    '''
    This function runs EsMeCaTa on the input files. Checks are made at the end of the 'Proteomes' and 'Annotation' processes 
    '''
    nb_cpu_available = multiprocessing.cpu_count()
    if nb_cpu_available <= 3:
        nb_cpu_available = 1
    if nb_cpu_available > 3:
        nb_cpu_available -= 2

    if update_ncbi:
        logger.info('Updating local NCBI database')
        ncbi = NCBITaxa()
        ncbi.update_taxonomy_database()

    input_location = esmecata_input
    esmecata_prot_out = os.path.join(esmecata_output_folder, 'esmecata_outputs_proteomes')
    esmecata_cluster_out = os.path.join(esmecata_output_folder, 'esmecata_outputs_clustering')
    esmecata_annots_out = os.path.join(esmecata_output_folder, 'esmecata_outputs_annots')

    ## EsMeCaTa
    count_check = False
    retries = 0

    path_proteomes = os.path.join(esmecata_prot_out, 'proteomes')
    while (not count_check) and (retries<20):
        try:
            retrieve_proteomes(input_location, esmecata_prot_out, option_bioservices=False)
        except Exception as exception:
            logger.critical(exception)
            pass

        count_check, nb_proteomes = proteome_check(esmecata_prot_out)
        retries+=1
    
    if retries >= 20:
        raise Exception("EsMeCaTa has failed 20 times in a row. A connexion error is likely. Aborting...")

    #Temporary solution to the UniProt troubles of EsMeCaTa: if a downloaded proteome is less than 50 kB, remove it (considered empty)
    filenames = [[entry.name, entry.stat().st_size]  for entry in sorted(os.scandir(path_proteomes),
                                                key=lambda x: x.stat().st_size, reverse=False)]
    removed_proteomes = []
    for file_stats in filenames:
        if file_stats[1]<=50:
            filename = file_stats[0]
            os.remove(os.path.join(path_proteomes, filename))
            removed_proteomes.append(filename.replace('.faa.gz', ''))

    proteome_tax_id_path = os.path.join(esmecata_prot_out, 'proteome_tax_id.tsv')
    df_proteome = pd.read_csv(proteome_tax_id_path, sep='\t')

    # Remove empty proteomes.
    lambda_remove_proteomes = lambda x: ','.join([proteome for proteome in x.split(',') if proteome not in removed_proteomes])
    df_proteome['proteome'] = df_proteome['proteome'].apply(lambda_remove_proteomes)
    # Remvoe empty rows.
    row_to_drop = df_proteome[df_proteome['proteome'] == ''].index

    df_proteome = df_proteome.drop(row_to_drop)
    df_proteome.to_csv(proteome_tax_id_path, sep='\t', index=False)
    # End of temporary solution.

    stat_number_clustering_filepath = os.path.join(esmecata_cluster_out, 'stat_number_clustering.tsv')
    if not os.path.exists(stat_number_clustering_filepath):
        if os.path.exists(esmecata_cluster_out):
            logger.info('Previous incomplete iteration of the clustering step found: deleting and starting over')
            shutil.rmtree(esmecata_cluster_out, ignore_errors=True)
        make_clustering(esmecata_prot_out, esmecata_cluster_out, nb_cpu=nb_cpu_available, mmseqs_options=None, clust_threshold=0.5, linclust=None, remove_tmp=True)
    else:
        logger.info('Clustering step already done, moving to annotation.')

    try:
        if eggnog_path is None:
            annotate_proteins(esmecata_cluster_out, esmecata_annots_out, uniprot_sparql_endpoint=None,
                            propagate_annotation=1, uniref_annotation=None, expression_annotation=None, bioservices=True)
        else:
            annotate_with_eggnog(esmecata_cluster_out, esmecata_annots_out, eggnog_path, nb_cpu=nb_cpu_available)
    except:
        pass

    ## Check
    retries = 0
    annotation_reference_folder = os.path.join(esmecata_annots_out, 'annotation_reference')
    reference_proteins_consensus_fasta_folder = os.path.join(esmecata_cluster_out, 'reference_proteins_consensus_fasta')
    proteomes_tax_id = os.path.join(esmecata_cluster_out, 'proteome_tax_id.tsv')
    check_outputs, empty_ids = check_annotation(reference_proteins_consensus_fasta_folder, annotation_reference_folder, proteomes_tax_id)
    logger.info('All annotations found: '+str(check_outputs))
    while (not check_outputs) and (retries<20):
        try:
            if eggnog_path is None:
                annotate_proteins(esmecata_cluster_out, esmecata_annots_out, uniprot_sparql_endpoint=None,
                                propagate_annotation=1, uniref_annotation=None, expression_annotation=None, bioservices=True)
            else:
                annotate_with_eggnog(esmecata_cluster_out, esmecata_annots_out, eggnog_path, nb_cpu=10)
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


def create_dataset_annotation_file(annotation_reference_folder, dataset_annotation_file_path, content="all"):
    """ From annotation reference folder, creates a file resuming the number of ECs for each observation names.

    Args:
        annotation_reference_folder (str): path to annotation reference folder
        dataset_annotation_file_path (str): path to output dataset annotation file
        content (str): indicates which data to parse (default 'EC' other possible value is 'GO' or 'all')

    Returns:
        dataset_annotation (dict): annotation dict: observation_name as key and EC number as value
    """
    if content not in ['EC', 'GO', 'all']:
        raise ValueError("Wrong content. Authorized values are 'EC', 'GO' or 'all.")

    dataset_annotation = {}
    total_annotations = []
    for annotation_file in os.listdir(annotation_reference_folder):
        annotation_file_name = os.path.splitext(annotation_file)[0]
        annotation_file_path = os.path.join(annotation_reference_folder, annotation_file)

        annotations = []
        with open(annotation_file_path, 'r') as open_annotation_input_file_path:
            csvreader = csv.DictReader(open_annotation_input_file_path, delimiter='\t')
            for line in csvreader:
                if content in ['EC', 'GO']:
                    intermediary_annots = line[content].split(',')
                if content == 'all':
                    intermediary_annots = line['GO'].split(',')
                    intermediary_annots.extend(line['EC'].split(','))
                annotations.extend(intermediary_annots)
        annotations = [annot for annot in annotations if annot != '']
        total_annotations.extend(annotations)
        dataset_annotation[annotation_file_name] = annotations

    total_annotations = list(set(total_annotations))

    with open(dataset_annotation_file_path, 'w') as dataset_annotation_file:
        csvwriter = csv.writer(dataset_annotation_file, delimiter='\t')
        csvwriter.writerow(['observation_name'] + total_annotations)
        for observation_name in dataset_annotation:
            occurrence_annotations = Counter(dataset_annotation[observation_name])
            observation_name_annotations = [occurrence_annotations[annot] for annot in total_annotations]
            csvwriter.writerow([observation_name] + observation_name_annotations)

    return dataset_annotation


def run_esmecata(label_filepath, abundance_filepath, output_folder, scaling='no scaling', esmecata_relaunch=None, eggnog_path=None, update_ncbi=None):
    ## Formatting step
    date_time_format = datetime.now()
    now_begin = date_time_format
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    logger.info('Input formatting...')
    dataset_full, otu_table_stripped, esmecata_input, esmecata_input_path, deepmicro_otu, label_file = formatting_step(label_filepath, abundance_filepath, output_folder, scaling=scaling)

    ####Time measurement####
    date_time_format = datetime.now()
    formatting_time = date_time_format - now_begin
    formatting_time_seconds = formatting_time.total_seconds()

    stopwatch_file = os.path.join(output_folder, 'stopwatch_.txt')
    f = open(stopwatch_file, "w")
    f.write("Formatting step length (s): "+str(formatting_time_seconds)+"\n")
    f.close()
    ########################

    #EsMeCaTa output folder
    esmecata_output_folder = os.path.join(output_folder, 'EsMeCaTa_outputs')
    annotation_reference_folder = os.path.join(esmecata_output_folder, 'esmecata_outputs_annots', 'annotation_reference')

    if os.path.isdir(esmecata_output_folder):
        if esmecata_relaunch is None:
            logger.info('An EsMeCaTa output has been found for your dataset. This output will be used for the rest of the pipeline. If you wish to re-launch EsMeCaTa, please remove the existing output before launching SPARTA.')
        else:
            logger.info('Re-launching EsMeCaTa over previous results')
            if os.path.isdir(annotation_reference_folder):
                #Removing previous 'annotations_reference' output to ensure the final step of EsMeCaTa does not malfunction
                shutil.rmtree(annotation_reference_folder, ignore_errors=True)
            esmecata_plus_check(esmecata_input_path, esmecata_output_folder, eggnog_path, update_ncbi)

    else:
        os.mkdir(esmecata_output_folder)
        logger.info('Launching EsMeCaTa')
        esmecata_plus_check(esmecata_input_path, esmecata_output_folder, eggnog_path, update_ncbi)

    ####Time measurement####
    date_time_esmecata = datetime.now()
    esmecata_time = date_time_esmecata - date_time_format
    esmecata_time_seconds = esmecata_time.total_seconds()

    f = open(stopwatch_file, "a")
    f.write("EsMeCaTa step length (s): "+str(esmecata_time_seconds)+"\n")
    f.close()
    ########################

    # Compute occurrences of annotations in organism
    functional_occurrence_filepath = os.path.join(output_folder, 'functional_occurrence.tsv')
    create_dataset_annotation_file(annotation_reference_folder, functional_occurrence_filepath, content="all")

    ## Calculating the scores of functional annotations
    SoFA_table_filepath = os.path.join(output_folder, 'SoFA_table.csv')
    sofa_table, deepmicro_sofa = sofa_calculation(annotation_reference_folder, SoFA_table_filepath, otu_table_stripped)

    ####Time measurement####
    date_time_score = datetime.now()
    score_calculation_time = date_time_score - date_time_esmecata
    score_calculation_time_seconds = score_calculation_time.total_seconds()

    f = open(stopwatch_file, "a")
    f.write("Functional score calculation step length (s): "+str(score_calculation_time_seconds)+"\n")
    f.close()
    ########################


    return SoFA_table_filepath, esmecata_input_path, functional_occurrence_filepath, otu_table_stripped