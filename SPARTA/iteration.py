import pandas as pd
import os
import numpy as np
import requests
import random


from tqdm import tqdm
from kneebow.rotor import Rotor
from Bio.ExPASy import Enzyme
from goatools import obo_parser
from sklearn.model_selection import train_test_split
from SPARTA.Deepmicro import run_exp

def inflexion_cutoff(datatable):
    
    datatable_ordered = abs(datatable).sort_values(by="Average", ascending=False)

    data = datatable_ordered["Average"].values

    data_rotor = []
    for i in range(len(data)):
        data_rotor.append([i+1,data[i]])

    rotor = Rotor()
    rotor.fit_rotate(data_rotor)

    elbow_idx = rotor.get_elbow_index()

    datatable_truncated = datatable_ordered["Average"].iloc[:elbow_idx+1]

    return datatable_truncated


def separate_test_train(labels_full, seed_split):
    if not labels_full.values.shape[1] > 1:
        label_flatten = labels_full.values.reshape((labels_full.values.shape[0]))
    else:
        raise Exception("Please submit the label file in the form of an ordered vector")

    indices = np.arange(len(label_flatten))

    y_train, y_test, sample_train, sample = train_test_split(label_flatten, indices, test_size=0.2, stratify=label_flatten, random_state=seed_split)

    return y_train, y_test, sample_train, sample

def create_test_train_subsets(full_dataset, indices):
    '''
    This function separates a dataset into a test and a train subset, determined by the input indices
    '''
    dataset_train = full_dataset.drop(indices, axis=0)
    dataset_test = full_dataset.loc[indices]

    return(dataset_train, dataset_test)

def run_deep_micro(set_test, set_train, label_test, label_train, dataset_name, iteration_nb, run_nb, data_dir, profile, classifiers=20, method='rf', var_ranking_method='gini', real_seed = 0, seed_valid = 0):
    samples_labels_splits = {}
    try:
        ##IMPORT SEPARATE TEST DATASET
        Xtest_ext = set_test.values.astype(np.float64)

        Ytest_ext = label_test.astype(int)

        best_auc = 0
        threshold_opt = 0
        perf_dict = {}
        repeats = int(classifiers)
        best_feature_records = []
        
        random.seed(seed_valid)
        vec_split = random.sample(range(1000),repeats)
        
        random.seed(real_seed)
        seed_rf = random.sample(range(1000),repeats)
        # hyper-parameter grids for classifiers
        if repeats > 1:
            for i in tqdm(range(repeats), desc="Training models for the "+profile+" profile (Run: "+str(run_nb)+", Iteration: "+str(iteration_nb)+")"):
                #logger.info('BEST FEATURE RECORDS:', best_feature_records)
                best_auc, threshold_opt, perf_dict, best_feature_records, labels_sets = run_exp(i, best_auc, threshold_opt, perf_dict, Xtest_ext, Ytest_ext, set_train, label_train, dataset_name, best_feature_records, data_dir, method, var_ranking_method, seed_rf[repeats-1], vec_split[i])
                samples_labels_splits[i] = labels_sets
                samples_labels_splits[i]['test_set'] = ','.join(set_test.index.tolist())
                samples_labels_splits[i]['test_set_labels'] = ','.join(label_test.astype(str))
        else:
            best_auc, threshold_opt, perf_dict, best_feature_records, labels_sets = run_exp(42, best_auc, threshold_opt, perf_dict, Xtest_ext, Ytest_ext, set_train, label_train, dataset_name, best_feature_records, data_dir, method, var_ranking_method, real_seed = 42, seed_valid = 84)
            samples_labels_splits[repeats] = labels_sets
            samples_labels_splits[repeats]['test_set'] = ','.join(set_test.index.tolist())
            samples_labels_splits[repeats]['test_set_labels'] = ','.join(label_test.astype(str))
        perf_df = pd.DataFrame.from_dict(perf_dict, orient='index', dtype=None, columns=['Best parameters', 'Validation set indices','Threshold', 'Training performance', 'Validation performance', 'Test performance'])


        return perf_df, best_feature_records, samples_labels_splits

    except OSError as error:
        raise Exception("DeepMicro encountered a problem.")

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

def find_relevant_reactions(dataframe, functional_occurrence_filepath):
    '''
    This function finds all annotations associated with each OTU
    '''
    found_reac = {}
    functional_occurrence = pd.read_csv(functional_occurrence_filepath, sep='\t', index_col=0)

    for organism, row in tqdm(functional_occurrence.iterrows(), desc="Linking annotations to OTUs...", total=len(functional_occurrence)):
        found_reac[organism] = [annot for annot in functional_occurrence.columns if row[annot] > 0]

    dataframe["Linked_annotations"] = dataframe["ID"].map(found_reac)
    return dataframe

def get_info_taxons(otu_db, esmecata_input, functional_occurrence_filepath):
    
    list_of_otus = list(otu_db.index)

    if esmecata_input is not None:
        annots_with_names = add_otu_names(list_of_otus, esmecata_input)
    else:
        annots_with_names = pd.DataFrame(list(zip(list_of_otus, list_of_otus)), columns=["ID","Name"])

    if functional_occurrence_filepath is not None:
        annots_with_names_and_associated_otus = find_relevant_reactions(annots_with_names, functional_occurrence_filepath)
    else:
        annots_with_names_and_associated_otus = annots_with_names

    return (annots_with_names_and_associated_otus)


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

def add_reaction_names(list_of_annots, output_folder):
    '''
    This function queries the OBO and ExPASy data banks to get the names of the IDs
    '''
    data_folder = os.path.join(output_folder, 'data')

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

    enzyme_dat_file = os.path.join(data_folder, 'enzyme.dat')
    #Check if we have enzyme file, if not download it
    if not(os.path.isfile(enzyme_dat_file)):
        url = "https://ftp.expasy.org/databases/enzyme/enzyme.dat"
        r = requests.get(url, allow_redirects=True)
        open(enzyme_dat_file, 'wb').write(r.content)

    go_basic_file = os.path.join(data_folder, 'go-basic.obo')
    #Same with GO file
    if(not os.path.isfile(go_basic_file)):
        go_obo_url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
        r = requests.get(go_obo_url, allow_redirects=True)
        open(go_basic_file, 'wb').write(r.content)

    go = obo_parser.GODag(go_basic_file)

    handle = open(enzyme_dat_file)
    records = Enzyme.parse(handle)

    reaction_names = []
    for id in tqdm(list_of_annots, desc='Checking out annotation names...'):
        reaction_names.append("Not Found")
        if "GO" in id:
            if id in go:
                go_term = go[id]
                reaction_names[-1] = go_term.name
        else:
            de_found = next((item["DE"] for item in records if item["ID"] == id), "Not Found")
            reaction_names[-1] = de_found

    dataframe = pd.DataFrame(list(zip(list_of_annots, reaction_names)), columns=["ID","Name"])
                             
    return dataframe

def find_relevant_otus(dataframe, functional_occurrence_filepath, otu_name_df):
    '''
    This function finds and names all OTUs associated with each annotation
    '''
    found_otu = {}
    found_otu_named = {}
    named_annotations = set(dataframe["ID"].tolist())
    functional_occurrence = pd.read_csv(functional_occurrence_filepath, sep='\t', index_col=0)

    for organism, row in tqdm(functional_occurrence.iterrows(), desc="Linking OTUs to annotations...", total=len(functional_occurrence)):
        organism_name_translated = otu_name_df[otu_name_df['observation_name'] == organism]['taxonomic_affiliation'].values[0]
        organism_name_translated_species = organism_name_translated.split(';')[-1]
        annotations_found = [single_annot for annot in functional_occurrence.columns if row[annot] > 0 for single_annot in [annot]*row[annot] if single_annot in named_annotations]
        for annotation in annotations_found:
            if annotation not in found_otu:
                found_otu[annotation] = [organism]
                found_otu_named[annotation] = [organism_name_translated_species]
            else:
                if organism not in found_otu[annotation]:
                    found_otu[annotation].append(organism)
                if organism_name_translated_species not in found_otu_named[annotation]:
                    found_otu_named[annotation].append(organism_name_translated_species)

    dataframe["Linked_OTUs"] = dataframe["ID"].map(found_otu)
    dataframe["Named_linked_OTUs"] = dataframe["ID"].map(found_otu_named)
    return dataframe


def get_info_annots(score_db, output_folder, otu_abundance_filepath, esmecata_input, functional_occurrence_filepath):
    
    list_of_annots = list(score_db.index)
    annots_with_names = add_reaction_names(list_of_annots, output_folder)

    if otu_abundance_filepath is not None and functional_occurrence_filepath is not None:
        annots_with_names_and_associated_otus = find_relevant_otus(annots_with_names, functional_occurrence_filepath, esmecata_input)
    else:
        annots_with_names_and_associated_otus = annots_with_names
        annots_with_names['Linked_OTUs'] = None
        annots_with_names['Named_linked_OTUs'] = None

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

def averaging_and_info_step(functional_profile_df, label_refs, output_folder, esmecata_input=None, functional_occurrence_filepath=None, otu_abundance_filepath=None):
    '''
    This script runs the steps to build an info database on the taxons and annotations from the dataset. Said databases will be used to add information to the final outputs.
    '''
    if otu_abundance_filepath is not None:
        otu_count_df = pd.read_csv(otu_abundance_filepath, sep='\t', index_col=0)
        info_taxons = get_info_taxons(otu_count_df, esmecata_input, functional_occurrence_filepath)
        avg_count_total = average_count_per_group(otu_count_df, label_refs)

    info_annots = get_info_annots(functional_profile_df, output_folder, otu_abundance_filepath, esmecata_input, functional_occurrence_filepath)
    avg_score_total = average_per_group(functional_profile_df, label_refs)

    df_max_scores = pd.DataFrame()
    df_max_counts = pd.DataFrame()

    for label in np.unique(label_refs.values):
        if otu_abundance_filepath is not None:
            avg_count_group = average_count_per_group(otu_count_df, label_refs, label)
            df_max_counts[label] = avg_count_group['average']
            info_taxons['Average presence in '+str(label)] = avg_count_group['average'].values

        avg_score_group = average_per_group(functional_profile_df, label_refs, label)
        df_max_scores[label] = avg_score_group['average']
        info_annots['Average presence in '+str(label)] = avg_score_group['average'].values
        
    
    if otu_abundance_filepath is not None:
        info_taxons['Representative_group'] = df_max_counts.idxmax(axis=1).values
        info_taxons['Average presence (total)'] = avg_count_total['average'].values

    else:
        info_taxons = None

    temp_records_sets_folder = os.path.join(output_folder, 'temp_records_sets')
    if not os.path.exists(temp_records_sets_folder):
        os.mkdir(temp_records_sets_folder)
    df_max_score_filepath = os.path.join(temp_records_sets_folder, 'df_max_scores.csv')
    df_max_scores.to_csv(df_max_score_filepath)
    info_annots['Representative_group'] = df_max_scores.idxmax(axis=1).values
    info_annots['Average presence (total)'] = avg_score_total['average'].values
    
    return info_annots, info_taxons


def run_iterate(functional_profile_filepath, label_filepath, run_output_folder, run_nb, nb_iterations, esmecata_input=None, esmecata_annotation_reference=None,
                otu_abundance_filepath=None, reference_test_sets_filepath=None,
                classifiers=20, method='rf', var_ranking_method='gini', seed_rf=0, seed_split=0, seed_valid=0):
    functional_profile_df = pd.read_csv(functional_profile_filepath, sep=',', index_col=0)

    label_file_df = pd.read_csv(label_filepath)
    label_file_df = label_file_df[functional_profile_df.columns].transpose()

    test_set_dict = {}
    bank_of_selections_annots = {}
    bank_of_selections_taxons = {}
    bank_of_performance_dfs_annots = {}
    bank_of_performance_dfs_taxons = {}
    bank_of_average_importances_annots = {}
    bank_of_average_importances_taxons = {}

    ## Calculating average presence of taxons and annotations per label, and collecting info about them
    info_annots, info_taxons = averaging_and_info_step(functional_profile_df, label_file_df, run_output_folder, esmecata_input, esmecata_annotation_reference, otu_abundance_filepath)
    info_annots.to_csv(run_output_folder+'/info_annots_check.csv')

    if reference_test_sets_filepath:
        #Get the test set references if they are given
        test_set_refs = pd.read_csv(reference_test_sets_filepath)
        test_labels = test_set_refs['Run_'+str(run_nb)].values

        sample_names = list(functional_profile_df.columns.values)

        labels_test_noshape = label_file_df.loc[test_labels]
        labels_test = labels_test_noshape.values.reshape((label_file_df.loc[test_labels].values.shape[0]))
        #logger.info('Label_test:', labels_test)
        train_samples = []
        for i in sample_names:
            if i not in test_labels:
                train_samples.append(i)

        labels_train_noshape = label_file_df.loc[train_samples]
        labels_train = labels_train_noshape.values.reshape((label_file_df.loc[train_samples].values.shape[0]))

    else:
        #Select a test subset
        labels_train, labels_test, sample_train_indices, test_indices = separate_test_train(label_file_df, seed_split)
        train_samples = label_file_df.iloc[sample_train_indices].index.tolist()

        test_labels = [functional_profile_df.columns[i] for i in test_indices]

    #Keeping track of the selected test sets (and writing them at every run, in case of a crash)
    test_set_dict = {}
    test_set_dict['Run_'+str(run_nb)] = test_labels
    test_set_df = pd.DataFrame.from_dict(test_set_dict)
    test_set_output_file = os.path.join(run_output_folder, 'Test_sets.csv')
    test_set_df.to_csv(test_set_output_file)

    if otu_abundance_filepath is not None:
        deepmicro_otu_iteration0 = pd.read_csv(otu_abundance_filepath, sep='\t', index_col=0).transpose()
    deepmicro_sofa_iteration0 = functional_profile_df.transpose()

    trained_classifiers_folder = os.path.join(run_output_folder, 'Trained_classifiers')
    if not os.path.exists(trained_classifiers_folder):
        os.mkdir(trained_classifiers_folder)
    classification_performances_folder = os.path.join(run_output_folder, 'Classification_performances')
    if not os.path.exists(classification_performances_folder):
        os.mkdir(classification_performances_folder)
    selected_variables_folder = os.path.join(run_output_folder, 'Selected_Variables')
    if not os.path.exists(selected_variables_folder):
        os.mkdir(selected_variables_folder)
    dataset_separation_folder = os.path.join(run_output_folder, 'Dataset_separation')
    if not os.path.exists(dataset_separation_folder):
        os.mkdir(dataset_separation_folder)

    random.seed(seed_rf)
    
    seed_rf_vec = random.sample(range(1000),nb_iterations)
    
    if otu_abundance_filepath is not None:
        deepmicro_otu_iteration = pd.DataFrame()
        deepmicro_otu_iteration = deepmicro_otu_iteration0
    deepmicro_sofa_iteration = deepmicro_sofa_iteration0
    
    #ITERATED:
    for iteration_number in range(nb_iterations):
    #Separate test and train subsets
        if otu_abundance_filepath is not None:
            otu_train, otu_test = create_test_train_subsets(deepmicro_otu_iteration, test_labels)
            otu_train = deepmicro_otu_iteration.loc[train_samples]
        annots_train, annots_test = create_test_train_subsets(deepmicro_sofa_iteration, test_labels)
        annots_train = deepmicro_sofa_iteration.loc[train_samples]

        #DeepMicro
        trained_classifiers_iteration_folder = os.path.join(trained_classifiers_folder, 'Iteration_'+str(iteration_number))
        if not os.path.exists(trained_classifiers_iteration_folder):
            os.mkdir(trained_classifiers_iteration_folder)

        # Run classification wtih DeepMicro.
        if otu_abundance_filepath is not None:
            perf_df_otu, best_feature_records_otu, taxon_training_validation_sets = run_deep_micro(otu_test, otu_train, labels_test, labels_train, 'test_OTU', iteration_number, run_nb, trained_classifiers_iteration_folder, "Taxonomic",
                                                                    classifiers, method, var_ranking_method, real_seed=seed_rf_vec[nb_iterations-1], seed_valid=seed_valid)
        perf_df_sofa, best_feature_records_sofa, function_training_validation_sets = run_deep_micro(annots_test, annots_train, labels_test, labels_train, 'test_Functions', iteration_number, run_nb, trained_classifiers_iteration_folder, "Functional",
                                                                  classifiers, method, var_ranking_method, real_seed=seed_rf_vec[nb_iterations-1], seed_valid=seed_valid)
        if otu_abundance_filepath is not None:
            taxon_dataset_separation_iteration_file = os.path.join(dataset_separation_folder, 'Taxonomic_samples_separation_Iteration_'+str(iteration_number)+'.csv')
            pd.DataFrame(taxon_training_validation_sets).to_csv(taxon_dataset_separation_iteration_file)

        function_dataset_separation_iteration_file = os.path.join(dataset_separation_folder, 'Annotation_samples_separation_Iteration_'+str(iteration_number)+'.csv')
        pd.DataFrame(function_training_validation_sets).to_csv(function_dataset_separation_iteration_file)

        classification_performances_iteration_folder = os.path.join(classification_performances_folder, 'Iteration_'+str(iteration_number))
        if not os.path.exists(classification_performances_iteration_folder):
            os.mkdir(classification_performances_iteration_folder)
        
        if otu_abundance_filepath is not None:
            taxonomic_performances_file = os.path.join(classification_performances_iteration_folder, 'Taxonomic_performances.csv')
            perf_df_otu.to_csv(taxonomic_performances_file)
        annotation_performances_file = os.path.join(classification_performances_iteration_folder, 'Annotation_performances.csv')
        perf_df_sofa.to_csv(annotation_performances_file)

        if otu_abundance_filepath is not None:
            if iteration_number not in bank_of_performance_dfs_taxons:
                bank_of_performance_dfs_taxons[iteration_number] = {}
            bank_of_performance_dfs_taxons[iteration_number][run_nb] = perf_df_otu
        if iteration_number not in bank_of_performance_dfs_annots:
            bank_of_performance_dfs_annots[iteration_number] = {}
        bank_of_performance_dfs_annots[iteration_number] [run_nb] = perf_df_sofa

        #If a variable ranking was made (not SVM): select, record and prepare datasets for next iteration 
        if method == "rf":
            if otu_abundance_filepath is not None:
                best_feature_records_otu_df = pd.concat(best_feature_records_otu, axis=1)
                best_feature_records_otu_df.columns = ['Model_'+str(i) for i in range(int(classifiers))]
                best_feature_records_otu_df.index = deepmicro_otu_iteration.columns
                #best_feature_records_otu_df.to_csv(pipeline_path+'/Meta-Outputs/'+data_ref_output_name+'/Run_'+str(run_nb)+'/Classification_performances/Iteration_'+str(iteration_number)+'/Best_feature_records_taxons.csv')
                
            best_feature_records_sofa_df = pd.concat(best_feature_records_sofa, axis=1)
            best_feature_records_sofa_df.columns = ['Model_'+str(i) for i in range(int(classifiers))]
            best_feature_records_sofa_df.index = deepmicro_sofa_iteration.columns
            #best_feature_records_sofa_df.to_csv(pipeline_path+'/Meta-Outputs/'+data_ref_output_name+'/Run_'+str(run_nb)+'/Classification_performances/Iteration_'+str(iteration_number)+'/Best_feature_records_annotations.csv')

            #Average Gini importances
            if otu_abundance_filepath is not None:
                best_feature_records_otu_df['Average'] = best_feature_records_otu_df.mean(axis=1)
                feature_importance_records_taxons_filepath = os.path.join(classification_performances_iteration_folder, 'Feature_importance_records_taxons.csv')
                best_feature_records_otu_df.to_csv(feature_importance_records_taxons_filepath)
                if iteration_number not in bank_of_average_importances_taxons:
                    bank_of_average_importances_taxons[iteration_number] = {}
                bank_of_average_importances_taxons[iteration_number][run_nb] = best_feature_records_otu_df['Average']
            best_feature_records_sofa_df['Average'] = best_feature_records_sofa_df.mean(axis=1)
            feature_importance_records_annotations_filepath = os.path.join(classification_performances_iteration_folder, 'Feature_importance_records_annotations.csv')
            best_feature_records_sofa_df.to_csv(feature_importance_records_annotations_filepath)
            if iteration_number not in bank_of_average_importances_annots:
                bank_of_average_importances_annots[iteration_number] = {}
            bank_of_average_importances_annots[iteration_number][run_nb] = best_feature_records_sofa_df['Average']

            #Rotor cutoff
            selected_variables_iteration_folder = os.path.join(selected_variables_folder, 'Iteration_'+str(iteration_number))
            if not os.path.exists(selected_variables_iteration_folder):
                os.mkdir(selected_variables_iteration_folder)

            if otu_abundance_filepath is not None:
                retained_otus = inflexion_cutoff(best_feature_records_otu_df)
            retained_annots = inflexion_cutoff(best_feature_records_sofa_df)

            if otu_abundance_filepath is not None:
                selection_plus_info_taxons = info_taxons[info_taxons['ID'].isin(list(retained_otus.index))]
                selection_plus_info_taxons['Average_importance'] = [best_feature_records_otu_df.loc[tax, 'Average'] for tax in retained_otus.index]
            selection_plus_info_annots = info_annots[info_annots['ID'].isin(list(retained_annots.index))]
            selection_plus_info_annots['Average_importance'] = [best_feature_records_sofa_df.loc[func, 'Average'] for func in retained_annots.index]

            # Write the selection files with info
            if otu_abundance_filepath is not None:
                signif_otus = []
                signif_otu_names = []

                for link_otu_list in selection_plus_info_annots['Linked_OTUs'].values:
                    signif_links = []
                    signif_links_named = []
                    if link_otu_list is not None:
                        
                        for otu in tqdm(link_otu_list, desc="link_otu_list: "+str(link_otu_list)+", type: "+str(type(link_otu_list))):
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

                if 'Linked_annotations' in selection_plus_info_taxons:
                    for link_annot_list in selection_plus_info_taxons['Linked_annotations'].values:
                        signif_links = []

                        for annot in link_annot_list:
                            if annot in retained_annots:
                                signif_links.append(annot)
                        signif_annots.append(signif_links)
                else:
                    for otu in retained_otus:
                        signif_annots.append([])

                selected_taxons_run_file = os.path.join(selected_variables_iteration_folder, 'Selected_taxons_run_'+str(run_nb)+'_iter_'+str(iteration_number)+'.csv')
                if 'Linked_annotations' in selection_plus_info_taxons:
                    selection_plus_info_taxons['Significant_linked_annotations'] = signif_annots
                selection_plus_info_taxons.to_csv(selected_taxons_run_file)

            selected_annotations_run_file = os.path.join(selected_variables_iteration_folder, 'Selected_annotations_run_'+str(run_nb)+'_iter_'+str(iteration_number)+'.csv')
            selection_plus_info_annots.to_csv(selected_annotations_run_file)

            # Store the selections for later incorporation to the Core/Meta sets
            if otu_abundance_filepath is not None:
                if iteration_number not in bank_of_selections_taxons:
                    bank_of_selections_taxons[iteration_number] = {}
                bank_of_selections_taxons[iteration_number][run_nb] = list(retained_otus.index)
            if iteration_number not in bank_of_selections_annots:
                bank_of_selections_annots[iteration_number] = {}
            bank_of_selections_annots[iteration_number][run_nb] = list(retained_annots.index)

            # Update deepmicro_sofa_iteration
            if otu_abundance_filepath is not None:
                deepmicro_otu_iteration = deepmicro_otu_iteration0[retained_otus.index]
            deepmicro_sofa_iteration = deepmicro_sofa_iteration0[retained_annots.index]
            
    return (test_set_dict, bank_of_selections_annots, bank_of_selections_taxons, bank_of_performance_dfs_annots, bank_of_performance_dfs_taxons, bank_of_average_importances_annots, bank_of_average_importances_taxons)
