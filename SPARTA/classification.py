import os
from datetime import datetime
import pandas as pd
import shutil
import logging
from SPARTA.iteration import run_iterate, averaging_and_info_step
from SPARTA.visualize import plot_classifs
from SPARTA.create_core_meta import extract_and_write_core_meta

logger = logging.getLogger(__name__)


def update_iteration_dict(general_dict, iteration_dict):
    """ Merge iteration dict from a run into general dict for all runs.
    Args:
        general_dict (dict): Dictionary containing all information from the different runs.
        iteration_dict (dict): Dictionary containing all results from all iterations of a run.
    Returns:
        general_dict (dict): Updated dictionary containing all information from the different runs.
    """
    for iteration_nb in iteration_dict:
        if iteration_nb not in general_dict:
            general_dict[iteration_nb] = iteration_dict[iteration_nb]
        else:
            general_dict[iteration_nb].update(iteration_dict[iteration_nb])
    return general_dict


def run_sparta_classification(functional_profile_filepath, label_filepath, output_folder, nb_runs, nb_iterations,
                            esmecata_input=None, functional_occurrence_filepath=None, otu_abundance_filepath=None, reference_test_sets_filepath=None,
                            classifiers=20, method='rf', var_ranking_method='gini', keep_temp=None):
    """ Run the classification part of SPARTA using either:
        - (1) a functional profile file (associating functions to samples) with a label file associating samples and labels.
        - (2) a functional profile file, a label file, esmecata input file, functional occurrence file and abundanc of organisms.

    Args:
        functional_profile_filepath (str): Path to the functional profile file.
        label_filepath (str): Path to the label file.
        output_folder (str): Path to the output folder.
        nb_runs (int): Number of runs to perform (a run consist of X iterations).
        nb_iterations (int): Number of iterations per run.
        esmecata_input (str): Path to esmecata input file.
        functional_occurrence_filepath (str): Path to functional occurrence file.
        otu_abundance_filepath (str): Path to the abundance file associating Organism and their abundances in samples.
        reference_test_sets_filepath (str): Path to a file indicating the test sets to used.
        classifiers (int): Number of classifiers to use in a Random Forests.
        method (str): Classifier to use (either rf or svm).
        var_ranking_method (str): Variable ranking method (gini or shap).
        keep_temp (bool): Bool to keep or not temporary files.
    """
    logger.info('SPARTA|classification| Begin classification.')
    date_time_now = datetime.now()
    ref_time = date_time_now
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    stopwatch_file = os.path.join(output_folder, 'stopwatch.txt')
    functional_profile_df = pd.read_csv(functional_profile_filepath, sep=',', index_col=0)

    label_file_df = pd.read_csv(label_filepath)
    label_file_df = label_file_df[functional_profile_df.columns].transpose()
    ## Calculating average presence of taxons and annotations per label, and collecting info about them.
    if esmecata_input is not None:
        esmecata_input = pd.read_csv(esmecata_input, sep='\t')
    info_annots, info_taxons = averaging_and_info_step(functional_profile_df, label_file_df, output_folder, esmecata_input, functional_occurrence_filepath, otu_abundance_filepath)

    nb_runs = int(nb_runs)
    nb_iterations = int(nb_iterations)
    bank_of_selections_annots = {}
    bank_of_selections_taxons = {}
    bank_of_performance_dfs_annots = {}
    bank_of_performance_dfs_taxons = {}
    test_set_dict = {}
    ## Launching the SPARTA runs
    for run_nb in range(1, nb_runs+1):
        logger.info('SPARTA|classification| Run number {0}.'.format(run_nb))
        run_output_folder = os.path.join(output_folder, 'Run_'+str(run_nb))
        if not os.path.exists(run_output_folder):
            os.mkdir(run_output_folder)
        # Launch the different iterations for the run.
        run_test_set_dict, run_bank_of_selections_annots, run_bank_of_selections_taxons, run_bank_of_performance_dfs_annots, run_bank_of_performance_dfs_taxons = run_iterate(functional_profile_filepath, label_filepath, run_output_folder,
                                                                                                                                                          run_nb, nb_iterations, esmecata_input, functional_occurrence_filepath,
                                                                                                                                                          otu_abundance_filepath, reference_test_sets_filepath,
                                                                                                                                                          classifiers, method, var_ranking_method)
        bank_of_selections_annots = update_iteration_dict(bank_of_selections_annots, run_bank_of_selections_annots)
        bank_of_selections_taxons = update_iteration_dict(bank_of_selections_taxons, run_bank_of_selections_taxons)
        bank_of_performance_dfs_annots = update_iteration_dict(bank_of_performance_dfs_annots, run_bank_of_performance_dfs_annots)
        bank_of_performance_dfs_taxons = update_iteration_dict(bank_of_performance_dfs_taxons, run_bank_of_performance_dfs_taxons)
        test_set_dict = update_iteration_dict(test_set_dict, run_test_set_dict)

        ####Time measurement####
        date_time_now = datetime.now()
        run_i_time = date_time_now - ref_time
        run_i_time_seconds = run_i_time.total_seconds()

        f = open(stopwatch_file, "a")
        f.write("SPARTA Run "+str(run_nb)+" length (s): "+str(run_i_time_seconds)+"\n")
        f.close()

        ref_time = date_time_now
        ########################

    visualisation_file = os.path.join(output_folder, 'median_OTU_vs_SoFA_(best_vs_best).png')
    best_selec_iter_annots, best_selec_iter_taxons = plot_classifs(bank_of_performance_dfs_annots, bank_of_performance_dfs_taxons, 'test', visualisation_file, otu_abundance_filepath)

    ##We only calculate Core and Meta selections if there was a variable selection (i.e: not SVM)
    if method == 'rf':
        logger.info('SPARTA|classification| Compute Core and Meta selections.')
        core_and_meta_outputs_folder = os.path.join(output_folder, 'Core_and_Meta_outputs')
        if not os.path.exists(core_and_meta_outputs_folder):
            os.mkdir(core_and_meta_outputs_folder)
        core_and_meta_outputs_all_iteration_folder = os.path.join(core_and_meta_outputs_folder, 'All_iterations')
        if not os.path.exists(core_and_meta_outputs_all_iteration_folder):
            os.mkdir(core_and_meta_outputs_all_iteration_folder)
        core_and_meta_outputs_best_iteration_folder = os.path.join(core_and_meta_outputs_folder, 'Best_iteration')
        if not os.path.exists(core_and_meta_outputs_best_iteration_folder):
            os.mkdir(core_and_meta_outputs_best_iteration_folder)

        df_perfs_and_selection_per_iter, warning_annots, warning_taxons = extract_and_write_core_meta(core_and_meta_outputs_folder, bank_of_selections_annots, bank_of_selections_taxons, bank_of_performance_dfs_annots,
                                                                                                      bank_of_performance_dfs_taxons, best_selec_iter_annots, best_selec_iter_taxons,
                                                                                                      info_annots, info_taxons, nb_runs, esmecata_input, functional_profile_df, label_file_df, otu_abundance_filepath)
        overall_selection_and_performance_metrics_filepath = os.path.join(output_folder, 'Overall_selection_and_performance_metrics.csv')
        pd.DataFrame.from_dict(df_perfs_and_selection_per_iter).to_csv(overall_selection_and_performance_metrics_filepath)

        if warning_annots:
            logger.info('SPARTA|classification| WARNING: functional classification performance is low (<0.6). The selected functional variables may not be accurate representatives of the classification task')
        
        if warning_taxons:
            logger.info('SPARTA|classification| WARNING: taxonomic classification performance is low (<0.6). The selected taxonomic variables may not be accurate representatives of the classification task')

    if not keep_temp:
        shutil.rmtree('/Outputs_temp/', ignore_errors=True)
    # pd.DataFrame.from_dict(bank_of_selections_annots).to_csv(pipeline_path+'/Meta-Outputs/'+data_ref_output_name+'/bank_of_selections_check.csv')
