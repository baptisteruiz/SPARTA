import os
from datetime import datetime
import pandas as pd
import shutil
import logging
from SPARTA.iteration import run_iterate, averaging_and_info_step
from SPARTA.visualize import plot_classifs
from SPARTA.create_core_meta import extract_and_write_core_meta

logger = logging.getLogger(__name__)

def run_sparta_classification(functional_profile_filepath, label_file, output_folder, nb_runs, nb_iterations,
                            esmecata_input=None, esmecata_annotation_reference=None, otu_abundance_filepath=None, reference_test_sets_filepath=None,
                            classifiers=3, method='rf', var_ranking_method='gini', keep_temp=None):
    date_time_now = datetime.now()
    ref_time = date_time_now
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    stopwatch_file = os.path.join(output_folder, 'stopwatch.txt')
    functional_profile_df = pd.read_csv(functional_profile_filepath, sep=',', index_col=0)

    label_file_df = pd.read_csv(label_file)
    label_file_df = label_file_df[functional_profile_df.columns].transpose()
    ## Calculating average presence of taxons and annotations per label, and collecting info about them.
    if esmecata_input is not None:
        esmecata_input = pd.read_csv(esmecata_input, sep='\t')
    info_annots, info_taxons = averaging_and_info_step(functional_profile_df, label_file_df, output_folder, esmecata_input, esmecata_annotation_reference, otu_abundance_filepath)

    nb_runs = int(nb_runs)
    nb_iterations = int(nb_iterations)
    bank_of_selections_annots = {}
    bank_of_selections_taxons = {}
    bank_of_performance_dfs_annots = {}
    bank_of_performance_dfs_taxons = {}
    test_set_dict = {}
    ## Launching the SPARTA runs
    for run_nb in range(1, nb_runs+1):
        run_output_folder = os.path.join(output_folder, 'Run_'+str(run_nb))
        if not os.path.exists(run_output_folder):
            os.mkdir(run_output_folder)
        run_test_set_dict, run_bank_of_selections_annots, run_bank_of_selections_taxons, run_bank_of_performance_dfs_annots, run_bank_of_performance_dfs_taxons = run_iterate(functional_profile_filepath, label_file, run_output_folder,
                                                                                                                                                          run_nb, nb_iterations, esmecata_input, esmecata_annotation_reference,
                                                                                                                                                          otu_abundance_filepath, reference_test_sets_filepath,
                                                                                                                                                          classifiers, method, var_ranking_method)
        bank_of_selections_annots.update(run_bank_of_selections_annots)
        bank_of_selections_taxons.update(run_bank_of_selections_taxons)
        bank_of_performance_dfs_annots.update(run_bank_of_performance_dfs_annots)
        bank_of_performance_dfs_taxons.update(run_bank_of_performance_dfs_taxons)
        test_set_dict.update(run_test_set_dict)

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
            logger.info('WARNING: functional classification performance is low (<0.6). The selected functional variables may not be accurate representatives of the classification task')
        
        if warning_taxons:
            logger.info('WARNING: taxonomic classification performance is low (<0.6). The selected taxonomic variables may not be accurate representatives of the classification task')

    if not keep_temp:
        shutil.rmtree('/Outputs_temp/', ignore_errors=True)
    # pd.DataFrame.from_dict(bank_of_selections_annots).to_csv(pipeline_path+'/Meta-Outputs/'+data_ref_output_name+'/bank_of_selections_check.csv')
