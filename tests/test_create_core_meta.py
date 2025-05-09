import shutil
import pandas as pd
import os

import subprocess

from sparta_pipeline.visualize import plot_classifs
from sparta_pipeline.create_core_meta import extract_and_write_core_meta

def get_run_bank_data(sparta_output_folder):
    bank_of_selections_annots = {}
    bank_of_selections_taxons = {}
    bank_of_performance_dfs_annots = {}
    bank_of_performance_dfs_taxons = {}
    bank_of_average_importances_annots = {}
    bank_of_average_importances_taxons = {}

    run_folders = [file_folder for file_folder in os.listdir(sparta_output_folder) if file_folder.startswith('Run_')]
    for run_folder in run_folders:
        run_number = run_folder.split('_')[1]
        # Get selected variables.
        selected_variables_folder_path = os.path.join(sparta_output_folder, run_folder, 'Selected_Variables')
        for iteration_folder in os.listdir(selected_variables_folder_path):
            iteration_number = int(iteration_folder.split('_')[1])
            if iteration_number not in bank_of_selections_annots:
                bank_of_selections_annots[iteration_number] = {}
            if iteration_number not in bank_of_selections_taxons:
                bank_of_selections_taxons[iteration_number] = {}
            iteration_folder_path = os.path.join(selected_variables_folder_path, iteration_folder)
            annotation_file = [os.path.join(iteration_folder_path, iteration_file) for iteration_file in os.listdir(iteration_folder_path) if 'annotation' in iteration_file][0]
            annotation_df = pd.read_csv(annotation_file, index_col=0)
            bank_of_selections_annots[iteration_number][run_number] = annotation_df['ID'].tolist()
            taxons_file = [os.path.join(iteration_folder_path, iteration_file) for iteration_file in os.listdir(iteration_folder_path) if 'taxons' in iteration_file]
            if len(taxons_file) > 0:
                taxons_file = taxons_file[0]
                taxons_df = pd.read_csv(taxons_file, index_col=0)
                bank_of_selections_taxons[iteration_number][run_number] = taxons_df['ID'].tolist()

        # Get the performance.
        classification_performance_folder_path = os.path.join(sparta_output_folder, run_folder, 'Classification_performances')
        for iteration_folder in os.listdir(classification_performance_folder_path):
            iteration_number = int(iteration_folder.split('_')[1])
            if iteration_number not in bank_of_performance_dfs_annots:
                bank_of_performance_dfs_annots[iteration_number] = {}
            if iteration_number not in bank_of_performance_dfs_taxons:
                bank_of_performance_dfs_taxons[iteration_number] = {}
            if iteration_number not in bank_of_average_importances_annots:
                bank_of_average_importances_annots[iteration_number] = {}
            if iteration_number not in bank_of_average_importances_taxons:
                bank_of_average_importances_taxons[iteration_number] = {}
            iteration_folder_path = os.path.join(classification_performance_folder_path, iteration_folder)
            annotation_performance = os.path.join(iteration_folder_path, 'Annotation_performances.csv')
            taxonomic_performance = os.path.join(iteration_folder_path, 'Taxonomic_performances.csv')
            bank_of_performance_dfs_annots[iteration_number][run_number] = pd.read_csv(annotation_performance, index_col=0)
            if os.path.exists(taxonomic_performance):
                bank_of_performance_dfs_taxons[iteration_number][run_number] = pd.read_csv(taxonomic_performance, index_col=0)
            annotation_importance = os.path.join(iteration_folder_path, 'Feature_importance_records_annotations.csv')
            bank_of_average_importances_annots[iteration_number][run_number] = pd.read_csv(annotation_importance, index_col=0)['Average']
            taxonomic_importance = os.path.join(iteration_folder_path, 'Feature_importance_records_taxons.csv')
            if os.path.exists(taxonomic_importance):
                bank_of_average_importances_taxons[iteration_number][run_number] = pd.read_csv(taxonomic_importance, index_col=0)['Average']

    return bank_of_selections_annots, bank_of_selections_taxons, bank_of_performance_dfs_annots, bank_of_performance_dfs_taxons, bank_of_average_importances_annots, bank_of_average_importances_taxons

def test_extract_and_write_core_meta():
    sparta_output_folder = os.path.join('expected', 'out')
    organism_abundance_filepath = os.path.join('input', 'test_taxon_profile.tsv')
    functional_profile_df = os.path.join('input', 'test_functional_profile.csv')
    esmecata_file_path = os.path.join('input', 'test_esmecata_input.tsv')
    esmecata_input = pd.read_csv(esmecata_file_path, sep='\t')
    nb_runs = 2
    info_annots = pd.read_csv(os.path.join(sparta_output_folder, 'info_annots_check.csv'))
    info_taxons = pd.read_csv(os.path.join(sparta_output_folder, 'info_taxons_check.csv'))

    core_and_meta_outputs_folder = 'output_folder'
    if not os.path.exists(core_and_meta_outputs_folder):
        os.mkdir(core_and_meta_outputs_folder)
    all_iterations_output = os.path.join(core_and_meta_outputs_folder, 'All_iterations')
    if not os.path.exists(all_iterations_output):
        os.mkdir(all_iterations_output)
    best_iteration_output = os.path.join(core_and_meta_outputs_folder, 'Best_iteration')
    if not os.path.exists(best_iteration_output):
        os.mkdir(best_iteration_output)

    bank_of_selections_annots, bank_of_selections_taxons, bank_of_performance_dfs_annots, bank_of_performance_dfs_taxons, bank_of_average_importances_annots, bank_of_average_importances_taxons = get_run_bank_data(sparta_output_folder)
    print(bank_of_selections_annots)
    print(bank_of_selections_taxons)
    print(bank_of_performance_dfs_annots)
    print(bank_of_performance_dfs_taxons)
    print(bank_of_average_importances_annots)
    print(bank_of_average_importances_taxons)
    visualisation_file = os.path.join(core_and_meta_outputs_folder, 'median_OTU_vs_SoFA_(best_vs_best).png')
    #visualisation_file_v2 = os.path.join(output_folder, 'median_OTU_vs_SoFA_(best_vs_best)_v2.png')
    best_selec_iter_annots, best_selec_iter_taxons = plot_classifs(bank_of_performance_dfs_annots, bank_of_performance_dfs_taxons, 'Test set', visualisation_file,
                                                                organism_abundance_filepath)

    df_perfs_and_selection_per_iter, warning_annots, warning_taxons = extract_and_write_core_meta(core_and_meta_outputs_folder, bank_of_selections_annots, bank_of_selections_taxons, bank_of_performance_dfs_annots,
                                                                                                bank_of_performance_dfs_taxons, bank_of_average_importances_annots, bank_of_average_importances_taxons,
                                                                                                best_selec_iter_annots, best_selec_iter_taxons,
                                                                                                info_annots, info_taxons, nb_runs, esmecata_input, functional_profile_df, organism_abundance_filepath)
    

if __name__ == "__main__":
    test_extract_and_write_core_meta()
