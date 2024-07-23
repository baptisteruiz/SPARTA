import shutil
import pandas as pd
import os

import subprocess

def test_cli_sparta_classification():
    functional_profile_filepath = os.path.join('input', 'test_functional_profile.csv')
    label_filepath = os.path.join('input', 'test_label.csv')
    output_folder = 'output_folder'
    run_nb = '1'
    nb_iterations = '1'
    nb_classifier = '1'
    subprocess.call(['sparta', 'classification', '-fp', functional_profile_filepath, '-l', label_filepath, '-o', output_folder, '-r', run_nb, '-i', nb_iterations, '-c', nb_classifier])

    # Get label for sample.
    label_file_df = pd.read_csv(label_filepath)
    functional_profile_df = pd.read_csv(functional_profile_filepath, index_col=0)
    label_file_df = label_file_df[functional_profile_df.columns].transpose()
    sample_to_label = label_file_df[0].to_dict()

    # Get test, training and validation labels.
    dataset_separation_filepath = os.path.join(output_folder, 'Run_1_MasterSeed_864', 'Dataset_separation', 'Annotation_samples_separation_Iteration_0.csv')
    dataset_separation = pd.read_csv(dataset_separation_filepath, index_col=0).to_dict()

    training_set_label = [str(sample_to_label[sample]) for sample in dataset_separation['1']['training_set'].split(',')]
    validation_set_label = [str(sample_to_label[sample]) for sample in dataset_separation['1']['validation_set'].split(',')]
    test_set_label = [str(sample_to_label[sample])  for sample in dataset_separation['1']['test_set'].split(',')]

    used_training_set_label = [label for label in dataset_separation['1']['training_set_labels'].split(',')]
    used_validation_set_label = [label for label in dataset_separation['1']['validation_set_labels'].split(',')]
    used_test_set_label = [label for label in dataset_separation['1']['test_set_labels'].split(',')]

    # There was an issue where used sample do not correspond to label ID (due to wrong index selection after splitting train/ validation sets).
    # This lead to mismatch between label for the X dataset and the y label.
    # These assertions ensure that both training and validation are performed with matching X dataset and y label.
    assert used_training_set_label == training_set_label
    assert used_validation_set_label == validation_set_label
    assert used_test_set_label == test_set_label

    shutil.rmtree(output_folder)


def test_cli_pipeline():
    pipeline_input_filepath = os.path.join('input', 'test_input_pipeline.txt')
    label_filepath = os.path.join('input', 'test_label.csv')
    input_folder = os.path.join('input', 'test_pipeline_esmecata')
    run_nb = '1'
    nb_iterations = '1'
    nb_classifier = '1'

    output_folder = 'output_folder'
    shutil.copytree(input_folder, output_folder)
    subprocess.call(['sparta', 'pipeline', '-p', pipeline_input_filepath, '-l', label_filepath, '-o', output_folder, '-r', run_nb, '-i', nb_iterations, '-c', nb_classifier])

    # Check esmecata output.
    # Check function table.
    output_sofa_table_filepath = os.path.join(output_folder, 'functional_occurrence.tsv')
    expected_sofa_table = os.path.join('input', 'test_functional_occurrence.tsv')
    df_expected = pd.read_csv(expected_sofa_table, index_col=0, sep='\t')
    computed_df = pd.read_csv(output_sofa_table_filepath, index_col=0, sep='\t')
    # Check that both files have same index (annotation IDs).
    assert set(df_expected.index) == set(computed_df.index)
    # Make the comparison on the same index.
    computed_df = computed_df.loc[df_expected.index]
    df_expected = df_expected[computed_df.columns]
    assert all(df_expected.compare(computed_df))

    # Check functional profile.
    output_sofa_table_filepath = os.path.join(output_folder, 'SoFA_table.csv')
    expected_sofa_table = os.path.join('expected', 'test_expected_SoFA_table_abundance_test.csv')
    df_expected = pd.read_csv(expected_sofa_table, index_col=0)
    computed_df = pd.read_csv(output_sofa_table_filepath, index_col=0)

    # Check that both files have same index (annotation IDs).
    assert set(df_expected.index) == set(computed_df.index)
    # Make the comparison on the same index.
    computed_df = computed_df.loc[df_expected.index]
    df_expected = df_expected[computed_df.columns]
    assert all(df_expected.compare(computed_df))

    # Check core results.
    expected_core_annots_table = os.path.join('expected', 'Core_annots_iteration_-1.csv')
    expected_core_annots_df = pd.read_csv(expected_core_annots_table)
    expected_core_taxons_table = os.path.join('expected', 'Core_taxons_iteration_-1.csv')
    expected_core_taxons_df = pd.read_csv(expected_core_taxons_table)

    computed_core_annots_table = os.path.join(output_folder, 'Core_and_Meta_outputs', 'Best_iteration', 'Core_annots_iteration_-1.csv')
    computed_core_annots_df = pd.read_csv(computed_core_annots_table)
    computed_core_taxons_table = os.path.join(output_folder, 'Core_and_Meta_outputs', 'Best_iteration', 'Core_taxons_iteration_-1.csv')
    computed_core_taxons_df = pd.read_csv(computed_core_taxons_table)

    assert expected_core_annots_df['ID'].to_list() == computed_core_annots_df['ID'].to_list()
    assert expected_core_taxons_df['ID'].to_list() == computed_core_taxons_df['ID'].to_list()

    shutil.rmtree(output_folder)