import shutil
import pandas as pd
import os

from SPARTA.iteration import run_iterate

def test_run_iterate():
    functional_profile_filepath = 'SoFA_table.csv'
    label_filepath = 'Label_abundance_test.csv'
    run_output_folder = 'output_folder'
    run_nb = 1
    nb_iterations = 1
    if not os.path.exists(run_output_folder):
        os.mkdir(run_output_folder)

    run_iterate(functional_profile_filepath, label_filepath, run_output_folder, run_nb, nb_iterations, classifiers=1)

    # Get label for sample.
    label_file_df = pd.read_csv(label_filepath)
    functional_profile_df = pd.read_csv(functional_profile_filepath, index_col=0)
    label_file_df = label_file_df[functional_profile_df.columns].transpose()
    sample_to_label = label_file_df[0].to_dict()

    # Get test set labels.
    test_set_filepath = os.path.join('output_folder', 'Test_sets.csv')
    test_set = pd.read_csv(test_set_filepath)
    test_set_label = [str(sample_to_label[sample]) for sample in test_set['Run_1'].tolist()]

    # Get training and validation labels.
    dataset_separation_filepath = os.path.join('output_folder', 'Dataset_separation', 'Iteration_0.csv')
    dataset_separation = pd.read_csv(dataset_separation_filepath, index_col=0).to_dict()

    training_set_label = [str(sample_to_label[sample]) for sample in dataset_separation['1']['training_set'].split(',')]
    validation_set_label = [str(sample_to_label[sample]) for sample in dataset_separation['1']['validation_set'].split(',')]
    used_training_set_label = [label for label in dataset_separation['1']['training_set_labels'].split(',')]
    used_validation_set_label = [label for label in dataset_separation['1']['validation_set_labels'].split(',')]

    # There was an issue where used sample do not correspond to label ID (due to wrong index selection after splitting train/ validation sets).
    # This lead to mismatch between label for the X dataset and the y label.
    # These assertions ensure that both training and validation are performed with matching X dataset and y label.
    assert used_training_set_label == training_set_label
    assert used_validation_set_label == validation_set_label


    #shutil.rmtree(run_output_folder)