import shutil
import pandas as pd
import os

import subprocess

def test_cli_sparta_classification():
    functional_profile_filepath = 'test_functional_profile.csv'
    label_filepath = 'test_label.csv'
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
    dataset_separation_filepath = os.path.join(output_folder, 'Run_1', 'Dataset_separation', 'Annotation_samples_separation_Iteration_0.csv')
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
