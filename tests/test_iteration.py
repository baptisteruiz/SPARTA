import shutil
import pandas as pd
import os

from SPARTA.iteration import run_iterate

def test_run_iterate():
    functional_profile_filepath = os.path.join('input', 'test_functional_profile.csv')
    label_filepath = os.path.join('input', 'test_label.csv')
    output_folder = 'output_folder'
    run_nb = 1
    nb_iterations = 1
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    run_iterate(functional_profile_filepath, label_filepath, output_folder, run_nb, nb_iterations, classifiers=1)
    
    # Get label for sample.
    label_file_df = pd.read_csv(label_filepath)
    functional_profile_df = pd.read_csv(functional_profile_filepath, index_col=0)
    label_file_df = label_file_df[functional_profile_df.columns].transpose()
    sample_to_label = label_file_df[0].to_dict()

    # Get test, training and validation labels.
    dataset_separation_filepath = os.path.join(output_folder, 'Dataset_separation', 'Annotation_samples_separation_Iteration_0.csv')
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


def test_run_iterate_reference_test_sets_filepath():
    functional_profile_filepath = os.path.join('input', 'test_functional_profile.csv')
    label_filepath = os.path.join('input', 'test_label.csv')
    reference_test_sets_filepath = os.path.join('input', 'test_reference_sets.csv')
    output_folder = 'output_folder'
    run_nb = 1
    nb_iterations = 1
    nb_classifiers = 1
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    run_iterate(functional_profile_filepath, label_filepath, output_folder, run_nb, nb_iterations, classifiers=nb_classifiers, reference_test_sets_filepath=reference_test_sets_filepath)

    expected_test_set = ['MH0037', 'MH0036', 'MH0021', 'MH0011', 'MH0009']
    expected_training_set = ['MH0003', 'MH0033', 'MH0032', 'MH0012', 'MH0038', 'MH0028', 'MH0031', 'MH0039', 'MH0024', 'MH0016', 'MH0030', 'MH0020', 'MH0026', 'MH0006']
    expected_validation_set = ['MH0014', 'MH0002', 'MH0025', 'MH0034', 'MH0035']

    # Get label for sample.
    label_file_df = pd.read_csv(label_filepath)
    functional_profile_df = pd.read_csv(functional_profile_filepath, index_col=0)
    label_file_df = label_file_df[functional_profile_df.columns].transpose()
    sample_to_label = label_file_df[0].to_dict()

    # Get test, training and validation labels.
    dataset_separation_filepath = os.path.join(output_folder, 'Dataset_separation', 'Annotation_samples_separation_Iteration_0.csv')
    dataset_separation = pd.read_csv(dataset_separation_filepath, index_col=0).to_dict()

    training_set = dataset_separation['1']['training_set'].split(',')
    validation_set = dataset_separation['1']['validation_set'].split(',')
    test_set = dataset_separation['1']['test_set'].split(',')

    training_set_label = [str(sample_to_label[sample]) for sample in training_set]
    validation_set_label = [str(sample_to_label[sample]) for sample in validation_set]
    test_set_label = [str(sample_to_label[sample])  for sample in test_set]

    used_training_set_label = [label for label in dataset_separation['1']['training_set_labels'].split(',')]
    used_validation_set_label = [label for label in dataset_separation['1']['validation_set_labels'].split(',')]
    used_test_set_label = [label for label in dataset_separation['1']['test_set_labels'].split(',')]

    # There was an issue where used sample do not correspond to label ID (due to wrong index selection after splitting train/ validation sets).
    # This lead to mismatch between label for the X dataset and the y label.
    # These assertions ensure that both training and validation are performed with matching X dataset and y label.
    assert used_training_set_label == training_set_label
    assert used_validation_set_label == validation_set_label
    assert used_test_set_label == test_set_label

    # Test that expected test, training and validation sets are always the same (for reproducibility).
    assert expected_test_set == test_set
    assert expected_training_set == training_set
    assert expected_validation_set == validation_set

    shutil.rmtree(output_folder)
