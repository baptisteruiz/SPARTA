# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 20:38:35 2024

@author: yannl
"""

import shutil
import pandas as pd
import os
import random

from SPARTA.iteration import run_iterate
from SPARTA.classification import run_sparta_classification


def test_seeding_run_iterate(seed_init=0):
    functional_profile_filepath = 'test_functional_profile.csv'
    label_filepath = 'test_label.csv'
    output_folder = 'output_folder' + "_seed" + str(seed_init)
    output_folder_verif = output_folder + "_verif"
    run_nb = 2
    nb_iterations = 1

    random.seed(seed_init)
    master_seeds = random.sample(range(1000), run_nb)
    random.seed(master_seeds[run_nb-1])
    seed_valid = random.sample(range(1000),1)
    seed_split = random.sample(range(1000),1)
    seed_rf = random.sample(range(1000),1)

    if not os.path.exists(output_folder_verif):
        os.mkdir(output_folder_verif)
        run_iterate(functional_profile_filepath, label_filepath, output_folder_verif, run_nb, nb_iterations, classifiers=2, reference_test_sets_filepath='test_reference_sets.csv',
                    seed_rf=seed_rf[0], seed_split=seed_split[0], seed_valid=seed_valid[0])

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    run_iterate(functional_profile_filepath, label_filepath, output_folder, run_nb, nb_iterations, classifiers=2, reference_test_sets_filepath='test_reference_sets.csv',
                seed_rf=seed_rf[0], seed_split=seed_split[0], seed_valid=seed_valid[0])

    # Expected sets with seed_init = 0
    # expected_test_set = ['MH0037', 'MH0036', 'MH0021', 'MH0011', 'MH0009']
    # expected_training_set = ['MH0028','MH0006','MH0014','MH0030','MH0032','MH0038','MH0012','MH0033','MH0020','MH0031','MH0003','MH0025','MH0026','MH0016']
    # expected_validation_set = ['MH0024','MH0002','MH0034','MH0039','MH0035']

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
    #  assert expected_test_set == test_set
    # assert expected_training_set == training_set
    # assert expected_validation_set == validation_set
    
    # Test whether the classification outputs are the same or not
    index = 0 # maybe try : nb_iterations-1 
    output_filepath = os.path.join(output_folder, 'Classification_performances', 'Iteration_' + str(index), 'Annotation_performances.csv')
    output_df = pd.read_csv(output_filepath)
    expected_output_filpath = os.path.join(output_folder, 'Classification_performances', 'Iteration_' + str(index), 'Annotation_performances.csv')
    expected_df = pd.read_csv(expected_output_filpath)

    assert all(output_df == expected_df)
    
    shutil.rmtree(output_folder)
    shutil.rmtree(output_folder_verif)

def test_seeding_run_sparta_classification(seed_init=0):
    functional_profile_filepath = 'test_functional_profile.csv'
    label_filepath = 'test_label.csv'
    output_folder = 'output_folder' + "_seed" + str(seed_init)

    run_nb = 2
    nb_iterations = 1

    run_sparta_classification(functional_profile_filepath, label_filepath, output_folder, run_nb, nb_iterations, classifiers=2, reference_test_sets_filepath=None, seed_init=seed_init)

    return True
