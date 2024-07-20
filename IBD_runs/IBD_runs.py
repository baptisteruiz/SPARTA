# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 06:53:59 2024

@author: yannl
"""

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


def IBD_runs(run_number, seed_init = 0):
    functional_profile_filepath = 'SoFA_table.csv'
    label_filepath = 'Label_abundance_IBD.csv'
    output_folder = 'output_folder' + "_seed" + str(seed_init) + "_run" + str(run_number)
    nb_iterations = 2
  
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    run_iterate(functional_profile_filepath, label_filepath, output_folder, run_number, nb_iterations, classifiers=20 , reference_test_sets_filepath="Test_sets_abundance_IBD.csv", seed_init = seed_init)

    # Get label for sample.
    label_file_df = pd.read_csv(label_filepath)
    functional_profile_df = pd.read_csv(functional_profile_filepath, index_col=0)
    label_file_df = label_file_df[functional_profile_df.columns].transpose()
    sample_to_label = label_file_df[0].to_dict()

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



    
    # shutil.rmtree(output_folder)
# functional_profile_filepath = 'SoFA_table.csv'
# label_filepath = 'Label_abundance_IBD.csv'
# # functional_profile_filepath = 'test_functional_profile.csv'
# # label_filepath = 'test_label.csv'
# functional_profile_df = pd.read_csv(functional_profile_filepath, sep=',', index_col=0)
# label_file_df = pd.read_csv(label_filepath)
# print(label_file_df)
# print("------\n")
# print(functional_profile_df.columns)
#label_file_df = label_file_df[functional_profile_df.columns]
seed_init = 12 # 11 with previous batch
random.seed(seed_init)
nruns = 2
for i in range(nruns):
    IBD_runs(i+1,random.sample(range(1000),1)[0])