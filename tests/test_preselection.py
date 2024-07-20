# -*- coding: utf-8 -*-
"""
Created on Sat Jul 20 11:49:20 2024

@author: yannl
"""

import shutil
import pandas as pd
import os
from csv_diff import load_csv, compare
import random
import numpy as np

from SPARTA.iteration import run_iterate
from SPARTA.classification import run_sparta_classification

def test_preselection_run_sparta_classification(seed_init = 0):
    functional_profile_filepath = 'test_functional_profile.csv'
    label_filepath = 'test_label.csv'
    output_folder = 'output_folder' + "_seed" + str(seed_init)
    selected_annots_filepath = "selected_annots.csv"
    run_nb = 3
    nb_iterations = 1
    
    run_sparta_classification(functional_profile_filepath, label_filepath, output_folder, run_nb, nb_iterations, classifiers=2, reference_test_sets_filepath=None, seed_init = seed_init, selected_annots_filepath = selected_annots_filepath)

    return True

seed_init = 57

# functional_profile_filepath = 'test_functional_profile.csv'
# functional_profile_df = pd.read_csv(functional_profile_filepath, sep=',', index_col=0)
# selected_annots_filepath = "selected_annots.csv"
# selected_annots = pd.read_csv(selected_annots_filepath)
# run_nb = 2

# selected_annots_run = selected_annots['Run_'+str(run_nb)].dropna()
# print(selected_annots_run.values)
# shape1 = functional_profile_df.shape
# print(shape1)
# deepmicro_sofa_iteration = functional_profile_df.loc[selected_annots_run.values].transpose()
# shape2 = deepmicro_sofa_iteration.shape
# print(shape2)

test_preselection_run_sparta_classification(seed_init)