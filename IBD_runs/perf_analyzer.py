# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 11:20:49 2024

@author: yannl
"""

import shutil
import pandas as pd
import os


def perf_analyzer():
    print("Analyzing runs")
    archive_path = ".\\tests_sets_from_paper\\"
    prefixed = [filename for filename in os.listdir(archive_path) if filename.startswith("IBD-output_folder_seed")]
    # index = 0
    # index1 = 1
    # iter_analyzed = "Run_" + str(index) + "\\"
    subpath = "Classification_performances\Iteration_"
    perf0 = pd.DataFrame()
    perf1 = pd.DataFrame()
    nrun_analyzed = 0
    for file in prefixed:
        print("Analyzing folder " + file)
        runs = [filename for filename in os.listdir(archive_path + file) if filename.startswith("Run")]
        for run in runs:
            data_classification_filepath0 = os.path.join(archive_path,file, run, subpath+ str(0) + "\\", 'Annotation_performances.csv')
            data_classification_filepath1 = os.path.join(archive_path,file, run, subpath+ str(1) + "\\", 'Annotation_performances.csv')
            if os.path.isfile(data_classification_filepath0) & os.path.isfile(data_classification_filepath1):
                data_classification = pd.read_csv(data_classification_filepath0, index_col=0).to_dict()
                perf0[file + "_" + run] = None
                perf0[file + "_" + run] = data_classification["Test performance"]
                data_classification = pd.read_csv(data_classification_filepath1, index_col=0).to_dict()
                perf1[file + "_" + run] = None
                perf1[file + "_" + run] = data_classification["Test performance"]
                nrun_analyzed=1+nrun_analyzed
    print("Done, analyzed " + str(nrun_analyzed) + " runs")
    perf0.to_csv("perf0.csv")    
    perf1.to_csv("perf1.csv")  
        
    
perf_analyzer()