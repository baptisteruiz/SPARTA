import pandas as pd
import numpy as np


import argparse

parser = argparse.ArgumentParser(description='A program that builds the scores of each functional annotation based on the prevalence of the OTUs matched with them')

parser.add_argument("-d", "--dataset_name", help="Name of the dataset used for this treatment")

parser.add_argument("-p", "--pipeline_path", help="Path to the whole pipeline folder")

args = parser.parse_args()

dataset_name = args.dataset_name
pipeline_path = args.pipeline_path


def absolute_to_relative(db):
    for col_num in db:
        col_sum = db[col_num].sum(axis = 0)
        db[col_num] = db[col_num].apply(lambda x: x/col_sum)
    return db
# dataset = pd.DataFrame([[100,2,0.001],[100,1,0.1],[101,1,0]])
dataset_path = pipeline_path+'/SoFA_calculation/outputs/'+dataset_name+'/'+dataset_name+'_stripped.tsv'
dataset = pd.read_csv(dataset_path, header=0, index_col= 0, sep = '\t')

 
df_scaled = absolute_to_relative(dataset)
#df_scaled = pd.DataFrame(df_scaled)

out_path = pipeline_path+'/SoFA_calculation/outputs/'+dataset_name+'/'+dataset_name+'_relative.tsv'
print(df_scaled)
df_scaled.to_csv(out_path, sep = '\t')

#print(max(dataset[10724]))