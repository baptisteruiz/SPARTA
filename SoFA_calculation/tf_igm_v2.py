import pandas as pd
import numpy as np
import math

import argparse

parser = argparse.ArgumentParser(description='A program that scales the calculated SoFAs based on the TF-IGM method')

parser.add_argument("-d", "--dataset_name", help="Name of the dataset used for this treatment")

parser.add_argument("-p", "--pipeline_path", help="Path to the whole pipeline folder")

args = parser.parse_args()

dataset_name = args.dataset_name
dataset_no_iter = args.dataset_name.split('_iteration')[0]
dataset_no_iter = dataset_no_iter.split('_test')[0]
pipeline_path = args.pipeline_path

def tf_igm(dataline, lbd):

    line_sorted = np.array(sorted(dataline, reverse = True))
    # print("line_sorted: ",line_sorted)
    fk1 = line_sorted[0]
    #print("fk1: ",fk1)
    multipliers_array = np.array([i+1 for i in range(len(line_sorted))])
    #print("multipliers_array: ",multipliers_array)
    fk_sum_list = line_sorted * multipliers_array
    #print("fk_sum_list: ",fk_sum_list)
    fk_sum = np.sum(fk_sum_list)
    #print("fk_sum: ",fk_sum)
    
    for data_index in range(len(dataline)):
        if fk_sum != 0:
            dataline[data_index] = math.sqrt(dataline[data_index]) * (1 + lbd * fk1 / fk_sum)
        
        else:
            dataline[data_index] = 0
        #print(dataline[data_index])
    
    return(dataline)


lbd = 7
# dataset = pd.DataFrame([[100,2,0.001],[100,1,0.1],[101,1,0]])
dataset_path = pipeline_path+'/SoFA_calculation/outputs/'+dataset_no_iter+'/entree_DeepMicro_'+dataset_name+'.csv'
dataset = pd.read_csv(dataset_path, header=None, sep = ',')

sum_vector = dataset.sum(axis=1)

freq_matrix = dataset.div(sum_vector, axis=0).fillna(0)

for i in dataset.columns.values:
    #print(i)
    score_test = tf_igm(freq_matrix.loc[:,i].values, lbd)
    dataset[i] = score_test

out_path = pipeline_path+'/DeepMicro/data/entree_DeepMicro_'+dataset_name+'_tfigm.csv'
dataset.to_csv(out_path, header=None, index = False)

#print(max(dataset[10724]))