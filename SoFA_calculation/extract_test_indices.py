import pandas as pd
import random
import numpy as np
from sklearn.model_selection import train_test_split
import csv
import argparse

parser = argparse.ArgumentParser(description='A program that raandomly extracts a portion of the datasets to be used as benchmarking tests')

parser.add_argument("-d", "--dataset_name", help="Name of the dataset used for this treatment")

parser.add_argument("-p", "--pipeline_path", help="Path to the whole pipeline folder")

args = parser.parse_args()


labels_full = pd.read_csv(args.pipeline_path+'/Inputs/Label_'+args.dataset_name+'.csv', sep=',', index_col=False, header=None)

print(labels_full)


if not labels_full.values.shape[1] > 1:
    label_flatten = labels_full.values.reshape((labels_full.values.shape[0]))
else:
    print('FileSpecificationError: The label file contains more than 1 column.')
    exit()

indices = np.arange(len(label_flatten))

y_train, y_test, sample_train, sample = train_test_split(label_flatten, indices, test_size=0.2, stratify=label_flatten)
print(y_test)

# sample_size = int(len(dataset_full.columns) / 5)
# sample = random.sample(list(dataset_full.columns), sample_size)

with open(args.pipeline_path+'/SoFA_calculation/outputs/'+args.dataset_name+'/test_indices.csv', 'w') as f:
    write = csv.writer(f)

    write.writerow(sample)


# dataset_sampled = dataset_full.iloc[:,sample]
# dataset_opposite = dataset_full.drop(sample, axis=1)

