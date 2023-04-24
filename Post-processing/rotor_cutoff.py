import pandas as pd
from kneebow.rotor import Rotor

import argparse

parser = argparse.ArgumentParser(description='A program that automatically detects and applies a cutoff to a table of ranked feature importances')

parser.add_argument("-d","--dataset_name", help="Name of the dataset used for this treatment")

parser.add_argument("-p", "--pipeline_path", help="Path to the whole pipeline folder")

args = parser.parse_args()

datatable = pd.read_csv(args.pipeline_path + '/Post-processing/outputs/' + args.dataset_name + '.csv', header = 0, index_col=None)


data = datatable["Average_importance"].values

data_rotor = []
for i in range(len(data)):
    data_rotor.append([i+1,data[i]])

rotor = Rotor()
rotor.fit_rotate(data_rotor)

elbow_idx = rotor.get_elbow_index()

datatable.loc[elbow_idx + 0.5] = ["CUTOFF"]*len(datatable.columns)
datatable = datatable.sort_index().reset_index(drop=True)

print(datatable)
print(args.pipeline_path + '/Post-processing/outputs/' + args.dataset_name + '.csv')
datatable.to_csv(args.pipeline_path + '/Post-processing/outputs/' + args.dataset_name + '.csv', index=None)