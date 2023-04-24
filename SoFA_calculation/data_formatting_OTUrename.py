import pandas as pd
import os
import argparse
pd.options.mode.chained_assignment = None  # default='warn'


parser = argparse.ArgumentParser()

parser.add_argument("-d", "--dataset_name", help="Name of the dataset used for this treatment")

parser.add_argument("-p", "--pipeline_path", help="Path to the whole pipeline folder")

args = parser.parse_args()


#Import the original file's description of the OTU taxonomies

path = args.pipeline_path + "/Inputs/"+args.dataset_name+".txt"
dataset_full = pd.read_csv(path, header = None, sep = "\t")
dataset_compos = dataset_full.loc[dataset_full[0].str.contains('k__')]


#Replace the separators to create compatibility with EsMeCaTa
new_cols = list(dataset_compos[0].values)

new_cols = [line.replace('k__','') for line in new_cols]
new_cols = [line.replace('|p__',';') for line in new_cols]
new_cols = [line.replace("|c__",";") for line in new_cols]
new_cols = [line.replace("|o__",";") for line in new_cols]
new_cols = [line.replace("|f__",";") for line in new_cols]
new_cols = [line.replace("|g__",";") for line in new_cols]
new_cols = [line.replace("|s__",";") for line in new_cols]
new_cols = [line.replace("_"," ") for line in new_cols]




#Formatting the dataset for DeepMicro and EsMeCaTa, and giving the OTUs standardised aliases (OTU_x)
otu_names = ['OTU_' + str(i+1) for i in range(len(new_cols))]
formatted_data = pd.DataFrame({'observation_name':otu_names, 'taxonomic_affiliation':new_cols}, index = None)

dataset_compos[0] = otu_names
dataset_compos.columns=list(dataset_full[dataset_full[0] == 'sampleID'].values)
dataset_compos = dataset_compos.rename(columns = {'sampleID':'OTU'})

if not os.path.exists(args.pipeline_path+"/SoFA_calculation/outputs"):
    os.makedirs(args.pipeline_path+"/SoFA_calculation/outputs")

if not os.path.exists(args.pipeline_path+"/SoFA_calculation/outputs/"+args.dataset_name):
    os.makedirs(args.pipeline_path+"/SoFA_calculation/outputs/"+args.dataset_name)


#Writing EsMeCaTa input
formatted_data.to_csv(args.pipeline_path+"/SoFA_calculation/outputs/"+args.dataset_name+"/"+args.dataset_name+".tsv", sep = "\t", index = None)
#Writing standardised version of the dataset
dataset_compos.to_csv(args.pipeline_path+"/SoFA_calculation/outputs/"+args.dataset_name+"/"+args.dataset_name+"_stripped.tsv", sep = "\t", index = None)

dataset_compos = dataset_compos.transpose()
dataset_compos.columns = dataset_compos.iloc[0]
dataset_compos = dataset_compos[1:]

#Writing DeepMicro input
dataset_compos.to_csv(args.pipeline_path+"/SoFA_calculation/outputs/"+args.dataset_name+"/entree_DeepMicro_"+args.dataset_name+"_OTU.csv", sep = ",", header = False, index=False)