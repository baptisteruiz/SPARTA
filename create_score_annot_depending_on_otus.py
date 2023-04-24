import pandas as pd
import os
import os.path
import numpy as np

from collections import Counter
from progress.bar import IncrementalBar

import argparse

parser = argparse.ArgumentParser(description='A program that builds the scores of each functional annotation based on the prevalence of the OTUs matched with them')

parser.add_argument("-d", "--dataset_name", help="Name of the dataset used for this treatment")

parser.add_argument("-p", "--pipeline_path", help="Path to the whole pipeline folder")

parser.add_argument("-s", "--scaling", help="Indicates the scaling method dictated by the user")

args = parser.parse_args()

dataset_name = args.dataset_name
print("dataset name:", dataset_name)

cwd = os.path.dirname(os.path.abspath(__file__))
os.chdir(cwd)

# Get the number of annotations for each OTU from EsMeCaTa outputs.

esmecata_output_path = args.pipeline_path + '/SoFA_calculation/outputs/'+dataset_name+'/esmecata_outputs_annot/annotation_reference'
otus_GOs = {}
otus_ECs = {}
all_gos = []
all_ecs = []

dir_list = os.listdir(esmecata_output_path)

bar = IncrementalBar('Number of annotations per OTU', max = len(dir_list))

for annot_file in dir_list:
    base_file = os.path.basename(annot_file)
    if base_file[0] != '.':
        base_filename = os.path.splitext(base_file)[0]
        annot_file_path = os.path.join(esmecata_output_path, annot_file)
        df = pd.read_csv(annot_file_path, sep='\t')
        df = df.replace(np.nan, '')
        go_series = df.GO.map(lambda x: [i.strip() for i in x.split(',')]).apply(pd.Series)
        if go_series.empty is False:
            otu_go_terms = list(go_series.stack())
            otu_go_terms = [go for go in otu_go_terms if go != '']
            otu_go_terms_counter = Counter(otu_go_terms)
            otus_GOs[base_filename] = otu_go_terms_counter
            all_gos.extend(otu_go_terms_counter.keys())

        ec_series = df.EC.map(lambda x: [i.strip() for i in x.split(',')]).apply(pd.Series)
        if ec_series.empty is False:
            otu_ec_numbers = list(ec_series.stack())
            otu_ec_numbers = [ec for ec in otu_ec_numbers if ec != '']
            otu_ec_numbers_counter = Counter(otu_ec_numbers)
            otus_ECs[base_filename] = otu_ec_numbers_counter
            all_ecs.extend(otu_ec_numbers_counter.keys())

    bar.next()


all_gos = list(set(all_gos))
all_ecs = list(set(all_ecs))

# Get the abundance of each OTU
if args.scaling == 'relative':
    abundance_sample = pd.read_csv(args.pipeline_path + '/SoFA_calculation/outputs/'+dataset_name+'/'+dataset_name+'_relative.tsv', sep='\t', header = 0, index_col= 0)
else:
    abundance_sample = pd.read_csv(args.pipeline_path + '/SoFA_calculation/outputs/'+dataset_name+'/'+dataset_name+'_stripped.tsv', sep='\t', header = 0, index_col= 0)
#abundance_sample.set_index('sampleID', inplace=True)
all_annots = all_gos + all_ecs
otu_annots = {}

bar2 = IncrementalBar('Checking out OTU abundances', max = len(abundance_sample.index))

for otu in abundance_sample.index:
    otu_annots[otu] = {}
    if otu in otus_ECs:
        otu_annots[otu].update(otus_ECs[otu])
    if otu in otus_GOs:
        otu_annots[otu].update(otus_GOs[otu])
    bar2.next()

# Compute for each annotation its abundance
annot_samples = pd.DataFrame(all_annots)
annot_samples.set_index(0, inplace=True)

bar3 = IncrementalBar('Calculating abundances', max = len(abundance_sample.columns))
for sample in abundance_sample.columns:
    otu_annots_dataframe = pd.DataFrame(otu_annots)
    for col in otu_annots_dataframe.columns: 
        otu_annots_dataframe[col] = otu_annots_dataframe[col] * abundance_sample[sample].loc[col]
    annot_samples[sample] = otu_annots_dataframe.sum(axis=1)
    bar3.next()

annot_samples.to_csv(args.pipeline_path + '/SoFA_calculation/outputs/'+dataset_name+'/score_annot_selon_otus_'+dataset_name+'.tsv', sep='\t')
annot_samples.transpose().to_csv(args.pipeline_path + '/SoFA_calculation/outputs/'+dataset_name+'/entree_DeepMicro_'+dataset_name+'.csv', sep=',', header = None, index = None)