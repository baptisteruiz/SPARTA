import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict


import argparse

parser = argparse.ArgumentParser(description='A program that averages the RF importance scores of the functional annotations, and associates them to OTUs')

parser.add_argument("-d","--dataset_name", help="Name of the dataset used for this treatment")

parser.add_argument("-p", "--pipeline_path", help="Path to the whole pipeline folder")

parser.add_argument("-r", "--data_ref", help="Full reference to the dataset's folder")

parser.add_argument("-n", "--nb_repeats", help="Number of times the pipeline has been repeated")

parser.add_argument("-i", "--nb_iterations", help="Number of iterations per run of the pipeline")

parser.add_argument("--tfigm", action="store_true", help="Indicates whether the pipeline focuses on TF-IGM recalculated data")
parser.add_argument("--scaling", action="store_true", help="Indicates whether the pipeline focuses on scaled data")

args = parser.parse_args()

dataset_name = args.dataset_name

pipeline_path = args.pipeline_path


###COMPILING INFO ABOUT ALL RUNS

median_perf_values_dict = {}

median_otu_test_perfs = defaultdict(list)
median_sofa_test_perfs = defaultdict(list)


for run_nb in range(1, int(args.nb_repeats) + 1):
    
    path_save = pipeline_path + '/Meta_Outputs/'+args.data_ref
    path_ref = path_save +'/'+args.data_ref+'_'+str(run_nb)

    for i in range(int(args.nb_iterations)):

        local_path = path_ref + '/Classification_performances/run_'+str(i+1)+'/SoFA_classif_perfs.csv'   

        dataframe_sofa = pd.read_csv(local_path, sep = ',')
        median_sofa_test_perfs[i].append(np.median(list(dataframe_sofa['Test performance'].values)))

        local_path_otu = path_ref + '/Classification_performances/run_'+str(i+1)+'/OTU_classif_perfs.csv'
        dataframe_otu = pd.read_csv(local_path_otu, sep = ',')
        median_otu_test_perfs[i].append(np.median(list(dataframe_otu['Test performance'].values)))
    
    median_perf_values_dict[dataset_name] = [median_sofa_test_perfs, median_otu_test_perfs]

###PLOTTING MEDIAN INFO

ticks = ['Original Data', 'First Selection', 'Second Selection', 'Third Selection', 'Fourth Selection']

sns.set_theme()
sns.set_context('paper')
sns.set_style('whitegrid')

fig,ax = plt.subplots(1,constrained_layout=True, figsize=(10,10))

best_runs = []
max_vals = []
for i, disease_name in enumerate(median_perf_values_dict.keys()):
    
    
    sofa_test_perfs = median_perf_values_dict[disease_name][0]

    otu_test_perfs = median_perf_values_dict[disease_name][1]
    
    median_sofas_perfs_tfigm = [np.mean(sofa_test_perfs[i]) for i in sofa_test_perfs.keys()]
    median_otu_perfs = [np.mean(otu_test_perfs[i]) for i in otu_test_perfs.keys()]
    

    max_ind_sofas = median_sofas_perfs_tfigm.index(max(median_sofas_perfs_tfigm))
    max_ind_otus = median_otu_perfs.index(max(median_otu_perfs))

    best_runs+=[max_ind_sofas] + [max_ind_otus]

    max_vals.append(max(sofa_test_perfs[max_ind_sofas]))
    max_vals.append(max(otu_test_perfs[max_ind_otus]))

    dataframe = pd.DataFrame.from_dict({'AUC score': sofa_test_perfs[max_ind_sofas]+otu_test_perfs[max_ind_otus]})
    dataframe['Profile']=['SoFAs' for i in range(len(sofa_test_perfs))]+['OTUs' for i in range(len(otu_test_perfs))]
    dataframe['Disease']=[disease_name.split('_')[1] for i in range(len(sofa_test_perfs)+len(otu_test_perfs))]
    dataframe['Best run']=[max_ind_sofas for i in range(len(sofa_test_perfs))] + [max_ind_otus for i in range(len(otu_test_perfs))]
    dataframe['Median']=[median_sofas_perfs_tfigm[max_ind_sofas] for i in range(len(sofa_test_perfs))] + [median_otu_perfs[max_ind_otus] for i in range(len(otu_test_perfs))]
    
    if i == 0:
        df_cumul = dataframe
    else:
        df_cumul = pd.concat([df_cumul, dataframe])

medians = sns.boxplot(
            y='Disease',
            x='AUC score',
            hue = 'Profile',
            data=df_cumul,
            ax=ax,
            palette = 'pastel',
            fliersize = 0)
sns.stripplot(y='Disease', x='AUC score', hue = 'Profile', data=df_cumul, jitter = True, dodge=True, ax=ax, palette = 'dark')



for xi,yi,tx, col in zip([1.05 for i in max_vals], (np.arange(len(dataframe['Median'].values))/2)-0.2, best_runs, ['blue','darkorange']*6):
     plt.text( xi, yi,tx, size=18, fontweight = 'bold', c=col)
# plt.text(1.01, -0.8, 'Optimal selection', size=10, fontweight = 'bold')

handles, labels = ax.get_legend_handles_labels()
l = plt.legend(handles[2:], labels[2:], loc='upper left', borderaxespad=0., fontsize = 15)

plt.yticks(size=15)
ax.set_ylabel('Disease',size=20)

plt.xticks(size=15)
ax.set_xlabel('Median AUC',size=20)

plt.savefig(path_save +'/median_OTU_vs_SoFA_(best_vs_best).png', dpi=300)

