import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from collections import Counter, defaultdict

def get_median_perfs_and_best_iter(bank_of_performance_dfs, median_classifs_per_iteration):
    '''
    This function calculates the median performances of each run for a given selection iteration, and returns the ones corresponding to the best iteration (i.e: best mean of median AUC performances)
    '''
    best_mean_perf = 0
    best_selec_iter = 0

    # Processing the mean performance of each run per level of iteration
    for iteration_lv in bank_of_performance_dfs.keys():
        iteration_lv_perfs = bank_of_performance_dfs[iteration_lv]

        for run_lv in iteration_lv_perfs.keys():
            run_df = iteration_lv_perfs[run_lv]
            test_perf_values = run_df['Test performance'].values
            median_test_perf_values = np.median(test_perf_values)
            median_classifs_per_iteration[iteration_lv].append(median_test_perf_values)
        
        #Finding and recording the best performing iteration
        mean_perf = np.mean(median_classifs_per_iteration[iteration_lv])
        if mean_perf > best_mean_perf:
            best_mean_perf = mean_perf
            best_selec_iter = iteration_lv

    return median_classifs_per_iteration[best_selec_iter], mean_perf, best_selec_iter

def plot_classifs(bank_of_performance_dfs_annots, bank_of_performance_dfs_taxons, dataset_name, output_file, otu_abundance_filepath=None):
    '''
    This function processes the classification performances obtained through the previous iterative process, and plots the best iteration's performances
    '''

    # Gather the best iteration and the corresponding median performances for each iterative level
    median_classifs_best_iteration_annots = defaultdict(list)
    median_classifs_best_iteration_annots, mean_perf_annots, best_selec_iter_annots = get_median_perfs_and_best_iter(bank_of_performance_dfs_annots, median_classifs_best_iteration_annots)
    best_iteration_indices = [best_selec_iter_annots]

    if otu_abundance_filepath is not None:
        median_classifs_best_iteration_taxons = defaultdict(list)
        median_classifs_best_iteration_taxons, mean_perf_taxons, best_selec_iter_taxons = get_median_perfs_and_best_iter(bank_of_performance_dfs_taxons, median_classifs_best_iteration_taxons)
        best_iteration_indices.append(best_selec_iter_taxons)
        uval, pval = mannwhitneyu(median_classifs_best_iteration_annots, median_classifs_best_iteration_taxons)  
    else:
        median_classifs_best_iteration_taxons = []
        best_selec_iter_taxons = None
    
    # Create a dataframe to be plotted
    performance_scores = median_classifs_best_iteration_annots + median_classifs_best_iteration_taxons
    profiles = ['Functional' for i in range(len(median_classifs_best_iteration_annots))]
    data_name = [dataset_name  for i in range(len(median_classifs_best_iteration_annots))]
    best_run = [best_selec_iter_annots for i in range(len(median_classifs_best_iteration_annots))]
    mean_perfs= [mean_perf_annots for i in range(len(median_classifs_best_iteration_annots))]


    if otu_abundance_filepath is not None:
        profiles = profiles + ['Taxonomic' for i in range(len(median_classifs_best_iteration_taxons))]
        data_name = data_name + [dataset_name  for i in range(len(median_classifs_best_iteration_taxons))]
        best_run = best_run + [best_selec_iter_taxons for i in range(len(median_classifs_best_iteration_taxons))]
        mean_perfs= mean_perfs + [mean_perf_taxons for i in range(len(median_classifs_best_iteration_taxons))]

    dataframe_to_plot = pd.DataFrame.from_dict({'AUC score': performance_scores, 'Profile': profiles, 'Disease': data_name, 'Best run': best_run, 'Mean Perfs': mean_perfs})

    # Plot the info
    fig,ax = plt.subplots(1,constrained_layout=True, figsize=(10,10))

    medians = sns.boxplot(
            y='Disease',
            x='AUC score',
            hue = 'Profile',
            data=dataframe_to_plot,
            ax=ax,
            palette = 'pastel',
            fliersize = 0)
    
    sns.stripplot(y='Disease', x='AUC score', hue = 'Profile', data=dataframe_to_plot, jitter = True, dodge=True, ax=ax, palette = 'dark')

    # Plot presentation
    max_val = max(performance_scores)

    plt.text(max_val + 0.01, -0.5, 'Optimal number\nof selections', size=15, fontweight = 'bold')
    
    for xi,yi,tx, col in zip([max_val + 0.03 for i in range(len(np.unique(profiles)))], (np.arange(len(dataframe_to_plot['Mean Perfs'].values))/2)-0.2, best_iteration_indices, ['blue','darkorange']):
        plt.text( xi, yi,tx, size=18, fontweight = 'bold', c=col)
    
    handles, labels = ax.get_legend_handles_labels()
    l = plt.legend(loc='upper left', borderaxespad=0., fontsize = 15)

    plt.yticks(size=15, rotation=90)
    ax.set_ylabel('Disease',size=20)

    plt.xticks(size=15)
    ax.set_xlabel('Median AUC',size=20)

    title = 'Classification_performances'
    if otu_abundance_filepath is not None:
        title += '\n(Mann-Whitney test p-value: '+str(pval)+')'

    ax.set_title(title)
    plt.savefig(output_file, dpi=300)

    return(best_selec_iter_annots, best_selec_iter_taxons)