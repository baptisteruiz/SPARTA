import pandas as pd
import numpy as np
from collections import defaultdict

def create_core_and_meta_dfs(iteration_selections_per_run, runs):
    '''
    This function creates, from a dictionary containing selected variables per iteration level and run:
        - a 'Core' dataframe containing the features that are present in all of the selections for a given iteration
        - a 'Meta' dataframe containing the features that are present in at least one of the selections for a given iteration, along with the number of runs in which it is selected
    '''
    
    full_list = sum([list(iteration_selections_per_run[i]) for i in iteration_selections_per_run.keys()],[])
    core_meta = {'core':[], 'meta':{'ID':[],'Count':[]}}

    for feature in np.unique(full_list):
        if full_list.count(feature) == runs:
            core_meta['core'].append(feature)
        else:
            core_meta['meta']['ID'].append(feature)
            core_meta['meta']['Count'].append(full_list.count(feature))
    
    core = pd.DataFrame(core_meta['core'], columns=['ID'])
    meta = pd.DataFrame(core_meta['meta']).sort_values(by='Count', ascending=False)

    return(core,meta)


def formatting_core_meta_outputs(info_df, core_df, meta_df, runs, zero_case=False):
    '''
    This function adds information inherited from upstream analyses to the core and meta dataframes
    '''
    
    if meta_df is None:
        meta_skip = True
    else:
        meta_skip = False
    if zero_case:
        core_info = info_df[info_df['ID'].isin(list(core_df.index))]
    else:
        core_info = info_df[info_df['ID'].isin(list(core_df['ID'].values))]

    if meta_skip:
        meta_info = None
        
    
    else:
        significance_category = []
        for count in meta_df['Count'].values:
            if count > runs*0.75:
                significance_category.append('Confident')
            else:
                significance_category.append('Candidate')

        meta_info = info_df[info_df['ID'].isin(list(meta_df['ID'].values))]
        meta_info['Significance_count'] = meta_df['Count'].values
        meta_info['Significance_category'] = significance_category

    return core_info, meta_info


def extract_core_associates(dataframe, core_list, esmecata_input=None):
    '''
    This function automatically extracts the core-significant from the lists of associated variables
    '''

    if 'Linked_annotations' in dataframe.columns:
        col_ref = 'Linked_annotations'
        new_col = 'Core_linked_annotaions'
        otu_links = False
    else:
        col_ref = 'Linked_OTUs'
        new_col = 'Significant_linked_OTUs'
        otu_links = True

    signif_vars = []
    for vars_assoc in dataframe[col_ref].values:
        vars_list = []
        for var in vars_assoc:
            if var in core_list:
                vars_list.append(var)
        signif_vars.append(vars_list)
    
    dataframe[new_col] = signif_vars

    if otu_links:
        signif_links_named = []
        for otu_list in signif_vars:
            named_links = []
            for otu in otu_list:
                otu_name_translated = esmecata_input[esmecata_input['observation_name'] == otu]['taxonomic_affiliation'].values[0]
                otu_name_translated_species = otu_name_translated.split(';')[-1]
                named_links.append(otu_name_translated_species)
            signif_links_named.append(named_links)
        dataframe['Named_significant_linked_OTUs'] = signif_links_named

    return dataframe


def extract_and_write_core_meta(path_core_meta, bank_of_selections_annots, bank_of_selections_taxons, bank_of_performance_dfs_annots, bank_of_performance_dfs_taxons, best_selec_iter_annots,
                                best_selec_iter_taxons, info_annots, info_taxons, runs, esmecata_input, otu_table_stripped, sofa_table, otu_abundance_filepath):
    '''
    This function launches the process to identify and write the Core and Meta taxons and annotations from all of the effectuated selections. It also records the 
    length of each selection category and classification performances per iteration, to give an overview ofthis information.
    '''

    df_perfs_and_selection_per_iter = defaultdict(defaultdict)
    warning_annots=False
    warning_taxons=False

    for iteration in bank_of_selections_annots.keys():
        
        iteration_selections_per_run_annots = bank_of_selections_annots[iteration]
        core_annots, meta_annots = create_core_and_meta_dfs(iteration_selections_per_run_annots, runs)
        core_annot_info, meta_annot_info = formatting_core_meta_outputs(info_annots, core_annots, meta_annots, runs, zero_case=False)

        if iteration == best_selec_iter_annots - 1:
            core_annots_opti, meta_annots_opti = core_annot_info, meta_annot_info


        if otu_abundance_filepath is not None:
            iteration_selections_per_run_taxons = bank_of_selections_taxons[iteration]
            core_taxons, meta_taxons = create_core_and_meta_dfs(iteration_selections_per_run_taxons, runs)
            core_taxons_info, meta_taxons_info = formatting_core_meta_outputs(info_taxons, core_taxons, meta_taxons, runs, zero_case=False)

            if iteration == best_selec_iter_taxons -1:
                core_taxons_opti, meta_taxons_opti = core_taxons_info, meta_taxons_info
            
            core_taxons_id = list(core_taxons['ID'].values)
            core_annots_id = list(core_annots['ID'].values)
            
            core_annot_info = extract_core_associates(core_annot_info, core_taxons_id, esmecata_input)
            core_taxons_info = extract_core_associates(core_taxons_info, core_annots_id, esmecata_input)

            meta_annot_info = extract_core_associates(meta_annot_info, core_taxons_id, esmecata_input)
            meta_taxons_info = extract_core_associates(meta_taxons_info, core_annots_id, esmecata_input)

            core_taxons_info.to_csv(path_core_meta+'/All_iterations/Core_taxons_iteration_'+str(iteration)+'.csv')
            meta_taxons_info.to_csv(path_core_meta+'/All_iterations/Meta_taxons_iteration_'+str(iteration)+'.csv')

            df_perfs_and_selection_per_iter['Iteration '+str(iteration)+'_Taxonomic']['Robust '] = core_taxons_info.shape[0]
            df_perfs_and_selection_per_iter['Iteration '+str(iteration)+'_Taxonomic']['Confident'] = core_taxons_info.shape[0] + meta_taxons_info[meta_taxons_info['Significance_category'] == 'Confident'].shape[0]
            df_perfs_and_selection_per_iter['Iteration '+str(iteration)+'_Taxonomic']['Candidate'] = core_taxons_info.shape[0] + meta_taxons_info.shape[0]


            iteration_median_perfs = []
            iteration_lv_perfs = bank_of_performance_dfs_taxons[iteration]

            for run_lv in iteration_lv_perfs.keys():
                run_df = iteration_lv_perfs[run_lv]
                test_perf_values = run_df['Test performance'].values
                median_test_perf_values = np.median(test_perf_values)
                iteration_median_perfs.append(median_test_perf_values)
            
            df_perfs_and_selection_per_iter['Iteration '+str(iteration)+'_Taxonomic']['Mean performance (AUC)'] = np.mean(iteration_median_perfs)

            if iteration == best_selec_iter_taxons and np.mean(iteration_median_perfs)<0.6:
                warning_taxons = True

        
        core_annot_info.to_csv(path_core_meta+'/All_iterations/Core_annots_iteration_'+str(iteration)+'.csv')
        meta_annot_info.to_csv(path_core_meta+'/All_iterations/Meta_annots_iteration_'+str(iteration)+'.csv')

        df_perfs_and_selection_per_iter['Iteration '+str(iteration)+'_Functional']['Robust '] = core_annot_info.shape[0]
        df_perfs_and_selection_per_iter['Iteration '+str(iteration)+'_Functional']['Confident'] = core_annot_info.shape[0] + meta_annot_info[meta_annot_info['Significance_category'] == 'Confident'].shape[0]
        df_perfs_and_selection_per_iter['Iteration '+str(iteration)+'_Functional']['Candidate'] = core_annot_info.shape[0] + meta_annot_info.shape[0]

        iteration_median_perfs = []
        iteration_lv_perfs = bank_of_performance_dfs_annots[iteration]

        for run_lv in iteration_lv_perfs.keys():
            run_df = iteration_lv_perfs[run_lv]
            test_perf_values = run_df['Test performance'].values
            median_test_perf_values = np.median(test_perf_values)
            iteration_median_perfs.append(median_test_perf_values)
        
        df_perfs_and_selection_per_iter['Iteration '+str(iteration)+'_Functional']['Mean performance (AUC)'] = np.mean(iteration_median_perfs)

        if iteration == best_selec_iter_annots and np.mean(iteration_median_perfs)<0.6:
            warning_annots = True


    if best_selec_iter_annots == 0:
        core_annots_opti, meta_annots_opti = formatting_core_meta_outputs(info_annots, sofa_table, None, runs, zero_case=True)

    if otu_abundance_filepath is not None:
        core_annots_opti_id = list(core_annots_opti['ID'].values)

        if best_selec_iter_taxons == 0:
            core_taxons_opti, meta_taxons_opti = formatting_core_meta_outputs(info_taxons, otu_table_stripped, None, runs, zero_case=True)
            core_taxons_opti_id = list(core_taxons_opti['ID'].values)
        

        else:
            core_taxons_opti_id = list(core_taxons_opti['ID'].values)
            meta_taxons_opti = extract_core_associates(meta_taxons_opti, core_annots_opti_id, esmecata_input)
        
        if best_selec_iter_annots != 0:
            meta_annots_opti = extract_core_associates(meta_annots_opti, core_taxons_opti_id, esmecata_input)
        
        
        core_annots_opti = extract_core_associates(core_annots_opti, core_taxons_opti_id, esmecata_input)
        core_taxons_opti = extract_core_associates(core_taxons_opti, core_annots_opti_id, esmecata_input)

        

        core_taxons_opti.to_csv(path_core_meta+'/Best_iteration/Core_taxons_iteration_'+str(best_selec_iter_taxons-1)+'.csv')
        if not (meta_taxons_opti is None):
            meta_taxons_opti.to_csv(path_core_meta+'/Best_iteration/Meta_taxons_iteration_'+str(best_selec_iter_taxons-1)+'.csv')

    core_annots_opti.to_csv(path_core_meta+'/Best_iteration/Core_annots_iteration_'+str(best_selec_iter_annots-1)+'.csv')
    if not (meta_annots_opti is None):
        meta_annots_opti.to_csv(path_core_meta+'/Best_iteration/Meta_annots_iteration_'+str(best_selec_iter_annots-1)+'.csv')
        
        
    return(df_perfs_and_selection_per_iter, warning_annots, warning_taxons)