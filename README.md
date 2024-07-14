# SPARTA (Shifting Paradigms to Annotation Representation from Taxonomy to identify Archetypes)

## Table of contents
- [SPARTA](#sparta-)
  - [Table of contents](#table-of-contents)
  - [Installation](#installation)

## Installation

### Requirements

- [pandas](https://pypi.org/project/pandas/): To read the input files and manage the matrix through all pipeline.
- [numpy](https://github.com/numpy/numpy): To manage the matrix through all pipeline.
- [scikit-learn](https://github.com/scikit-learn/scikit-learn): To perform the classification, the search for hyperparameters and the computation of performance.
- [joblib](https://github.com/joblib/joblib): used to save classifiers.
- [kneebow](https://github.com/georg-un/kneebow): To detect the elbow of the variable importance curve to decide a threshold.
- [shap](https://github.com/shap/shap): used as an optional alternative to gini to rank variable importance.
- [tqdm](https://github.com/tqdm/tqdm): used to create progress bar.
- [scipy](https://github.com/scipy/scipy): used to compute statistical tests between iteration and function/organism classifications.
- [matplotlib](https://github.com/matplotlib/matplotlib): used to plot comparison of performance between functions and organism classifications.
- [seaborn](https://github.com/mwaskom/seaborn): used to plot comparison of performance between functions and organism classifications.
- [goatools](https://github.com/tanghaibao/goatools): used to handle Gene Ontology Terms.
- [biopython](https://github.com/biopython/biopython): used to handle Enzyme Commission numbers.
- [requests](https://pypi.org/project/requests/): For the REST queries to download GO and EC databases.
- [esmecata](https://github.com/AuReMe/esmecata): To infer functions occurrences (EC numbers and GO Terms) from taxonomic affiliations.

### Classification installation

SPARTA can be installed using pip:

````sh
git clone https://github.com/baptisteruiz/SPARTA.git
cd SPARTA
pip install -e .
````

This will download all dependencies required for the `classification` subcommand of SPARTA (the main part of SPARTA peforming and the Random Forests and then selecting informative features).
For example, to redo the analysis perform in the article, you can install SPARTA this way and use the provided input files in the INPUTS folder.

### Pipeline and esmecata installation

To use the `pipeline` or `esmecata` subcommands (to create functional profile from taxonomic affiliations), you need to install `mmseqs` and potentially `eggnog-mapper`.  This can be done with `conda`:

```sh
conda install mmseqs2 eggnog-mapper -c conda-forge -c bioconda
```

To work eggnog-mapper requires its database, you have to setup it and install [its database](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.12#storage-requirements), refer to the [setup part of the doc](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.12#setup).

## HOW TO USE:

Run the main.py script from the command line, with the following arguments:
    
### REQUIRED:
    
"-d","--dataset_name": Name of the dataset

### OPTIONAL:

    - "-t", "--treatment" (str, default=None): Data treatment for the functional table (can be: 'tf_igm', default: no treatment)
    - "-s", "--scaling" (str, default=None): Scaling method to apply to the taxonomic table (can be: 'relative', default: no scaling)
    - "-i", "--iterations" (int, default=5): Number of iterations of the method
    - "-f", "--forests" (int, default=20): Amount of trained classifiers per iteration of the command
    - "-r", "--runs" (int, default=10): Amount of pipeline runs
    - "-e", "--esmecata" (bool, default=True): Launch EsMeCaTa within the pipeline
    - "--eggnog" (default=False): Path to the eggnog database for the EsMeCaTa pipeline. If not given, the pipeline will be launched with the 'UniProt' workflow by default.
    - "--annotations_only" (default=False): This is a flag that signals that the input is a functional table. If True, all steps involving taxonomic tables will be skipped, and SPARTA will iteratively classify and select on the given functional table alone.
    - "--reference_test_sets" (default=False): This is a flag that allows the user to give their own test sets to be used during classification.
    - "--esmecata_relaunch" (default=False): This is a flag that allows the user to force a re-run of the EsMeCaTa pipeline over an already existing output. This is particularly useful if a previous run of the pipeline was botched at this step.
    
EXAMPLE:
python main.py -d abundance_testing -s relative -t tf_igm -i 3 -r 4

## INPUTS:
All input files should be stored inside the 'Inputs' folder before launching the pipeline.
### REQUIRED:

    - Abundance profile:
    
        Must be a .txt file with tabular separation
        
        First row must be named 'sampleID', with the subsequent sample names heading each column
        
        If the input is a taxonomic abundance file:
        
            Metadata can also be included within the file, it will not be taken account of
            
            Each row must be named after an OTU, and give its abundance per sample
            
            The OTU names should be in the format: "k__kingdom|p__phylum|c__class|o__order|f__family|g__genus|s__species"
            
        If the input is a functional abundance file:
        
            Please do not include metadata in the file
            
            Each row must be named after an annotation (ID: GO term or EC number), and give its abundance per sample
            
            If the input is functional, please raise the --annotations_only flag. No other profile will be calculated.
            
    - Labels:
    
        Must be a .csv file named according to the corresponding abundance file (i.e: if the abundance data is named 'abundance_tryout.txt', the associated labels must be named 'Label_abundance_tryout.csv')
        
        Must be a vector of the individuals' labels (numerical category descriptor), respecting the order in which said samples are given in the taxonomic abundance table.

### CAN BE REQUIRED:
    
    - Test sets to be used (only if the --reference_test_sets flag is used):
    
        Name the file according to the corresponding input (i.e: if the abundance data is named 'abundance_tryout.txt', the associated test set dataframe must be named 'Test_sets_abundance_tryout.csv')
        
        The 'Test_sets.csv' output file from a previous iteration can be used for reproductibility
        
        Column names: 'Run_i' for i in the required amount of runs (make sure it matches the -i argument) 
        
        List the sample IDs to be used as test subsets for each run in the columns
        
        An index column is required 
        

## OUTPUTS DESCRIPTION:

    Meta_Outputs:
       └── run reference (name of the dataset + transformation arguments + date of beginning):
            └── median_OTU_vs_SoFA_(best_vs_best).png: graphical representation of the classification performances (median ROC AUC per run) at the optimal selective iteration for both taxonomic and functional profiles. Both performance 
            |       distributions are compared statistically by a Mann-Whitney U-test, the p-value of which is given in the figure's title. The optimal selection levels for both profiles are also given.
            └── Test_sets.csv: sample IDs used as test sets for each run of the pipeline.
            └── SoFA_table_[dataset_name].csv: functional profile calculated from the taxonomic input
            └── Run_n (for n in the amount of pipeline repetitions)
            |    ├──  Classification_performances
            |    |     └── Iteration_i (for all i iterations of the method
            |    |           ├── Taxonomic_performances.csv: Performance metrics of each RF model trained during this iteration on the training, validation and test subsets of the taxonomic data, and their optimal found parameters
            |    |           └── Annotation_performances.csv: Performance metrics of each RF model trained during this iteration on the training, validation and test subsets of the functional data, and their optimal found parameters
            |    ├── Selected_outputs
            |    |     └── Iteration_i (for all i iterations of the method
            |    |           ├── Selected_annotations_run_n_iter_i.csv : List of all functional annotations selected at this run and iteration, ranked by decreasing importance
            |    |           └── Selected_taxons_run_n_iter_i.csv : List of all taxons selected at this run and iteration, ranked by decreasing importance
            |    └── Trained Classifiers
            |       └── [dataset_name]_[Functions/OTU]_saved_classifier.joblib : best performing classifiers, as determined by performance on validation set, in this iteration.
            └── Core_and_Meta_outputs
                ├── All_iterations
                |      ├── Core_annots_iteration_i (for i in the amount of iterations of the process per repetition): list of the robust (Core) annotations at selection level i
                |      ├── Meta_annots_iteration_i (for i in the amount of iterations of the process per repetition): list of the candidate and confident (Meta) annotations at selection level i (the 'Significance_count' column gives the amount of classifiers that count the feature as significant)
                |      ├── Core_taxons_iteration_i (for i in the amount of iterations of the process per repetition): list of the robust (Core) taxons at selection level i
                |      └── Meta_taxons_iteration_i (for i in the amount of iterations of the process per repetition): list of the candidate and confident (Meta) taxons at selection level i (the 'Significance_count' column gives the amount of classifiers that count the feature as significant)
                |    For each variable in these sublists, the following extra information is given:
                |         - The name of the taxon/annotation
                |         - The associated counterpart variables (if annotation: the taxons that express it; if taxon: the annotation that it expresses), and among them, those that are 'Robust' at the same level of iteration
                |         - Their average expression in each label category, and the group it is representative of (highest average expression)
                |         - If Meta: the amount of times the variable has been selected over all runs of the pipeline at the given iteration level, and which category ('Candidate' or 'Confident')
                |
                └── Best_iterations
                        └── Same outputs as 'All iterations', but only for the level of selection that gives the best classification performances for the functional and taxonomic profiles

### Other outputs:
    - EsMeCaTa_outputs: outputs of the EsMeCaTa pipeline, only the results of the 'annotation' step are kept for storage efficiency. These outputs can be re-used from one application of the pipeline to a dataset to another.
    - data: enzyme and GO OBO databases, for reference. Only downloaded once.

## STEPS OF THE PIPELINE:

1) Formatting:
   
This step shapes the inputs of the pipeline into a form that is compatible with the next steps of the pipeline, notably those based on external tools (EsMeCaTa and DeepMicro), and removes eventual 
metadata from the input. Importation and handling of the data is made through the pandas library.

Main function: formatting_step(dataset_name, data_ref_output_name, pipeline_path)

    Inputs: dataset_full: the full original abundance table
            dataset_name: the original name of the dataset
            dataset_ref_output_name: the name used as reference for the work in progress (name of the dataset + transformation arguments + date of beginning)
            pipeline path: path to main.py

    Outputs: otu_table_stripped: a metadata-less version of the original abundance table
             esmecata_input: a taxonomic description of all taxons, used to launch the esmecata pipeline afterwards
             deepmicro_otu: a transposed version of otu_table_stripped, compatible as a taxonomic input for DeepMicro (or functional input if '--annotations_only')

    This function also applies the eventual formatting steps to be applied to the taxonomic data (conversion to relative abundances, for example)
    If the '--annotation only' flag was raised, this step will be applied to the functional profile given as input

    Calls:
        pre_formatting(dataset_full, dataset_name, data_ref_output_name, pipeline_path)
            
            Outputs: formatted_data: the original information from the full taxonomic table, with normalized names (OTU_i) and with all metadata removed
                     dataset_compos: data frame associating taxonomic unit references to their taxonomy; compatible input for EsMeCaTa

2) EsMeCaTa:
   
This step recovers the functional annotation of the taxonomic units given as input, through the application of the EsMeCaTa pipeline (https://github.com/AuReMe/esmecata.git).
This step will be skipped entirely if results for a dataset with the same name are found, unless the '--esmecata_relaunch' flag is raised, to avoid redundant calculation. 
It will also not occur if the '--annotation only' flag was raised.

Main function: run_esmecata_plus_check(esmecata_input, pipeline_path, data_ref_output_name, dataset_name, eggnog_path=False)

    Inputs: esmecata_input: a taxonomic description of all taxons (taxonomy described with ';' separators)
            pipeline path: path to main.py
            dataset_name: the original name of the dataset
            dataset_ref_output_name: the name used as reference for the work in progress (name of the dataset + transformation arguments + date of beginning)
            eggnog_path: value of the '--eggnog' option, if given; corresponds to the path to the eggnog database to be used by the eggnog-mapper tool
            

    Outputs: None (hard-written, see the outputs file tree)

    This function calls EsMeCaTa's steps individually, and checks the outputs after the proteome and annotations step (most likely to have issues with an HTTPS connexion error).
    The amount of cores used for clustering and eggnog annotation is determined via the multiprocessing.cpu_count() method. If more than 3 CPUs are available, 2 CPUs less than the available total will be used, 
    to avoid potentially exceeding the machine's computational capacity. If 3 or less CPUs are available, only 1 of them will be used by SPARTA in this step. 
    
    'Proteome': called with default parameters, aside from 'option_bioservices' which is set to 'True' because it generates less instability. Results are checked with the 
        proteome_check() function: for each OTU given as input, do we have a downloaded fasta file (checks the contents of the 'proteomes' folder with the os library, and 
        cross-references with the contents of the 'proteome_tax_id.tsv' file). If the verification fails, this step can be relaunched up to 20 times. If the amount of retries 
        exceeds this number, the operation is aborted and an exception is raised.

    'Clustering': called with default paramaters. Re-running this step over previous results can cause crashes, therefore if an existing incomplete output directory is found, 
        it will be deleted and rewritten.
    
    'Annotation': The option is given here to use the eggnog-mapper tool to perform the annotation. Default parameters are used. Results are checked with the check_annotation() 
        function: check if all OTUs with a reference in the 'consensus_fasta' folder an output written in the 'annotations_reference output folder'. Similarly to the 'Proteomes' step,
        this process can be iterated up to 20 times.

3) Score calculation:
   
This step involves the calculation of a functional profile, from the initial taxonomic abundances and the outputs of EsMeCaTa. The profile is exported both in the form of a DataFrame, presented in a similar 
manner to the original output, and in a transposed form compatible with DeepMicro
If requested, transformation of the functional scores via TF-IGM normalization will be oprtated at this step.

Main function: sofa_calculation(pipeline_path, dataset_name, data_ref_output_name, otu_table_stripped)

    Inputs: pipeline path: path to main.py
            dataset_name: the original name of the dataset
            dataset_ref_output_name: the name used as reference for the work in progress (name of the dataset + transformation arguments + date of beginning)
            otu_table_stripped: a metadata-less version of the original taxonomic abundance table
            

    Outputs: sofa_table: a taxonomic description of the microbian communities described by the original taxonomic abundance table. Abundance of taxonomic units are replaced by scores of functional annotations, 
    for the GO terms and EC numbers associated to the taxonomic profile by EsMeCaTa.
             deepmicro_sofa: a transposed version of sofa_table, compatible as a functional input for DeepMicro
        
5) Averaging scores and fetching information:

In this step, descriptive tables are created for each taxonomic unit and/or functional annotation in our profiles. These descriptions notably tell: the average persence of a feature in each labeled catgory, the 
category that expresses it the most on average, and the total average expression of the feature over the whole profile. For taxons, the detailed taxonomy and the annotations it expresses are fetched. For each 
functional annotation, the name of the annotation is fetched from a database (OBO PURLs for GO terms (http://purl.obolibrary.org/obo/go/go-basic.obo), and the Expasy database for EC numbers 
(https://ftp.expasy.org/databases/enzyme/enzyme.dat)), and a list of all OTUs that express it is gathered.
The created descriptive tables are used for reference in post-processing steps.

Main function: averaging_and_info_step(otu_table_stripped, sofa_table, label_file, esmecata_input, pipeline_path, dataset_name)

    Inputs: otu_table_stripped: a metadata-less version of the original taxonomic abundance table
            sofa_table: a taxonomic description of the microbian communities
            label_file: the original label input
            esmecata_input: a taxonomic description of all taxons (taxonomy described with ';' separators)
            pipeline path: path to main.py
            dataset_name: the original name of the dataset

    Outputs: info_annots: a database of information about functional annotations (average abundance per profile, associated taxonomic units, annotation names)
             info_taxons: a database of information about taxonomic units (average abundance per profile, associated functional annotations, detailed taxonomy)

    Recovery of annotation names involves the downloading of the two aforementioned datasets, which are stored in the 'data' folder. These datasets are interrogated using the following tools:
        - goatools' obo_parser function for GO terms (https://github.com/tanghaibao/goatools)
        - the Enzyme module of Bio.ExPASy (https://biopython.org/docs/1.75/api/Bio.ExPASy.html)

6) Iterated classification and selection:
   
This step is repeated as many times as indicated by the '-r' argument. It consists in classifying individuals based on their functional and, if given, taxonomic profiles. This takes place in three steps:

The following is only done once per run:
A) Setting aside a test set:

    Main functions: 
    - separate_test_train(label_file)

        Inputs: label_file: the original label input

        Outputs: y_train: labels associated to the individuals of the training set
                 y_test: labels associated to the individuals of the test set
                 sample: indices of the individuals selected to be set aside as a test set

        This function selects individuals to be set aside as a test set, using the sklearn library's train_test_split (https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.train_test_split.html) 
        function, with a test_size paramater of 0.2.
    
    - create_test_train_subsets(full_dataset, indices)

        Inputs: full_dataset: dataset to be separated (functional or taxonomic)
                indices: indices of the individuals selected to be set aside as a test set
        
        Outputs: dataset_train: subset of the input dataset not corresponding to the test individuals' indices
                 dataset_test: subset of the input dataset corresponding to the test individuals' indices
        
        This is applied to both the functional and taxonomic profiles.

The following steps are then iterated as many times as required by the value of '-i', for each profile independently.
B) Classification (DeepMicro-based):

    Main function: run_exp(seed, best_auc, threshold_opt, perf_dict, Xtest_ext, Ytest_ext, Y_data, label_data, dataset_name, best_feature_records, data_dir)

        Inputs: seed: random seed value for the classification tasks and the selection of validation subsets. If the process is to be iterated ('-i' > 1), each iteration's seed will be the iteration number 
                    (0,1,2,3...). If not, the seed value is set as 42.
                best_auc: best ROC AUC obtained to date. Set to 0 at the beginning of the iteration.
                threshold_opt: optimal distinction threshold value for the best Random Forest classifier to date, as determined through the Youden J. statistic (doi: 10.1002/sim.2993). Initialized 
                    at 0.
                perf_dict: dictionary containing all classification performances obtained for each forest in this iteration. Initialized as an empty dictionary.
                Xtest_ext: subset of the input dataset corresponding to the test individuals' indices
                Ytest_ext: labels associated to the individuals of the test set
                Y_data: subset of the input dataset not corresponding to the test individuals' indices
                label_data: labels associated to the individuals of the training set
                dataset_name: the original name of the dataset
                best_feature_records: list containing all of the variables' Gini importance scores obtained for each forest in this iteration. Initialized as an empty list.
                data_dir: path to the folder where the best trained classifiers of this iteration will be saved.
        
        Outputs: best_auc: updated value of the best ROC AUC obtained to date
                 threshold_opt: updated value of the optimal distinction threshold value for the best Random Forest classifier to date
                 perf_dict: updated dictionary containing all classification performances obtained for each forest in this iteration
                 best_feature_records: updated list containing all of the variables' Gini importance scores obtained for each forest in this iteration
        
        This function is taken and adapted from the implementation of the DeepMicro tool (https://github.com/minoh0201/DeepMicro), and is repeated as many times within an iteration as demanded by the '-f' argument.
        During this process, a validation subset is set aside randomly from the training subset using the sklearn library's train_test_split 
        (https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.train_test_split.html) function, with a test_size paramater of 0.2, and using the random seed given as input.
        Following this procedure, a Random Forest model is trained, with parameters optimized via sklearn's GridSearchCV method (https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.GridSearchCV.html).
        The Random Forest model is also trained through sklearn's implementation (https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html). The hyper-parameters for GridSearch 
        fine-tuning are as follows:
            - 'n_estimators' in range(100, 1001, 200)
            - 'max_features' in ['sqrt', 'log2']
            - 'min_samples_leaf' in [1, 2, 3, 4, 5]
            - 'criterion': 'gini'
        Each time a model is trained, its performance metrics (ROC AUC score, accuracy, recall, precision and F1 score) on training, validation and test sets are recorded. The model with the best performance on the validation set
        during the iteration is exported. The Gini importance metric of each variable for the trained model are also exported.
        Predictions are made through the application of an optimal classification threshold, determined through the maximization of the Youden J. statistic, to probabilistic predictions by the model. The optimal thresholds found 
        for each model are exported, alongside their performance metrics.

C) Variable selection:

    Main function: inflexion_cutoff(datatable) 

        Inputs: datatable: table containing variables, and their average Gini importance score over all Random Forest models trained this iteration. 

        Outputs: datatable_truncated: a truncated version of the input, retaining only the variables that are above or equal to the selection index, calculated to be set at the Gini importance scores' inflexion 
                    point.
        
        From the variables' Gini importance scores averaged over all iterations of the Random Forest training process, an automatic variable selection is operated. This consists in:
            - ranking the variables by decreasing Gini importance score 
            - calculating the inflexion point of the curve of decreasing importance scores, using kneebow's Rotor and get_elbow_index methods (https://github.com/georg-un/kneebow/tree/master)
            - cutting off all variables below the obtained index

        A selection of the retained variables can then be operated on the taxonomic and functional profiles, to be used as input for step B) as many times as required by the '-i' input.

Classification performances and lists of selected variables for all runs and iterations are passed on to the subsequent steps of the pipeline.

7) Plotting classification performances:

Main function: plot_classifs(bank_of_performance_dfs_annots, bank_of_performance_dfs_taxons, dataset_name, pipeline_path, data_ref_output_name, args_passed)

    Inputs: bank_of_performance_dfs_annots: dictionary containing all classification performances obtained on the functional profile
            bank_of_performance_dfs_taxons: dictionary containing all classification performances obtained on the taxonomic profile
            dataset_name: the original name of the dataset
            pipeline_path: path to main.py
            data_ref_output_name: the name used as reference for the work in progress (name of the dataset + transformation arguments + date of beginning)
    
    Returns: best_selec_iter_annots: the selective iteration that yields the best classification performance for the functional profile
             best_selec_iter_taxons: the selective iteration that yields the best classification performance for the taxonomic profile

    This function processes the classification performances obtained through the previous iterative process, and plots the best iteration's performances. The best iteration is defined as the one that maximises the 
    mean median ROC AUC (i.e: for each run, the median AUC o each iteration is calculated. These median performances are then averaged per iteration level over all runs, and the iteration level that maximises this 
    value is retained as optimal). The optimal selection levels found in this step are passed onward.

    Plotting is done using the matplotlib and seaborn libraries. A statistical comparison of the median performances per run at optimal iteration for each profile is also made at this point, using scipy's 
    implementation of the Mann-Whitney U-test (https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html). The test's p-value is indicated on the plot.

8) Measuring the robustness of variables over all runs:

Main function: extract_and_write_core_meta(path_core_meta, bank_of_selections_annots, bank_of_selections_taxons, best_selec_iter_annots, best_selec_iter_taxons, info_annots, info_taxons, runs, esmecata_input, otu_table_stripped, sofa_table)

    Inputs: path_core_meta: path to the folder in which the robust sublists will be saved
            bank_of_selections_annots: dictionary of annotations selected at every iterative level for each run
            bank_of_selections_taxons: dictionary of taxons selected at every iterative level for each run
            best_selec_iter_annots: the selective iteration that yields the best classification performance for the functional profile
            best_selec_iter_taxons: the selective iteration that yields the best classification performance for the taxonomic profile
            info_annots: a database of information about functional annotations (average abundance per profile, associated taxonomic units, annotation names)
            info_taxons: a database of information about taxonomic units (average abundance per profile, associated functional annotations, detailed taxonomy)
            runs: amount of runs effectuated by the pipeline (given by the '-r' argument)
            esmecata_input: a taxonomic description of all taxons (taxonomy described with ';' separators)
            otu_table_stripped: a metadata-less version of the original abundance table
            sofa_table: a taxonomic description of the microbian communities described by the original taxonomic abundance table.

    Outputs: None (hard-written, see the outputs file tree)

    This function compares the sublists obtained at each level of iteration for each run of the pipeline, and counts the common occurences. Unanimously selected variables are labeled as 'Core' (or 'Robust), and 
    those that appear at least once are labeled 'Meta' ('Confident' if they are selected in 75% of the runs or more, 'Candidate' otherwise). The obtained lists are enriched with information, gathered previously 
    in step 4) notably. A sublist of significant linked counterparts is established for each variable, referencing the linked counterparts (taxons that express an annotation, or annotations expressed by a taxon)
    that are listed as 'Robust' at the same level of iteration. Taxons are converted to their taxonomic name rather than their 'OTU_i' identifier in these lists at this point.
    
    The detailed Core and Meta lists for all iterations are written in the 'All_iterations' output folder. The sublists obtained on the optimal iterations for taxonomic and functional profiles, as defined in step 
    6), are saved separately in the 'Best_iteration' output folder. These lists reference each other when establishing significant linked counterparts (i.e: the significant annotations linked to a taxon in these 
    sublists is 'Robust' at the optimal iteration for the taxonomic profile).
    
