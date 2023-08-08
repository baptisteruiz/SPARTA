# SPARTA (Shifting Paradigms to Annotation Representation from Taxonomy to identify Archetypes)
### Note: please refrain from using '\_test' in the input datasets' names, as it would be confusing during the process

List of necessary Python dependencies:

\*: will not be maintained in future updates, but is necessary in the pipeline's current format

  - esmecata
  - pandas
  - numpy
  - sklearn
  - unittest
  - matplotlib
  - tensorflow*
  - keras*
  - importlib
  - joblib
  - seaborn
  - progress
  - goatools
  - biopython
  - requests
  - operator
  - kneebow
  - argparse




#### SETUP:
The following input files must be added to the "Inputs" folder:

  - OTU relative abundances:
    - Must be a .txt file with tabular separation
    - metadata can also be included within the file, it will not be taken account of
    - the OTU names should be in the format: "k__kingdom|p__phylum|c__class|o__order|f__family|g__genus|s__species"
  - Labels:
    - Must be a .csv file named according to the corresponding OTU abundance file (i.e: if the OTU data is named 'abundance_tryout.txt', the associated labels must be named 'Label_abundance_tryout.csv')
    - Must be shaped like a vector of '0's and '1's, designating which group each individual from the dataset belongs to.
  - EsMeCaTa annotation outputs (only if option e is False):
    - Must be the 'annotation_reference' output folder of EsMeCaTa's 'annotation' workflow
    - The folder must be renamed into the following format: 'annotation_reference_abundance_tryout'

  
  
 #### RUNNING THE SCRIPT:
 Run the run_abundance_to_score_v6.sh script, with the following arguments:
  - Mandatory:
    - d: name of the dataset
  - Optional:
    - t: data treatment (can be: 'tf_igm' or 'scaling', default: no treatment)
    - s: scaling method (can be: 'relative', default: no scaling)
    - i: number of iterations of the method (default: 1 iteration)
    - f: amount of trained classifiers per iteration of the command (default: 20 forests)
    - r: amount of pipeline runs (default: 1 run)
    - e: launch EsMeCaTa within the pipeline (default: True)
 
 To launch a quick (~1hr) test, use command:

    ./sparta_main.sh -d abundance_tryout -t tf_igm -s relative -i 3 -f 3 -r 3 -e False
    
 #### OUTPUT:
 A temporary "Outputs" folder will be generated during the process.
 The "Meta_Outputs" directory will be automatically created if it does not already exist. It will be organised as such:
 
 ````
Meta_Outputs:
    └── X_options (X being the name of the original OTU count file, and options being the -t and -s optional arguments given at launch):
        └── X_options_n (for n in the amount of pipeline repetitions)
          ├──  Classification_performances
          |     ├── Data_leak_free_X_performances.txt: Performance scores (AUC) of the best performing RF model on the independent test SoFA datasets
          |     ├── Data_leak_free_X_OTU_performances.txt: Performance scores (AUC) of the best performing RF model on the independent test OTU datasets
          |     └── run_i (for all i iterations of the method
          |           ├── OTU_classif_perfs.txt: Performance metrics of each RF model trained during this iteration on the internal verification subset of the OTU data
          |           └── SoFA_classif_perfs.txt: Performance metrics of each RF model trained during this iteration on the internal verification subset of the SoFA data 
          ├── Selection_outputs
          |     └── run_i (for all i iterations of the method
          |           ├── OTU_X.csv : List of all OTUs ranked by decreasing importance, with the importance cutoff integrated
          |           ├── scores_X.csv : List of all functional annotations ranked by decreasing importance, with the importance cutoff integrated
          |           ├── X_output_comparison_signif_rf.txt : temporary: useless data that I'm going to get rid of eventually...
          |           └── X_options_correlation_plots.png : Visual representation of each feature's correlation score with its counterparts
          └── SoFAs
                ├── OTU_table_stripped.tsv : Original OTU count table with new OTU names and without metadata
                ├── SoFA_table.tsv : SoFA profiles calculated for each sample
                └── OTU_name_table-ref.tsv : list of the new normalised OTU names, and their phylogenetic correspondance

 ````

#### REPRODUCTION:
The original datasets and labels are given in the default "Inputs" folder. To run a reproduction of the paper's experiments, run the command:

  ./sparta_main.sh -d dataset -t tf_igm -i 5 -r 10
  
for each of the included datasets ('abundance_IBD', 'abundance_Obesity', 'abundance_Cirrhosis', 'abundance_Colorectal', 'abundance_T2D', 'abundance_WT2D')
