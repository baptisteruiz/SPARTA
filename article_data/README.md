# REPRODUCTION OF THE ARTICLE RESULTS:

In order to reproduce the results presented in the SPARTA article, the user may run the following commands, using the files given in the article_data directory

## EsMeCaTa:

Using the files included in each 'dataset' directory in the 'esmecata_test' folder, as well as the associated 'annotation_reference' folders for an exact match of the EsMeCaTa outputs :

`sparta esmecata -p abundance_dataset.txt -t tf_igm --esmecata-results local/path/to/annotation_reference -o output_folder_esmecata`

## Classification:

Using the contents of each 'dataset' directory in the 'classification_test' folder:

`sparta classification -fp dataset/dataset_functional_profile.csv -l dataset/dataset_label.csv -ta dataset/dataset_esmecata.tsv -fo dataset/dataset_annotation.tsv -tp dataset/dataset_taxonomic_profile.tsv --reference_test_sets dataset/dataset_test_sets.csv -o output_folder_classification`

## Full pipeline:

Using the contents of each 'dataset' directory in the 'pipeline_test' folder:

`sparta pipeline -p abundance_dataset.txt -l dataset/dataset_label.csv --reference_test_sets dataset/dataset_test_sets.csv -t tf_igm --esmecata-results local/path/to/annotation_reference -o output_folder_pipeline`


