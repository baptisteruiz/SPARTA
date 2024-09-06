import pandas as pd
import os
import shutil

from sparta_pipeline.esmecata import sofa_calculation, run_esmecata

def test_sofa_calculation():
    esmecata_annotation_reference = os.path.join('input', 'annotation_reference')
    output_sofa_table_filepath = os.path.join('computed_sofa_table.csv')
    otu_table_stripped = os.path.join('input', 'test_taxon_profile.tsv')
    otu_table_stripped_df = pd.read_csv(otu_table_stripped, index_col=0, sep='\t')

    sofa_calculation(esmecata_annotation_reference, output_sofa_table_filepath, otu_table_stripped_df, treatment='tf_igm')

    expected_sofa_table = os.path.join('expected', 'test_expected_SoFA_table_abundance_test.csv')
    df_expected = pd.read_csv(expected_sofa_table, index_col=0)
    computed_df = pd.read_csv(output_sofa_table_filepath, index_col=0)

    # Check that both files have same index (annotation IDs).
    assert set(df_expected.index) == set(computed_df.index)
    # Make the comparison on the same index.
    computed_df = computed_df.loc[df_expected.index]
    assert all(df_expected.compare(computed_df))

    os.remove(output_sofa_table_filepath)


def test_esmecata_results():
    esmecata_annotation_reference = os.path.join('input', 'annotation_reference')
    input_file = os.path.join('input', 'test_input_pipeline.txt')
    output_folder = 'output_folder'
    run_esmecata(input_file, output_folder, treatment=None, scaling='no scaling', esmecata_results_path=esmecata_annotation_reference)

    output_sofa_table_filepath = os.path.join('output_folder', 'SoFA_table.csv')
    expected_sofa_table = os.path.join('expected', 'test_expected_SoFA_table_abundance_test.csv')
    df_expected = pd.read_csv(expected_sofa_table, index_col=0)
    computed_df = pd.read_csv(output_sofa_table_filepath, index_col=0)

    # Check that both files have same index (annotation IDs).
    assert set(df_expected.index) == set(computed_df.index)
    # Make the comparison on the same index.
    computed_df = computed_df.loc[df_expected.index]
    assert all(df_expected.compare(computed_df))

    shutil.rmtree(output_folder)


test_sofa_calculation()
test_esmecata_results()
