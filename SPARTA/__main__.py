import os
import shutil

import pandas as pd

from collections import defaultdict

import logging
from datetime import datetime
import sys
from tqdm import trange
from tqdm.contrib.logging import logging_redirect_tqdm

import argparse

from SPARTA import __version__ as VERSION
from SPARTA.esmecata import run_esmecata
from SPARTA.classification import run_sparta_classification
MESSAGE = '''
A program that averages the RF importance scores of the functional annotations, and associates them to OTUs.
'''

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
"""
logger = logging.getLogger(__name__)

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    with logging_redirect_tqdm():
        for i in trange(9):
            if i == 4:
                logger.info("console logging redirected to `tqdm.write()`")
    # logging restored
"""
def main():
    '''
    This function formats the inputs, formats the output directories, and launches the iterative seection and classification process r times.
    '''
    now_begin = datetime.now()
    date_time = now_begin.strftime("%d%m%Y%H%M")

    parser = argparse.ArgumentParser(
        'sparta',
        description=MESSAGE + ' For specific help on each subcommand use: esmecata {cmd} --help'
    )
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s ' + VERSION + '\n')

    # Pipeline arguments.
    parent_parser_d = argparse.ArgumentParser(add_help=False)
    parent_parser_d.add_argument("-d","--dataset_name", help="Name of the dataset", required=True)
    parent_parser_t = argparse.ArgumentParser(add_help=False)
    parent_parser_t.add_argument("-t", "--treatment", default=None, help="Data treatment for the functional table (can be: 'tf_igm', default: no treatment)")
    parent_parser_s = argparse.ArgumentParser(add_help=False)
    parent_parser_s.add_argument("-s", "--scaling", default=None, help="Scaling method to apply to the taxonomic table (can be: 'relative', default: no scaling)")
    parent_parser_i = argparse.ArgumentParser(add_help=False)
    parent_parser_i.add_argument("-i", "--iterations", default=5, help="Number of iterations of the method (default: 5 iterations)")
    parent_parser_c = argparse.ArgumentParser(add_help=False)
    parent_parser_c.add_argument("-c", "--classifiers", default=20, help="Amount of trained classifiers per iteration of the command (default: 20)")
    parent_parser_r = argparse.ArgumentParser(add_help=False)
    parent_parser_r.add_argument("-r", "--runs", default=10, help="Amount of pipeline runs (default: 10 runs)")
    parent_parser_m = argparse.ArgumentParser(add_help=False)
    parent_parser_m.add_argument("-m", "--method", default="rf", help="Classifying method to be run (default: Random Forest (rf). Can be: svm)")
    parent_parser_v = argparse.ArgumentParser(add_help=False)
    parent_parser_v.add_argument("-v", "--variable_ranking", default="gini", help="Method for Random Forest variable importance ranking (default: gini. Can be: shap)")

    # EsMeCaTa arguments.
    parent_parser_ta = argparse.ArgumentParser(add_help=False)
    parent_parser_ta.add_argument("-fp", "--taxon_abundance", help="Taxonomic profile matrix file (tsv) having samples as columns and organisms as row, organisms as descrbied as taxonomic affilaitions.", required=True)
    parent_parser_e = argparse.ArgumentParser(add_help=False)
    parent_parser_e.add_argument("--eggnog", default=None, help="Path to the eggnog database for the EsMeCaTa pipeline. If not given, the pipeline will be launhed with the 'UniProt' workflow by default.")
    parent_parser_annotations_only = argparse.ArgumentParser(add_help=False)
    parent_parser_annotations_only.add_argument("--annotations_only", default=False, action='store_true', help="This is a flag that signals that the input is a functional table. If True, all steps involving taxonomic tables will be skipped, and SPARTA will iteratively classify and select on the given functional table alone.")
    parent_parser_reference_test_sets = argparse.ArgumentParser(add_help=False)
    parent_parser_reference_test_sets.add_argument("--reference_test_sets", default=False, action='store_true', help="This option allows the user to give their own test sets to be used during classification.")
    parent_parser_esmecata_relaunch = argparse.ArgumentParser(add_help=False)
    parent_parser_esmecata_relaunch.add_argument("--esmecata_relaunch", default=False, action='store_true', help="This option allows the user to force a re-run of the EsMeCaTa pipeline over an already existing output. This is particularly useful if a previous run of the pipeline was botched at this step.")
    parent_parser_keep_temp = argparse.ArgumentParser(add_help=False)
    parent_parser_keep_temp.add_argument("--keep_temp", default=False, action='store_true', help="This option allows the user to keep the contents of the 'Outputs_temp' folder at the end of the run.")
    parent_parser_update_ncbi = argparse.ArgumentParser(add_help=False)
    parent_parser_update_ncbi.add_argument("--update_ncbi", default=False, action='store_true', help="This option allows the user to force an update of the local NCBI database (taxdump.tar.gz).")

    # Classification arguments.
    # Required arguments for functional profile.
    parent_parser_fp = argparse.ArgumentParser(add_help=False)
    parent_parser_fp.add_argument("-fp", "--functional_profile", help="Functional profile matrix file (tsv) having samples as columns and functions as row.", required=True)
    parent_parser_label = argparse.ArgumentParser(add_help=False)
    parent_parser_label.add_argument("-l", "--label", help="Label matrix (tsv) associating samples with label.", required=True)
    parent_parser_o = argparse.ArgumentParser(add_help=False)
    parent_parser_o.add_argument("-o", "--output", help="Output folder.", required=True)
    # Optional arguments for taxon profiles analysis.
    parent_parser_tp = argparse.ArgumentParser(add_help=False)
    parent_parser_tp.add_argument("-tp", "--taxon_profile", help="Taxonomic profile matrix file (tsv) having samples as columns and organisms as row.", required=False, default=None)
    parent_parser_fo = argparse.ArgumentParser(add_help=False)
    parent_parser_fo.add_argument("-fo", "--functional_occurrence", help="Functions associated with organisms.", required=False, default=None)
    parent_parser_ei = argparse.ArgumentParser(add_help=False)
    parent_parser_ei.add_argument("-ei", "--esmecata_input", help="Input files for esmecata.", required=False, default=None)

    # subparsers
    subparsers = parser.add_subparsers(
        title='subcommands',
        description='valid subcommands:',
        dest='cmd')

    pipeline_parser = subparsers.add_parser(
        'pipeline',
        help='Run all SPARTA pipeline with esmecata.',
        parents=[
            parent_parser_e,
            parent_parser_esmecata_relaunch,
            parent_parser_keep_temp, parent_parser_update_ncbi
            ],
        allow_abbrev=False)

    esmecata_parser = subparsers.add_parser(
        'esmecata',
        help='Run functional profile prediction with esmecata.',
        parents=[
            parent_parser_ta, parent_parser_label, parent_parser_o, parent_parser_t, parent_parser_s,
            parent_parser_i, parent_parser_c, parent_parser_r,
            parent_parser_m, parent_parser_v,
            parent_parser_e, parent_parser_annotations_only,
            parent_parser_reference_test_sets, parent_parser_esmecata_relaunch,
            parent_parser_keep_temp, parent_parser_update_ncbi
            ],
        allow_abbrev=False)

    classification_parser = subparsers.add_parser(
        'classification',
        help='Classify functions from functional profile and label files.',
        parents=[
            parent_parser_fp, parent_parser_label, parent_parser_o, parent_parser_t,
            parent_parser_s, parent_parser_i, parent_parser_c,
            parent_parser_r, parent_parser_m, parent_parser_v,
            parent_parser_reference_test_sets, parent_parser_keep_temp,
            parent_parser_fo, parent_parser_ei, parent_parser_tp
            ],
        allow_abbrev=False)

    pd.options.mode.chained_assignment = None

    args_passed = parser.parse_args()
    # If no argument print the help.
    if len(sys.argv) == 1 or len(sys.argv) == 0:
        parser.print_help()
        sys.exit(1)

    if not os.path.exists(args_passed.output):
        os.mkdir(args_passed.output)
    formatter = logging.Formatter('%(message)s')
    log_file_path = os.path.join(args_passed.output, f'sparta_{args_passed.cmd}.log')
    file_handler = logging.FileHandler(log_file_path, 'w+')
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    # set up the default console logger
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    if args_passed.cmd == 'pipeline':
        dataset_name = args_passed.dataset_name
        pipeline_path = os.getcwd()
        runs = int(args_passed.runs)
        iterations = int(args_passed.iterations)

    if args_passed.method != 'rf':
        logger.info('Only 1 iteration will be effectuated with classification method '+args_passed.method)
        iterations = 1

    test_set_dict = {}
    bank_of_selections_annots = defaultdict(defaultdict)
    bank_of_selections_taxons = defaultdict(defaultdict)
    bank_of_performance_dfs_annots = defaultdict(defaultdict)
    bank_of_performance_dfs_taxons = defaultdict(defaultdict)

    if args_passed.cmd == 'classification':
        # taxon_profile 'abundance_test_stripped.tsv'
        run_sparta_classification(args_passed.functional_profile, args_passed.label, args_passed.output, args_passed.runs, args_passed.iterations,
                            args_passed.esmecata_input, args_passed.functional_occurrence, args_passed.taxon_profile, args_passed.reference_test_sets,
                            args_passed.classifiers, args_passed.method, args_passed.variable_ranking, args_passed.keep_temp)
    elif args_passed.cmd == 'esmecata':
        run_esmecata(args_passed.label, args_passed.taxon_abundance, args_passed.output, args_passed.scaling, args_passed.esmecata_relaunch, args_passed.eggnog, args_passed.update_ncbi)
    elif args_passed.cmd == 'pipeline':
        functional_profile_path, esmecata_input_path, functional_occurrence_path, otu_table_stripped = run_esmecata(args_passed.label, args_passed.taxon_abundance, args_passed.output,
                                                                                                                    args_passed.scaling, args_passed.esmecata_relaunch, args_passed.eggnog, args_passed.update_ncbi)

        run_sparta_classification(functional_profile_path, args_passed.label, args_passed.output, runs, iterations,
                                esmecata_input_path, functional_occurrence_path, otu_table_stripped, args_passed.reference_test_sets,
                            args_passed.classifiers, args_passed.method, args_passed.variable_ranking, args_passed.keep_temp)
