import pandas as pd
import json
from collections import defaultdict
from io import StringIO

import os

from sparta_pipeline.visualize import get_median_perfs_and_best_iter

def test_median_perfs_and_best_iter():
    
    ### Getting the test inputs and expected outputs
    working_dir = os.getcwd()
    reference_test_inputs_filepath = os.path.join(working_dir, 'input', 'test_get_median_perfs_and_best_iter.json')
    reference_test_outputs_filepath = os.path.join(working_dir, 'expected', 'test_expected_median_perfs_and_best_iter_outputs.json')

    ## Loading inputs
    bank_of_performance_dfs_annots = defaultdict(defaultdict)

    json_pds = json.load(open(reference_test_inputs_filepath))
    for i in range(5):
        for r in range(10):
            bank_of_performance_dfs_annots[i][r] = pd.read_json(StringIO(json_pds[str(i)][str(r)]))

    ## Loading outputs
    expected_outputs = json.load(open(reference_test_outputs_filepath))

    ### Running the function on test inputs

    median_classifs_best_iteration_annots = defaultdict(list)
    median_classifs_best_iteration_annots, mean_perf_annots, best_selec_iter_annots = get_median_perfs_and_best_iter(bank_of_performance_dfs_annots, median_classifs_best_iteration_annots)

    ### Checking results

    assert median_classifs_best_iteration_annots == expected_outputs['median_classifs_best_iteration_annots']
    assert mean_perf_annots == expected_outputs['mean_perf_annots']
    assert best_selec_iter_annots == expected_outputs['best_selec_iter_annots']


test_median_perfs_and_best_iter()