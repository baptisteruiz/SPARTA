import os
import pandas as pd
import requests
import json

from Bio.ExPASy import Enzyme
from tqdm import tqdm
from goatools import obo_parser

ROOT = os.path.dirname(__file__)
DATA_ROOT = os.path.join(ROOT, 'data')


def create_annotation_mapping_files():
    '''
    This function queries the OBO and ExPASy data banks to get the names of the IDs
    '''
    data_folder = DATA_ROOT

    # Check if we have the ./data directory already
    if(not os.path.isfile(data_folder)):
        # Emulate mkdir -p (no error if folder exists)
        try:
            os.mkdir(data_folder)
        except OSError as e:
            if(e.errno != 17):
                raise e
    else:
        raise Exception('Data path (' + data_folder + ') exists as a file. '
                    'Please rename, remove or change the desired location of the data path.')

    enzyme_dat_file = os.path.join(data_folder, 'enzyme.dat')
    #Check if we have enzyme file, if not download it
    if not(os.path.isfile(enzyme_dat_file)):
        url = "https://ftp.expasy.org/databases/enzyme/enzyme.dat"
        r = requests.get(url, allow_redirects=True)
        open(enzyme_dat_file, 'wb').write(r.content)

    go_basic_file = os.path.join(data_folder, 'go-basic.obo')
    #Same with GO file
    if(not os.path.isfile(go_basic_file)):
        go_obo_url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
        r = requests.get(go_obo_url, allow_redirects=True)
        open(go_basic_file, 'wb').write(r.content)

    annotation_ids = []
    annotation_names = []
    annotation_alternative_ids = []
    go = obo_parser.GODag(go_basic_file)
    for go_id in go:
        go_term = go[go_id]
        annotation_names.append(go_term.name)
        annotation_ids.append(go_id)
        annotation_alternative_ids.append(','.join([alt_id for alt_id in go_term.alt_ids]))

    with open(enzyme_dat_file) as handle:
        enzyme_expasy_version = ''
        for index, line in enumerate(handle):
            if index == 5:
                enzyme_expasy_version = line

    with open(enzyme_dat_file) as handle:
        records = Enzyme.parse(handle)
        for item in records:
            annotation_names.append(item["DE"])
            annotation_ids.append(item["ID"])
            annotation_alternative_ids.append('')

    dataframe = pd.DataFrame(list(zip(annotation_ids, annotation_names, annotation_alternative_ids)), columns=["ID","Name", "Alternative ID"])
    annotation_mapping_filepath = os.path.join(DATA_ROOT, 'mapping_annotation_names.csv')
    dataframe.to_csv(annotation_mapping_filepath, index=False)

    metadata = {}
    metadata['GO version'] = go.version
    metadata['ENZYME Expasy version'] = enzyme_expasy_version
    annotation_medata_filepath = os.path.join(DATA_ROOT, 'mapping_annotation_names.json')
    with open(annotation_medata_filepath, 'w') as ouput_file:
        json.dump(metadata, ouput_file, indent=4)

    os.remove(enzyme_dat_file)
    os.remove(go_basic_file)