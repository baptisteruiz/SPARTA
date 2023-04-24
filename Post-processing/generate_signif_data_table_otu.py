import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='A program that extracts the significant functional annotations from the original dataset')

parser.add_argument("-d","--dataset_name", help="Name of the dataset used for this treatment")

parser.add_argument("-p", "--pipeline_path", help="Path to the whole pipeline folder")

parser.add_argument("-r", "--data_ref", help="Full reference to the dataset's folder")

parser.add_argument("-i", "--iteration", help="Reference for the method's iteration")

args = parser.parse_args()


database_no_headers = pd.read_csv(args.pipeline_path + '/SoFA_calculation/outputs/'+ args.dataset_name +'/'+args.dataset_name+'_stripped.tsv', sep = '\t')#, header = None)
database_no_headers.columns.values[0] = 'ID'
id_list = pd.read_csv(args.pipeline_path + '/Outputs/' + args.data_ref + '/Selection_outputs/run_'+args.iteration+'/OTU_'+args.dataset_name+'.csv', sep=',', header = 0)
signif_list = id_list.iloc[:id_list.loc[id_list['ID']=='CUTOFF'].index[0],:]
result = pd.merge(signif_list["ID"], database_no_headers, on="ID")

result.to_csv(args.pipeline_path + '/SoFA_calculation/outputs/'+ args.dataset_name +'/'+args.dataset_name+'_iteration_'+str(int(args.iteration)+1)+'_stripped.tsv', sep = '\t', header=True, index=False)
result = result.transpose()
result = result.drop('ID', axis = 0)
result.to_csv(args.pipeline_path + '/SoFA_calculation/outputs/'+ args.dataset_name +'/entree_DeepMicro_'+args.dataset_name+'_iteration_'+str(int(args.iteration)+1)+'_OTU.csv', header=False, index=False)