#!/bin/bash

#SBATCH --job-name=Colorectal_Construction_fichier_melange_comptages_annot_with_esmecata
#SBATCH --output=output_Colorectal_fichier_melange_comptages_annot.txt

## test command example:
## chmod +x  sparta_main.sh
##./sparta_main.sh -d abundance_tryout -t tf_igm -s relative -i 3 -f 3 -r 3 -e False

# . /local/env/envconda.sh
# conda activate esmecata

scale=None
iterations=1
repeats=1
forests=20
esmecata_run=True

while getopts "d:t:s:i:f:r:e:" arg; do
    case $arg in
        d) dataset_name=$OPTARG;;
        t) treatment=$OPTARG;;
        s) scale=$OPTARG;;
        i) iterations=$OPTARG;;
		f) forests=$OPTARG;;
		r) repeats=$OPTARG;;
		e) esmecata_run=$OPTARG;;
    esac
done


esmecata_plus_check () {

	## EsMeCaTa
  	esmecata proteomes -i $PWD/SoFA_calculation/outputs/$dataset_name/$dataset_name.tsv -o $PWD/SoFA_calculation/outputs/$dataset_name/esmecata_outputs --remove-tmp

	esmecata clustering -i $PWD/SoFA_calculation/outputs/$dataset_name/esmecata_outputs -o $PWD/SoFA_calculation/outputs/$dataset_name/esmecata_outputs_clustering --remove-tmp

	esmecata annotation -i $PWD/SoFA_calculation/outputs/$dataset_name/esmecata_outputs_clustering -o $PWD/SoFA_calculation/outputs/$dataset_name/esmecata_outputs_annot 
		
	## Check
	
	while IFS=$'\t' read -r rec_column1 rec_column2; do
  		if [ ! -f "$PWD/SoFA_calculation/outputs/$dataset_name/esmecata_outputs_annot/annotation_reference/${rec_column1}.tsv" ]
  			then
  				echo "${rec_column1} : EsMeCaTa go brrrrrrrrrrrrrrrrrrrrr"
  				esmecata_plus_check
  		fi; 
	done < <(tail -n +2 $PWD/SoFA_calculation/outputs/$dataset_name/$dataset_name.tsv)	
	
	## Clean
	
	rm -r $PWD/SoFA_calculation/outputs/$dataset_name/esmecata_outputs
	rm -r $PWD/SoFA_calculation/outputs/$dataset_name/esmecata_outputs_clustering
}



if [ $dataset_name == "" ]
	then
    		echo "Please enter a name for the dataset"
    
elif ([ "$treatment" == "scaling" ] || [ "$treatment" == "tf_igm" ] || [ "$treatment" == "" ]) && ( [ "$scale" == "relative" ] || [ "$scale" == "None" ])
	then
		path=$PWD
		
		
		######## Creating correct output folder
		data_ref="${dataset_name}"

		if [ "$treatment" != "" ]
			then
				data_ref="${dataset_name}_${treatment}"
			fi

		if [ "$scale" != "None" ]
			then
			
			data_ref="${data_ref}_${scale}"
			fi

		if [ ! -d "$PWD/Meta_Outputs" ]
			then
				mkdir $PWD/Meta_Outputs
			fi
		
		if [  -d "$PWD/Meta_Outputs/${data_ref}" ]
			then
				rm -r $PWD/Meta_Outputs/${data_ref}
			fi
		mkdir $PWD/Meta_Outputs/${data_ref}

		for (( repeat_nb=1; repeat_nb<=$repeats; repeat_nb++ )); do

			it=1
		
			if [ ! -d "$PWD/Outputs" ]
				then
					mkdir $PWD/Outputs
				fi
			
			if [ ! -d "$PWD/Outputs/$data_ref" ]
				then
					mkdir $PWD/Outputs/$data_ref
				fi
			########## STEP 1: Build annotations scores #############
			
			if [ ! -d "$PWD/Outputs/$data_ref/SoFAs" ]
				then
					mkdir $PWD/Outputs/$data_ref/SoFAs
				fi


			#python -c 'from ete3 import NCBITaxa;ncbi = NCBITaxa();ncbi.update_taxonomy_database()'

			if [ ! -d "$PWD/SoFA_calculation/outputs" ]
				then
					mkdir $PWD/SoFA_calculation/outputs
				fi
			
			if [ ! -d "$PWD/SoFA_calculation/outputs/$dataset_name" ]
				then
					mkdir $PWD/SoFA_calculation/outputs/$dataset_name
				fi	
			
			echo 'Formatting entries...'
			python $PWD/SoFA_calculation/data_formatting_OTUrename.py -d "$dataset_name" -p "$PWD"

			if [ "$scale" == "relative" ]
				then
					python SoFA_calculation/relative_transform.py -d "$dataset_name" -p "$PWD"
				fi


			if [ "$esmecata_run" == True ]
				then
				esmecata_plus_check
			
			else
				if [ ! -d $PWD/SoFA_calculation/outputs/${dataset_name}/esmecata_outputs_annot ]
					then
						mkdir $PWD/SoFA_calculation/outputs/${dataset_name}/esmecata_outputs_annot
					fi
				
				if [ -d $PWD/SoFA_calculation/outputs/${dataset_name}/esmecata_outputs_annot/annotation_reference ]
					then
						rm -r $PWD/SoFA_calculation/outputs/${dataset_name}/esmecata_outputs_annot/annotation_reference

					fi

				cp -r $PWD/Inputs/annotation_reference_${dataset_name} $PWD/SoFA_calculation/outputs/${dataset_name}/esmecata_outputs_annot/annotation_reference

				fi
			
			if [ ! -d $PWD/DeepMicro/data ]
				then
					mkdir $PWD/DeepMicro/data
				fi

			echo 'Building functional scores...'
			python SoFA_calculation/create_score_annot_depending_on_otus.py -d "$dataset_name" -p "$PWD" -s "$scale"
			
			echo 'Separate test and train...'
			
			python SoFA_calculation/extract_test_indices.py -d "$dataset_name" -p "$PWD"
			
			
			if [ "$treatment" == "tf_igm" ]
				then
					python SoFA_calculation/tf_igm_v2.py -d "$dataset_name" -p "$PWD"
					python SoFA_calculation/deep_micro_subset_separation.py -d "$dataset_name" -p "$PWD" -r "$data_ref" --tfigm

			elif [ "$treatment" == "scaling" ]
				then
					python SoFA_calculation/scaling_data.py -d "$dataset_name" -p "$PWD"
					python SoFA_calculation/deep_micro_subset_separation.py -d "$dataset_name" -p "$PWD" -r "$data_ref" --scaling
				
			else
				python SoFA_calculation/deep_micro_subset_separation.py -d "$dataset_name" -p "$PWD" -r "$data_ref"
				
				fi
			
			
			
			### Move appropriate files
			
			cp "$PWD/SoFA_calculation/outputs/$dataset_name/${dataset_name}.tsv" "$PWD/Outputs/$data_ref/SoFAs/OTU_name_table_ref.tsv"
			cp "$PWD/SoFA_calculation/outputs/$dataset_name/${dataset_name}_stripped.tsv" "$PWD/Outputs/$data_ref/SoFAs/OTU_table_stripped.tsv"
			cp "$PWD/SoFA_calculation/outputs/$dataset_name/score_annot_selon_otus_${dataset_name}.tsv" "$PWD/Outputs/$data_ref/SoFAs/SoFA_table.tsv"
			
			########## STEP 2: DeepMicro analysis #########
			
			if [ ! -d "$PWD/Outputs/$data_ref/Classification_performances" ]
				then
					mkdir $PWD/Outputs/$data_ref/Classification_performances
				fi
			
			if [ ! -d "$PWD/Outputs/$data_ref/Classification_performances/run_1" ]
				then
					mkdir $PWD/Outputs/$data_ref/Classification_performances/run_1
				fi
			
			cd $PWD/DeepMicro

			if [ "$treatment" == "tf_igm" ]
				then
					python DMmodif_export_test_ver.py -r $forests -cd entree_DeepMicro_${dataset_name}_tfigm.csv -cl Label_$dataset_name.csv -m rf 
					python test_set_measurement.py -d "${dataset_name}" -p "$path" -r "$data_ref" -i "1" --tfigm
					
			elif [ "$treatment" == "scaling" ]
				then
					python DMmodif_export_test_ver.py -r $forests -cd entree_DeepMicro_${dataset_name}_scaled.csv -cl Label_$dataset_name.csv -m rf 
					python test_set_measurement.py -d "${dataset_name}" -p "$path" -r "$data_ref" -i "1" --scaling
			else
					python DMmodif_export_test_ver.py -r $forests -cd entree_DeepMicro_${dataset_name}.csv -cl Label_$dataset_name.csv -m rf 
					python test_set_measurement.py  -d"${dataset_name}" -p "$path" -r "$data_ref" -i "1"
				fi
			
		
			
			python DMmodif_export_test_ver.py -r $forests -cd entree_DeepMicro_${dataset_name}_OTU.csv -cl Label_$dataset_name.csv -m rf 
			python test_set_measurement.py -d "${dataset_name}_OTU" -p "$path" -r "$data_ref" -i "1"

			cd $path
			
			python Post-processing/heatmap_generation.py -d "${dataset_name}" -p "$PWD" -r "$data_ref"
			### Move appropriate files
			
			cp "$PWD/DeepMicro/results/Performance_dataframe_entree_DeepMicro_${dataset_name}_OTU.csv" "$PWD/Outputs/$data_ref/Classification_performances/run_1/OTU_classif_perfs.csv"
			
			if [ "$treatment" == "tf_igm" ]
				then
					cp "$PWD/DeepMicro/results/Performance_dataframe_entree_DeepMicro_${dataset_name}_tfigm.csv" "$PWD/Outputs/$data_ref/Classification_performances/run_1/SoFA_classif_perfs.csv"				
			elif [ "$treatment" == "scaling" ]
				then
					cp "$PWD/DeepMicro/results/Performance_dataframe_entree_DeepMicro_${dataset_name}_scaled.csv" "$PWD/Outputs/$data_ref/Classification_performances/run_1/SoFA_classif_perfs.csv"
			else
					cp "$PWD/DeepMicro/results/Performance_dataframe_entree_DeepMicro_${dataset_name}.csv" "$PWD/Outputs/$data_ref/Classification_performances/run_1/SoFA_classif_perfs.csv"
				fi
				
			###### STEP 3: Importance Score Analysis #######
			
			if [ ! -d "$PWD/Outputs/$data_ref/Selection_outputs" ]
				then
					mkdir $PWD/Outputs/$data_ref/Selection_outputs
				fi
			
			if [ ! -d "$PWD/Outputs/$data_ref/Selection_outputs/run_1" ]
				then
					mkdir $PWD/Outputs/$data_ref/Selection_outputs/run_1
				fi
			
			
			if [ "$treatment" == "tf_igm" ]
				then
					python $PWD/Post-processing/get_param_scores.py -d "$dataset_name" -p "$PWD" --tfigm
			elif [ "$treatment" == "scaling" ]
				then
					python $PWD/Post-processing/get_param_scores.py -d "$dataset_name" -p "$PWD" --scaling
			else
					python $PWD/Post-processing/get_param_scores.py -d "$dataset_name" -p "$PWD"
				fi

			python $PWD/Post-processing/get_param_scores_OTU.py -d "$dataset_name" -p "$PWD"

			if [ "$treatment" == "tf_igm" ]
				then
					python $PWD/Post-processing/rotor_cutoff.py -d "${dataset_name}_tfigm_annots_rf" -p "$PWD" 
			elif [ "$treatment" == "scaling" ]
				then
					python $PWD/Post-processing/rotor_cutoff.py -d "${dataset_name}_scaled_annots_rf" -p "$PWD" 
			else
					python $PWD/Post-processing/rotor_cutoff.py -d "${dataset_name}_annots_rf" -p "$PWD" 
				fi
				
			python $PWD/Post-processing/rotor_cutoff.py -d "${dataset_name}_OTU_rf" -p "$PWD"

			if [ "$treatment" == "tf_igm" ]
				then
					python $PWD/Post-processing/averaging_per_group.py -d "$dataset_name" -p "$PWD" --tfigm
					python $PWD/Post-processing/compare_rf_signif_comparison_sick_ctrl.py -d "$dataset_name" -p "$PWD" -r "$data_ref" -i 1 --tfigm
			elif [ "$treatment" == "scaling" ]
				then	
					python $PWD/Post-processing/averaging_per_group.py -d "$dataset_name" -p "$PWD" --scaling
					python $PWD/Post-processing/compare_rf_signif_comparison_sick_ctrl.py -d "$dataset_name" -p "$PWD" -r "$data_ref" -i 1 --scaling
			else
					python $PWD/Post-processing/averaging_per_group.py -d "$dataset_name" -p "$PWD"
					python $PWD/Post-processing/compare_rf_signif_comparison_sick_ctrl.py -d "$dataset_name" -p "$PWD" -r "$data_ref" -i 1
				fi

			######STEP 4: iteration of the method
			
			while [ $it -lt $iterations ]
			do
				python $PWD/Post-processing/generate_signif_data_table.py -d "$dataset_name" -p "$PWD" -r "$data_ref" -i "$it"
				python $PWD/Post-processing/generate_signif_data_table_otu.py -d "$dataset_name" -p "$PWD" -r "$data_ref" -i "$it"
				((it=it+1))
				
				echo "Iteration: ${it}"
				
				if [ "$treatment" == "tf_igm" ]
					then
						python SoFA_calculation/tf_igm_v2.py -d "${dataset_name}_iteration_${it}" -p "$PWD"
						python SoFA_calculation/deep_micro_subset_separation.py -d "${dataset_name}_iteration_${it}" -p "$PWD" -r "$data_ref" --tfigm

				elif [ "$treatment" == "scaling" ]
					then
						python SoFA_calculation/scaling_data.py -d "${dataset_name}_iteration_${it}" -p "$PWD"
						python SoFA_calculation/deep_micro_subset_separation.py -d "${dataset_name}_iteration_${it}" -p "$PWD" -r "$data_ref" --scaling
				
				else
					python SoFA_calculation/deep_micro_subset_separation.py -d "${dataset_name}_iteration_${it}" -p "$PWD" -r "$data_ref"
				
					fi
				
				
					### Move appropriate files
				# cp "$PWD/SoFA_calculation/outputs/$dataset_name/entree_DeepMicro_${dataset_name}_iteration_${it}.csv" "$PWD/DeepMicro/data"
				# cp "$PWD/SoFA_calculation/outputs/$dataset_name/entree_DeepMicro_${dataset_name}_iteration_${it}_test.csv" "$PWD/DeepMicro/data"
				# cp "$PWD/SoFA_calculation/outputs/$dataset_name/entree_DeepMicro_${dataset_name}_iteration_${it}_OTU.csv" "$PWD/DeepMicro/data"
				# cp "$PWD/SoFA_calculation/outputs/$dataset_name/entree_DeepMicro_${dataset_name}_iteration_${it}_OTU_test.csv" "$PWD/DeepMicro/data"
				
				########## Iterating STEP 2: DeepMicro analysis #########
				cd $PWD/DeepMicro
				if [ "$treatment" == "tf_igm" ]
					then
						python DMmodif_export_test_ver.py -r $forests -cd entree_DeepMicro_${dataset_name}_iteration_${it}_tfigm.csv -cl Label_$dataset_name.csv -m rf 
						python test_set_measurement.py -d "${dataset_name}_iteration_${it}" -p "$path" -r "$data_ref" -i "$it" --tfigm
				elif [ "$treatment" == "scaling" ]
					then
						python DMmodif_export_test_ver.py -r $forests -cd entree_DeepMicro_${dataset_name}_iteration_${it}_scaled.csv -cl Label_$dataset_name.csv -m rf 
						python test_set_measurement.py -d "${dataset_name}_iteration_${it}" -p "$path" -r "$data_ref" -i "$it" --scaling
				else
						python DMmodif_export_test_ver.py -r $forests -cd entree_DeepMicro_${dataset_name}_iteration_${it}.csv -cl Label_$dataset_name.csv -m rf 
						python test_set_measurement.py -d "${dataset_name}_iteration_${it}" -p "$path" -r "$data_ref" -i "$it"
					fi
					
				python DMmodif_export_test_ver.py -r $forests -cd entree_DeepMicro_${dataset_name}_iteration_${it}_OTU.csv -cl Label_$dataset_name.csv -m rf 
				python test_set_measurement.py -d "${dataset_name}_iteration_${it}_OTU" -p "$path" -r "$data_ref" -i "$it"

				cd $path

				python Post-processing/heatmap_generation.py -d "${dataset_name}_iteration_${it}" -p "$PWD" -r "$data_ref"
				
				if [ ! -d "$PWD/Outputs/$data_ref/Classification_performances/run_${it}" ]
				then
					mkdir $PWD/Outputs/$data_ref/Classification_performances/run_${it}
				fi
			
				cp "$PWD/DeepMicro/results/Performance_dataframe_entree_DeepMicro_${dataset_name}_iteration_${it}_OTU.csv" "$PWD/Outputs/$data_ref/Classification_performances/run_${it}/OTU_classif_perfs.csv"
				
				if [ "$treatment" == "tf_igm" ]
					then
						cp "$PWD/DeepMicro/results/Performance_dataframe_entree_DeepMicro_${dataset_name}_iteration_${it}_tfigm.csv" "$PWD/Outputs/$data_ref/Classification_performances/run_${it}/SoFA_classif_perfs.csv"				
				elif [ "$treatment" == "scaling" ]
					then
						cp "$PWD/DeepMicro/results/Performance_dataframe_entree_DeepMicro_${dataset_name}_iteration_${it}_scaled.csv" "$PWD/Outputs/$data_ref/Classification_performances/run_${it}/SoFA_classif_perfs.csv"
				else
						cp "$PWD/DeepMicro/results/Performance_dataframe_entree_DeepMicro_${dataset_name}_iteration_${it}.csv" "$PWD/Outputs/$data_ref/Classification_performances/run_${it}/SoFA_classif_perfs.csv"
					fi
				
				###### Iterating STEP 3: Importance Score Analysis #######
				
				if [ ! -d "$PWD/Outputs/$data_ref/Selection_outputs/run_${it}" ]
					then
						mkdir $PWD/Outputs/$data_ref/Selection_outputs/run_${it}
					fi
				
				if [ "$treatment" == "tf_igm" ]
					then
						python $PWD/Post-processing/get_param_scores.py -d "${dataset_name}_iteration_${it}" -p "$PWD" --tfigm
				elif [ "$treatment" == "scaling" ]
					then
						python $PWD/Post-processing/get_param_scores.py -d "${dataset_name}_iteration_${it}" -p "$PWD" --scaling
				else
						python $PWD/Post-processing/get_param_scores.py -d "${dataset_name}_iteration_${it}" -p "$PWD"
					fi

				python $PWD/Post-processing/get_param_scores_OTU.py -d "${dataset_name}_iteration_${it}" -p "$PWD"
			
			
				if [ "$treatment" == "tf_igm" ]
					then
						python $PWD/Post-processing/rotor_cutoff.py -d "${dataset_name}_iteration_${it}_tfigm_annots_rf" -p "$PWD"
				elif [ "$treatment" == "scaling" ]
					then
						python $PWD/Post-processing/rotor_cutoff.py -d "${dataset_name}_iteration_${it}_scaled_annots_rf" -p "$PWD"
				else
						python $PWD/Post-processing/rotor_cutoff.py -d "${dataset_name}_iteration_${it}_annots_rf" -p "$PWD" 
					fi
					
				python $PWD/Post-processing/rotor_cutoff.py -d "${dataset_name}_iteration_${it}_OTU_rf" -p "$PWD"
				
				
				if [ "$treatment" == "tf_igm" ]
					then
						python $PWD/Post-processing/compare_rf_signif_comparison_sick_ctrl.py -d "${dataset_name}_iteration_${it}" -p "$PWD" -r "$data_ref" -i "$it" --tfigm
				elif [ "$treatment" == "scaling" ]
					then	
						python $PWD/Post-processing/compare_rf_signif_comparison_sick_ctrl.py -d "${dataset_name}_iteration_${it}" -p "$PWD" -r "$data_ref" -i "$it" --scaling
				else
						python $PWD/Post-processing/compare_rf_signif_comparison_sick_ctrl.py -d "${dataset_name}_iteration_${it}" -p "$PWD" -r "$data_ref" -i "$it"
					fi
			done

			mv $PWD/Outputs/$data_ref $PWD/Meta_Outputs/${data_ref}/${data_ref}_${repeat_nb}
			cp $PWD/SoFA_calculation/outputs/${dataset_name}/test_indices.csv $PWD/Meta_Outputs/${data_ref}/${data_ref}_${repeat_nb}/test_indices.csv
			

		done
	
	python $PWD/Post-processing/plot_otu_vs_sofa.py -d "${dataset_name}" -p "$PWD" -r "$data_ref" -n "$repeats" -i "$iterations"
	python $PWD/Post-processing/export_core_and_meta.py -d "${dataset_name}" -p "$PWD" -r "$data_ref" -n "$repeats" -i "$iterations"
	rm -r $PWD/Outputs
	
				
else
	echo "Please enter a valid value for the arguments:"
	echo "-t: scaling, tf_igm"
	echo "-s: relative"
	fi



