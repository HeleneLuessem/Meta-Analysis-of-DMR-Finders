name="run.sh"
path_to_config_file="BS_config.tsv"
path_to_data_file="BS_data_prep.csv" #_prep.csv"
path_to_scripts_folder="01_scripts/BSmooth/scripts"
path_to_data_folder="02_data"
path_to_results_folder="03_results"

printHelp(){
echo -e "" >&2
echo -e "Does DMR calling from starting with bed files to converting the BSmooth outputs to the standard format. Choose which steps to perform" >&2
echo -e "" >&2
echo -e "" >&2
echo -e "Usage: $name <options>" >&2
echo -e "" >&2
echo -e " Mandatory options:" >&2
echo -e "  -i FILE\tFlag to indicate that the input files should be converted to BSmooth input format" >&2
echo -e "  -o FILE\tFlag to indicate that DMRs should be called with BSmooth" >&2
echo -e "  -o FILE\tFlag to indicate that the BSmooth output should be converted to the stnadard outpur format" >&2
echo -e "" >&2
}

convert=FALSE
callDMRs=FALSE
standardize=FALSE

while getopts "ido" opt
do
	case "$opt" in
		i) convert=TRUE;;
		d) callDMRs=TRUE;;
		o) standardize=TRUE;;
	esac
done

# (1) Convert bed input files to BSmooth input files
if [ "$convert" = "TRUE" ]; then
	#mkdir -p $path_to_data_folder/BSmooth/converted
        rm -f $path_to_data_folder/BSmooth/group_A.txt
        rm -f $path_to_data_folder/BSmooth/group_B.txt
        touch $path_to_data_folder/BSmooth/group_A.txt
        touch $path_to_data_folder/BSmooth/group_B.txt
	ending="_BSmooth_input.tsv"
	while IFS=$',' read -r -a data
	do
		file_name=${data[0]##*/}
		file_name_no_ending=${file_name%.*}
		line=${data[0]}
		if [ ${line:0:1} != "#" ]; 
		then
			# Create files listing the group's inputs
			if [ "${data[1]}" = "B" ]; then
			#	echo "$path_to_data_folder/BSmooth/converted/$file_name_no_ending$ending" >> $path_to_data_folder/BSmooth/group_B.txt
				echo $line >> $path_to_data_folder/BSmooth/group_B.txt
			else
				if [ "${data[1]}" = "A" ]; then
			#		echo "$path_to_data_folder/BSmooth/converted/$file_name_no_ending$ending" >> $path_to_data_folder/BSmooth/group_A.txt
					echo $line >> $path_to_data_folder/BSmooth/group_A.txt
				else
					echo "Error: Group Name different from A and B!"
					exit 1
				fi
			fi
		fi
	done < $path_to_data_file
	echo -e "\n\nBSmooth input data is prepared.\n"	
fi

# (2) Run BSmooth
if [ "$callDMRs" = "TRUE" ]; then
	mkdir -p $path_to_results_folder/BSmooth
	time Rscript $path_to_scripts_folder/BSmooth.R $path_to_data_folder/BSmooth/group_A.txt $path_to_data_folder/BSmooth/group_B.txt $path_to_results_folder/BSmooth/BSmooth_DMRs_raw.tsv
	echo -e "DMRs are called\n" 
fi


# (3) Convert BSmooth output files to standard format
if [ "$standardize" = "TRUE" ]; then
	# TODO (1) Store result to $path_to_results_folder/BSmooth/RnBeads_DMRs_std.tsv
	awk -vOFS='\t' '$1 != "start"{print ($2, $3, $4, $8, $14, $15, $13, $11, $12, $16)}'  $path_to_results_folder/BSmooth/BSmooth_DMRs_raw.tsv >  $path_to_results_folder/BSmooth/BSmooth_DMRs_std.tsv
	echo -e "Results are now in standard format: Chr Start End #CpGs meanMet1 meanMet2 MetDiff QualityMeasure\n"
fi 

touch BSmooth_dummy.txt
