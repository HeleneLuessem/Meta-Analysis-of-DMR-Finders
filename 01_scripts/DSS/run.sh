# Pre-assuptions for bed Files:

#       * Sorted
#       * Strands (+,-) are Merged
#       * Coverage Filtered
#       * Extra Chromosomes removed

name="run.sh"
path_to_config_file="config.tsv"
path_to_data_file="data_prep.csv"
path_to_scripts_folder="01_scripts/DSS/scripts"
path_to_data_folder="02_data"
path_to_results_folder="03_results"

printHelp(){
echo -e "" >&2
echo -e "Does DMR calling from starting with bed files to converting the DSS outputs to the standard format. Choose which steps to perform" >&2
echo -e "" >&2
echo -e "" >&2
echo -e "Usage: $name <options>" >&2
echo -e "" >&2
echo -e " Mandatory options:" >&2
echo -e "  -i FILE\tFlag to indicate that the input files should be converted to DSS input format" >&2
echo -e "  -o FILE\tFlag to indicate that DMRs should be called with DSS" >&2
echo -e "  -o FILE\tFlag to indicate that the DSS output should be converted to the stnadard outpur format" >&2
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

# (1) Convert bed input files to DSS input files
if [ "$convert" = "TRUE" ]; then
	mkdir -p $path_to_data_folder/DSS/converted
	rm -f $path_to_data_folder/DSS/group_A.txt
        rm -f $path_to_data_folder/DSS/group_B.txt
	touch $path_to_data_folder/DSS/group_A.txt
        touch $path_to_data_folder/DSS/group_B.txt
	ending="_DSS_input.tsv"
	while IFS=$',' read -r -a data
	do
		file_name=${data[0]##*/}
		file_name_no_ending=${file_name%.*}
		line=${data[0]}
		if [ ${line:0:1} != "#" ]; then 
			time python $path_to_scripts_folder/convert.py -i ${data[0]} -o "$path_to_data_folder/DSS/converted/$file_name_no_ending$ending" 
			# Create files listing the group's inputs
			if [ "${data[1]}" = "B" ]; then
				echo "$path_to_data_folder/DSS/converted/$file_name_no_ending$ending" >> $path_to_data_folder/DSS/group_B.txt
			else 
				if [ "${data[1]}" = "A" ]; then
					echo "$path_to_data_folder/DSS/converted/$file_name_no_ending$ending" >> $path_to_data_folder/DSS/group_A.txt
				else
					echo "Error: Group Name different from A and B!"
					exit 1
				fi
			fi
		fi
	done < $path_to_data_file
	echo -e "\n\nDSS input data is prepared.\n"
fi

# (2) Run DSS
if [ "$callDMRs" = "TRUE" ]; then
	mkdir -p $path_to_results_folder/DSS
	time Rscript $path_to_scripts_folder/DSS.R $path_to_data_folder/DSS/group_A.txt $path_to_data_folder/DSS/group_B.txt $path_to_results_folder/DSS/DSS_DMRs_raw.tsv
        echo -e "DMRs are called\n"
fi

# (3) Convert DSS output files to standard format
if [ "$standardize" = "TRUE" ]; then
	time awk -vOFS="\t" '$10 < 0 {$10 = $10*(-1)} {print  ($2, $3, $4, $6, $7, $8, $9, $10) }' $path_to_results_folder/DSS/DSS_DMRs_raw.tsv > $path_to_results_folder/DSS/DSS_DMRs_std.tsv
#	rm all.Rout
	echo -e "Results are now in standard format: Chr Start End #CpGs meanMet1 meanMet2 MetDiff QualityMeasure\n"
fi

touch DSS_dummy.txt
