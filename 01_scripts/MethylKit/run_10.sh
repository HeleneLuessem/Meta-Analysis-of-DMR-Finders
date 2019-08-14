name="run.sh"
path_to_config_file="../../config.tsv"
path_to_data_file="../../data_prep_10.csv" #"data_test"
path_to_scripts_folder="scripts"
path_to_data_folder="../../02_data"
path_to_results_folder="../../03_results"

printHelp(){
echo -e "" >&2
echo -e "Does DMR calling from starting with bed files to converting the MethylKit outputs to the standard format. Choose which steps to perform" >&2
echo -e "" >&2
echo -e "" >&2
echo -e "Usage: $name <options>" >&2
echo -e "" >&2
echo -e " Mandatory options:" >&2
echo -e "  -i FILE\tFlag to indicate that the input files should be converted to MethylKit input format" >&2
echo -e "  -d FILE\tFlag to indicate that DMRs should be called with MethylKit" >&2
echo -e "  -o FILE\tFlag to indicate that the MethylKit output should be converted to the stnadard outpur format" >&2
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

# (1) Convert bed input files to MethylKit input files
if [ "$convert" = "TRUE" ]; then
	mkdir -p $path_to_data_folder/MethylKit/converted
	ending="_MethylKit_input.tsv"
	while IFS=$',' read -r -a data
	do
		file_name=${data[0]##*/}
		file_name_no_ending=${file_name%.*}
		line=${data[0]}
		if [ ${line:0:1} != "#" ]; 
		then
		        awk -vOFS='\t' 'BEGIN {print "chrBase", "chr", "base", "strand", "coverage", "freqC", "freqT"} {print $1"."$2, $1, $2, $6, $5, $4*100, 100-($4*100)}' ${data[0]} > $path_to_data_folder/MethylKit/converted/$file_name_no_ending$ending
		fi
	done < $path_to_data_file
	echo -e "\n\nMethylKit input data is prepared.\n"	
fi

# (2) Run MethylKit
if [ "$callDMRs" = "TRUE" ]; then
	mkdir -p $path_to_results_folder/MethylKit
	# Add all input files to a string
	input_files=""
	ending="_MethylKit_input.tsv"
	num_A=0
	num_B=0
	while IFS=$',' read -r -a data     
	do
		file_name=${data[0]##*/}
		file_name_no_ending=${file_name%.*}
		line=${data[0]}
		if [ ${line:0:1} != "#" ];
		then
			input_files="$input_files $path_to_data_folder/MethylKit/converted/$file_name_no_ending$ending"
#			echo "${data[1]}"
#			echo -e "$path_to_data_folder/MethylKit/converted/$file_name_no_ending$ending\n\n"
			# Get number of samples for group A and B
			if [ ${data[1]} == "A" ];
			then
				num_A=$((num_A+1))
			fi
			if [ ${data[1]} == "B" ];
                        then
                                num_B=$((num_B+1))
                        fi
		fi
	done < $path_to_data_file
	Rscript $path_to_scripts_folder/MethylKit.R $num_A $num_B $path_to_results_folder/MethylKit/MethylKit_DMRs_raw.tsv $input_files
	echo -e "DMRs are called\n"
fi


# (3) Convert MethylKit output files to standard format TODO: Remove p-value filtering
if [ "$standardize" = "TRUE" ]; then
	# (1) Find p-value cutoff in config file, default is 0.05
	#p_val_cutoff=0.05
	#while IFS=$'\t', read -r -a config
        #do
        #        case "${config[0]}" in
        #                "P-value Cutoff") p_val_cutoff=${config[1]};;
	#	esac
	#done < $path_to_config_file
	# (2) Store result to $path_to_results_folder/MethylKit/RnBeads_DMRs_std.tsv
	#awk -vOFS="\t" -vpval=$p_val_cutoff '$1 != "chr" && $6<=pval{print ($2, $3, $4, "NA", "NA", "NA", $8, $6, $7)}' "$path_to_results_folder/MethylKit/MethylKit_DMRs_raw.tsv" > "$path_to_results_folder/MethylKit/MethylKit_DMRs_std.tsv"
	
	
	awk -vOFS='\t' '{print($2, $3, $4, "NA", "NA", "NA", $8, $6, $7)}' "$path_to_results_folder/MethylKit/MethylKit_DMRs_raw.tsv" > "$path_to_results_folder/MethylKit/MethylKit_DMRs_std.tsv"
	echo -e "Results are now in standard format: Chr Start End #CpGs meanMet1 meanMet2 MetDiff QualityMeasure\n"
fi 
