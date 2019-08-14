# Pre-assuptions for bed Files:
#	* Sorted
#	* Strands (+,-) are Merged
#	* Coverage Filtered
#	* Extra Chromosomes removed

name="run.sh"
path_to_config_file="../../config.tsv"
path_to_data_file="../../data_prep_10.csv"

path_to_scripts_folder="scripts"
path_to_data_folder="../../02_data"
path_to_results_folder="../../03_results"

printHelp(){
echo -e "" >&2
echo -e "Does DMR calling from starting with bed files to converting the Metilene outputs to the standard format. Choose which steps to perform" >&2
echo -e "" >&2
echo -e "" >&2
echo -e "Usage: $name <options>" >&2
echo -e "" >&2
echo -e " Mandatory options:" >&2
echo -e "  -i FILE\tFlag to indicate that the input files should be converted to Metilene input format" >&2
echo -e "  -d FILE\tFlag to indicate that DMRs should be called with Metilene" >&2
echo -e "  -o FILE\tFlag to indicate that the Metilene output should be converted to the stnadard outpur format" >&2
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

# (1) Convert bed input files to Metilene input files
if [ "$convert" = "TRUE" ]; then
	mkdir -p $path_to_data_folder/Metilene/converted	
	ending="_Metilene_input.bed"
	group_A_string=""
	group_B_string=""
	while IFS=$',' read -r -a data
	do
		file_name=${data[0]##*/}
                file_name_no_ending=${file_name%.*}
                line=${data[0]}
		if [ ${line:0:1} != "#" ]; then
			time awk -vOFS="\t" '{print ($1, $2, $3, $4)}' ${data[0]} > "$path_to_data_folder/Metilene/converted/$file_name_no_ending$ending"
			if [ "${data[1]}" = "A" ]; then
				group_A_string="$group_A_string,$path_to_data_folder/Metilene/converted/$file_name_no_ending$ending"
                        else
                                if [ "${data[1]}" = "B" ]; then
                                        group_B_string="$group_B_string,$path_to_data_folder/Metilene/converted/$file_name_no_ending$ending"
                                else
                                        echo "Error: Group Name different from A and B!"
                                        exit 1
                                fi
                        fi
		fi
	done < $path_to_data_file
	# Chop of first comma
	group_A_string="${group_A_string:1}"
        group_B_string="${group_B_string:1}"
	# Run perl script to create Metilene input Matrix
	touch log1
	time perl $path_to_scripts_folder/merge_to_matrix.pl --in1 $group_A_string --in2 $group_B_string --out "$path_to_data_folder/Metilene/converted/metilene_input_matrix.csv" --h1 A --h2 B &> log1
	echo -e "\n\nMetilene input data is prepared.\n"
fi

# (2) Run Metilene
if [ "$callDMRs" = "TRUE" ]; then
	mkdir -p $path_to_results_folder/Metilene
	# Set Parameters to Default
	dis_merge=300
	minCG=10
	min_diff=0.1
	num_threads=1
	group_A="A"
	group_B="B"
	min_nonmissing_g1=80
	min_nonmissing_g2=80
	stringency_filter=0.7
	# Overwrite Parameters
	echo $path_to_config_file
	while IFS=$'\t', read -r -a config
	do
                case "${config[0]}" in
                        "Maximum Distance between two CpGs in one DMR")	dis_merge=${config[1]};;
                        "Minimum Number of CpGs")        		minCG=${config[1]};;
                        "Minimum Mean Meth Difference")        		min_diff=${config[1]};;
                        "Number of Threads")        			num_threads=${config[1]};;
                        "Minimal Number of Non-Missing Values in g1")   min_nonmissing_g1=${config[1]};;
                        "Minimal Number of Non-Missing Values in g2")   min_nonmissing_g2=${config[1]};;
                        "Stringency of Valley Filter")        		stringency_filter=${config[1]};;
                esac
        done < $path_to_config_file

	echo -e "\nPARAMETERS: \n\n"
	echo -e "Maximum Distance between two CpGs in one DMR: $dis_merge\n"
	echo -e "Minimum Number of CpGs: $minCG\n"
        echo -e "Minimum Mean Meth Difference: $min_diff\n"
        echo -e "Number of Threads: $num_threads\n"
        echo -e "Minimal Number of Non-Missing Values in g1: $min_nonmissing_g1\n"
        echo -e "Minimal Number of Non-Missing Values in g2: $min_nonmissing_g2\n"
        echo -e "Stringency of Valley Filter: $stringency_filter\n\n"

	# Run Metilene
	touch log2
	time metilene --maxdist $dis_merge --mincpgs $minCG --minMethDiff $min_diff --threads $num_threads --groupA $group_A --groupB $group_B --minNoA $min_nonmissing_g1 --minNoB $min_nonmissing_g2 --valley $stringency_filter "$path_to_data_folder/Metilene/converted/metilene_input_matrix.csv" > "$path_to_results_folder/Metilene/Metilene_DMRs_raw.tsv" 2> log2
	echo -e "DMRs are called\n"
fi

# (3) Convert Metilene output files to standard format
if [ "$standardize" = "TRUE" ]; then
	awk -vOFS="\t" '{print ($1, $2, $3, $6, $9, $10, $5, $7, $8, $4)}' "$path_to_results_folder/Metilene/Metilene_DMRs_raw.tsv" > "$path_to_results_folder/Metilene/Metilene_DMRs_std.tsv"
	echo -e "Results are now in standard format: Chr Start End #CpGs meanMet1 meanMet2 MetDiff p (2D KS) p (MWU)q-value\n"	
fi 


