
# (1) Sorting
# (2) Merging +,- Strands
# (3) Coverage Filtering
# (4) Removal of extra Chromosomes

name="pre_process.sh"
path_to_config_file="../config.tsv"
path_to_data_file="../data.csv"

path_to_results_folder="../00_inputs"

printHelp(){
echo -e "" >&2
echo -e "Modifies bed files so that the following is given:\n\t(1) Bed files are sorted\n\t(2) + and - strands are merged\n\t(3) Bed files are filtered for sequencing coverage (default 10)\n\t(4) Extra chromosomes are removed" >&2
echo -e "" >&2
echo -e "" >&2
echo -e "Usage: $name <options>" >&2
echo -e "" >&2
echo -e " Mandatory options:" >&2
echo -e "  -m FILE\tFlag to indicate that the input files should not be merged" >&2
echo -e "  -r FILE\tFlag to indicate that the input files should not be freed from extra chromosomes" >&2
echo -e "" >&2
}

toMerge=TRUE
toRemoveChr=TRUE

while getopts "mr" opt
do
 case "$opt" in
  m) toMerge=FALSE;;
  r) toRemoveChr=FALSE;;
 esac
done


# Read in data file
minCov=10
while IFS=$'\t' read -r -a config
do
	case "${config[0]}" in
		"Minimum Read Coverage")        minCov=${config[1]};;
	esac
done < $path_to_config_file

# Pre-process Data
while IFS=$',' read -r -a data
do
	line=${data[0]}
        if [ ${line:0:1} != "#" ]; then
		line_num_before=$(wc -l $line | cut -d' ' -f1)
		ending="_pre_processed.bed"
		file_name=${data[0]##*/}
                file_name_no_ending=${file_name%.*}
		echo $line
		if [ "$toMerge" = "TRUE" ]; 
		then
			if [ "$toRemoveChr" = "TRUE" ];
			then
				# (1) Sort (2) Merge (3) Coverage Filter (4) Chr Filter
                                python merge_CpGs.py -i $line -o temp
				awk -vOFS="\t" -vcov=$minCov '$5 >= cov && ($1 == "chr1" || $1 == "chr2" || $1 == "chr3" || $1 == "chr4" || $1 == "chr5" || $1 == "chr6" || $1 == "chr7" || $1 == "chr8" || $1 == "chr9" || $1 == "chr10" || $1 == "chr11" || $1 == "chr12" || $1 == "chr13" || $1 == "chr14" || $1 == "chr15" || $1 == "chr16" || $1 == "chr17" || $1 == "chr18" || $1 == "chr19" || $1 == "chr20" || $1 == "chr21" || $1 == "chr22" || $1 == "chrX" || $1 == "chrY" || $1 == "1" || $1 == "2" || $1 == "3" || $1 == "4" || $1 == "5" || $1 == "6" || $1 == "7" || $1 == "8" || $1 == "9" || $1 == "10" || $1 == "11" || $1 == "12" || $1 == "13" || $1 == "14" || $1 == "15" || $1 == "16" || $1 == "17" || $1 == "18" || $1 == "19" || $1 == "20" || $1 == "21" || $1 == "22" || $1 == "X" || $1 == "Y")   {print ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11)}' temp | bedtools sort -i  > "$path_to_results_folder/$file_name_no_ending$ending"
				rm temp
				echo -e "(1) Sorted\n(2) Merged\n(3) Coverage Filtered (>= $minCov)\n(4) Chromosomes Filtered\n\n"
				line_num_after=$(wc -l $path_to_results_folder/$file_name_no_ending$ending | cut -d' ' -f1)
				echo -e "# of CpGs: $line_num_before --> $line_num_after"
			else
				# (1) Sort (2) Merge (3) Coverage Filter
				python merge_CpGs.py -i $line -o temp
				awk -vOFS="\t" -vcov=$minCov '$5 >= cov {print ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11)}' temp | bedtools sort -i  > "$path_to_results_folder/$file_name_no_ending$ending"
				rm temp
				echo -e "(1) Sorted\n(2) Merged\n(3) Coverage Filtered (>= $minCov)\n\n"
				line_num_after=$(wc -l $path_to_results_folder/$file_name_no_ending$ending | cut -d' ' -f1)
                                echo -e "# of CpGs: $line_num_before --> $line_num_after"

			fi
		else
			if [ "$toRemoveChr" = "TRUE" ];
			then
				# (1) Sort (3) Coverage Filter (4) Chr Filter
				awk -vOFS="\t" -vcov=$minCov '$5 >= cov && ($1 == "chr1" || $1 == "chr2" || $1 == "chr3" || $1 == "chr4" || $1 == "chr5" || $1 == "chr6" || $1 == "chr7" || $1 == "chr8" || $1 == "chr9" || $1 == "chr10" || $1 == "chr11" || $1 == "chr12" || $1 == "chr13" || $1 == "chr14" || $1 == "chr15" || $1 == "chr16" || $1 == "chr17" || $1 == "chr18" || $1 == "chr19" || $1 == "chr20" || $1 == "chr21" || $1 == "chr22" || $1 == "chrX" || $1 == "chrY" || $1 == "1" || $1 == "2" || $1 == "3" || $1 == "4" || $1 == "5" || $1 == "6" || $1 == "7" || $1 == "8" || $1 == "9" || $1 == "10" || $1 == "11" || $1 == "12" || $1 == "13" || $1 == "14" || $1 == "15" || $1 == "16" || $1 == "17" || $1 == "18" || $1 == "19" || $1 == "20" || $1 == "21" || $1 == "22" || $1 == "X" || $1 == "Y")   {print ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11)}' $line | bedtools sort -i  > "$path_to_results_folder/$file_name_no_ending$ending"
				echo -e "(1) Sorted\n(2) Coverage Filtered (>= $minCov)\n(3) Chromosomes Filtered\n\n"
				line_num_after=$(wc -l $path_to_results_folder/$file_name_no_ending$ending | cut -d' ' -f1)
                                echo -e "# of CpGs: $line_num_before --> $line_num_after"
			else
				# (1) Sort (3) Coverage Filter
				awk -vOFS="\t" -vcov=$minCov '$5 >= cov {print ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11)}' $line | bedtools sort -i  > "$path_to_results_folder/$file_name_no_ending$ending"
				echo -e "(1) Sorted\n(2) Coverage Filtered (>= $minCov)\n\n"
				line_num_after=$(wc -l $path_to_results_folder/$file_name_no_ending$ending | cut -d' ' -f1)
                                echo -e "# of CpGs: $line_num_before --> $line_num_after"

			fi
		fi
	fi
done < $path_to_data_file

