name="run.sh"
path_to_data_file="data.csv"
#path_to_data_file="../../data_head.csv"
ath_to_scripts_folder="01_scripts/RnBeads/scripts"
path_to_data_folder="02_data"
path_to_results_folder="03_results"
path_to_sample_annotation="$path_to_data_folder/RnBeads/sample_annotation.csv"


printHelp(){
echo -e "" >&2
echo -e "Does DMR calling from starting with bed files to converting the Rnbeads outputs to the standard format. Choose which steps to perform" >&2
echo -e "" >&2
echo -e "" >&2
echo -e "Usage: $name <options>" >&2
echo -e "" >&2
echo -e " Mandatory options:" >&2
echo -e "  -i FILE\tFlag to indicate that the input files should be converted to RnBeads input format" >&2
echo -e "  -o FILE\tFlag to indicate that DMRs should be called with Rnbeads" >&2
echo -e "  -o FILE\tFlag to indicate that the RnBeads output should be converted to the stnadard outpur format" >&2
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

# (1) Convert bed input files to RnBeads input files
if [ "$convert" = "TRUE" ]; then		
	# RnBeads input format as given in the BS files
	mkdir -p $path_to_data_folder/RnBeads/converted
	
	rm -f $path_to_sample_annotation
	touch $path_to_sample_annotation
	echo "bedFile,SampleID,Sample_Group" >>  $path_to_sample_annotation
	
	while IFS=$',' read -r -a data
        do
		if [ ${data:0:1} != "#" ]; then
			first=${data[0]}
			group=${data[1]}
			IFS="/"
			array=( $first )
			file_name="${array[${#array[@]}-1]}"
			file_name_no_ending=${file_name%.*}
			ending="_RnBeads_input.bed"
			cp "$first" "$path_to_data_folder/RnBeads/converted/$file_name_no_ending$ending"
			unset IFS
			if [ "$group" = "A" ]; then
				echo "$file_name_no_ending$ending,SS1,$group" >> $path_to_sample_annotation
			else
				echo "$file_name_no_ending$ending,SS2,$group" >> $path_to_sample_annotation

			fi
		fi
	done < $path_to_data_file
	echo -e "\n\nRnbeads input data is prepared.\n"	
fi

# (2) Run RnBeads
if [ "$callDMRs" = "TRUE" ]; then
	mkdir -p $path_to_results_folder/RnBeads
	Rscript $path_to_scripts_folder/RnBeads.R
	# Create Raw Files
	#cp $path_to_results_folder/RnBeads/reports/differential_methylation_data/diffMethTable_region_cmp1_genes.csv $path_to_results_folder/RnBeads/RnBeads_DMRs_genes_raw.tsv
	#cp $path_to_results_folder/RnBeads/reports/differential_methylation_data/diffMethTable_region_cmp1_promoters.csv $path_to_results_folder/RnBeads/RnBeads_DMRs_promoters_raw.tsv
	#cp $path_to_results_folder/RnBeads/reports/differential_methylation_data/diffMethTable_region_cmp1_tiling.csv $path_to_results_folder/RnBeads/RnBeads_DMRs_tiling_raw.tsv
	#cp $path_to_results_folder/RnBeads/reports/differential_methylation_data/diffMethTable_region_cmp1_cpgislands.csv $path_to_results_folder/RnBeads/RnBeads_DMRs_cpgislands_raw.tsv
	#rm -rf $path_to_results_folder/RnBeads/reports
	echo -e "DMRs are called\n"
fi


# (3) Convert Rnbeads output files to standard format
if [ "$standardize" = "TRUE" ]; then
	# (1) Store result to $path_to_results_folder/RnBeads/RnBeads_DMRs_std.tsv
	awk -v OFS='\t' -F "," '$1 != "id"{print $2, $3, $4, $14, $7, $8, $9, $11, $13}' $path_to_results_folder/RnBeads/RnBeads_DMRs_genes_raw.tsv > $path_to_results_folder/RnBeads/RnBeads_DMRs_genes_std.tsv
	awk -v OFS='\t' -F "," '$1 != "id"{print $2, $3, $4, $14, $7, $8, $9, $11, $13}' $path_to_results_folder/RnBeads/RnBeads_DMRs_promoters_raw.tsv > $path_to_results_folder/RnBeads/RnBeads_DMRs_promoters_std.tsv
	awk -v OFS='\t' -F "," '$1 != "id"{print ($2, $3, $4, $12, $5, $6, $7, $9, $11)}' $path_to_results_folder/RnBeads/RnBeads_DMRs_tiling_raw.tsv > $path_to_results_folder/RnBeads/RnBeads_DMRs_tiling_std.tsv
	awk -v OFS='\t' -F "," '$1 != "id"{print ($2, $3, $4, $12, $5, $6, $7, $9, $11)}' $path_to_results_folder/RnBeads/RnBeads_DMRs_cpgislands_raw.tsv > $path_to_results_folder/RnBeads/RnBeads_DMRs_cpgislands_std.tsv
	#cat $path_to_results_folder/RnBeads/RnBeads_DMRs_genes_std.tsv $path_to_results_folder/RnBeads/RnBeads_DMRs_promoters_std.tsv $path_to_results_folder/RnBeads/RnBeads_DMRs_tiling_std.tsv $path_to_results_folder/RnBeads/RnBeads_DMRs_cpgislands_std.tsv > $path_to_results_folder/RnBeads/RnBeads_DMRs_std.tsv
	#rm -rf $path_to_results_folder/RnBeads/reports 
	echo -e "Results are now in standard format: Chr Start End #CpGs meanMet1 meanMet2 MetDiff QualityMeasure\n"
fi

touch RnBeads_dummy.txt
