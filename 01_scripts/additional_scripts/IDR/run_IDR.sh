path_to_scripts_folder="scripts"

mkdir -p temp
TMPWD="temp"

# Score (higher for more significant)
# Also non-sig DMRs
# QM at $8

printHelp() {
	echo -e "Usage: `basename $0` INPUTFILE_1 ... INPUTFILE_N"
	echo -e ""
	echo -e " Mandatory:"
	echo -e "  -i INPUTFILE\tFiles contianing DMR outputs. TODO Requirements"
}

while getopts ":hi:m:" opt;
do
	case "$opt" in
		i) inputs+=( "${OPTARG}" ) ;;
		*) printHelp; exit 1 ;;
	esac
done


# (1) Generate merged Regions

cmd="bedtools merge -i <(sort -k1,1 -k2,2n -k3,3n `for s in "${inputs[@]}"
do
	if [[ "$s" == *gz ]]
	then
		echo "<(zcat $s | cut -f 1-3)"
	else
		echo "<(cat $s | cut -f 1-3)"
        fi
done | tr '\n' ' '`) > $TMPWD/mergedPeaks.bed"

eval $cmd



# (2) Convert Scores to Ranks and put behind merged regions per tool

mkdir -p "$TMPWD/inputs_ranks"

echo -e "Location\tQM\tTool" > "$TMPWD/all_tools"

for val in "${inputs[@]}";
do
	filename="$(echo $val | rev | cut -d'/' -f1 | rev)"
	echo $filename
	Rscript $path_to_scripts_folder/convert_scores_to_ranks.R $val "$TMPWD/inputs_ranks/Ranks_$filename"

	bedtools intersect -a "$TMPWD/mergedPeaks.bed" -b "$TMPWD/inputs_ranks/Ranks_$filename" -wa -wb | awk -v var="$filename" '{print ($1"@"$2"@"$3"\t"$11"\t"var)}' >> "$TMPWD/all_tools"
done

rm "temp_ordered"

Rscript $path_to_scripts_folder/generate.R

awk '{print $2}' "$TMPWD/idr_Input" | awk -F"@" -vOFS='\t' '{gsub("\"", ""); print $1, $2, $3}' > "$TMPWD/locations"

awk -vOFS='\t' '{print ($3, $4, $5, $6, $7)}' "$TMPWD/idr_Input" > "$TMPWD/qms"
paste "$TMPWD/locations"  "$TMPWD/qms" | sort -k1,1V -k2,2n > "$TMPWD/Merged_Regions_Scores_1.tsv"
tail -n +2 "$TMPWD/Merged_Regions_Scores_1.tsv" > "$TMPWD/Merged_Regions_Scores.tsv"
rm "$TMPWD/Merged_Regions_Scores_1.tsv"

# (3) Replace NA (0) by ranks < 0 and run IDR

Rscript $path_to_scripts_folder/run_IDR.R

awk -vOFS='\t' '$1 != "V4"{print($1, $2, $3, $4, $5, $7, $8)}' IDR_out.tsv > t1


paste t1 temp/Merged_Regions_Scores.tsv | awk -vOFS='\t' '$6 <= 0.05{print($8, $9, $10, $11, $12, $13, $14, $15, $6, $7)}' > IDR_out.tsv

rm t1



