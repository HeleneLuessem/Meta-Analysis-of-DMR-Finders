#!/bin/bash

bedtools=bedtools

outPrefix="vennDiagram"
format="svg"
#declare -a allowedFormats=('svg' 'png' 'pdf')

printHelp() {
	echo -e "Usage: `basename $0` INPUTFILE_1 ... INPUTFILE_N"
	echo -e ""
	# echo -e " Mandatory:"
	# echo -e "  -i INPUTFILE\tcan be specified multiple times"
	echo -e " Optional:"
	echo -e "  -o OUTPREFIX\tprefix for result files (default: $outPrefix)"
	echo -e "  -t FILE\tname of output table"
	echo -e "  -d FILE\tname of diagram file"
	echo -e "  -f format\tformat of the picture: svg, png or pdf. (default: $format)"
}

mkdir -p tmp
TMPWD=tmp

while getopts ":ho:t:d:f:" opt
do
	case "$opt" in
		h) printHelp; exit 1 ;;
		o) outPrefix="$OPTARG" ;;
		t) tableName="$OPTARG" ;;
		d) diagramName="$OPTARG" ;;
		f) format="$OPTARG" ;;
#  i) i+=( "${OPTARG}" ) ;;
		*) printHelp; exit 1 ;;
	esac
done

shift $((OPTIND-1))

if (( $# < 2 ))
then
	echo "Too few input files given to calculate an intersection"
	echo ""
	printHelp
	exit 1
fi

cmd="$bedtools merge -i <(sort -k1,1 -k2,2n -k3,3n `for s in $*
do
	if [[ "$s" == *gz ]]
	then
		echo "<(zcat $s | cut -f 1-3)"
	else
		echo "<(cat $s | cut -f 1-3)"
	fi
done | tr '\n' ' '`) > $TMPWD/mergedPeaks.bed"

eval $cmd

echo "The merged set contains `cat $TMPWD/mergedPeaks.bed |wc -l` peaks" 1>&2

#awk -vOFS='\t' '{print ($1, $2, $3, -99999)}' mergedPeaks.bed > mergedPeaks2.bed

for s in $*
do
	echo "$s" 1>&2
	if [[ "$s" == *gz ]]
	then
		zcat $s
	else
		cat $s 
	fi | cut -f 1-3 | $bedtools intersect -c -a $TMPWD/mergedPeaks.bed -b - | awk -vOFS='\t' '$NF==0 {$NF=0} $NF>0 {$NF=$NF} {print}' > $TMPWD/mergedPeaks2.bed && mv $TMPWD/mergedPeaks2.bed $TMPWD/mergedPeaks.bed
done

mv $TMPWD/mergedPeaks.bed merged_count_DMRs_per_merged.bed
rm -rf $TMPWD
