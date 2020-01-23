#!/bin/bash

name="run.sh"
path_to_config_file="config.tsv"
path_to_data_file="data.csv"

printHelp(){
	echo -e "" >&2
	echo -e "This script will dynamically create a snakemake pipeline and execute it. All required information is stored in config.tsv and data.csv." >&2
	echo -e "" >&2
	echo -e "" >&2
	echo -e "Usage: $name" >&2
	echo -e "" >&2
}

#Tools to run
bSmooth=FALSE
dss=FALSE
methylKit=FALSE
metilene=FALSE
rnBeads=FALSE

#Pros-Porcess Paramters
PP_minlen=5
PP_mindiff=0.05
PP_minCpG=5

# Start to create directory structure
mkdir -p 02_data
mkdir -p 02_data/prep
mkdir -p 03_results

cp config.tsv RB_config.tsv

# Read in relevant parameters from config file
while IFS=$'\t', read -r -a config
	do
		case "${config[0]}" in
			"BSmooth")	bSmooth=${config[1]}; 	t1="\n\t\t\"BSmooth_dummy.txt\",";;
			"DSS")		dss=${config[1]}; 	t2="\n\t\t\"DSS_dummy.txt\",";;
			"MethylKit")	methylKit=${config[1]}; t3="\n\t\t\"MethylKit_dummy.txt\",";;
			"Metilene")	metilene=${config[1]};	t4="\n\t\t\"Metilene_dummy.txt\",";;
			"RnBeads")	rnBeads=${config[1]};	t5="\n\t\t\"RnBeads_dummy.txt\"";;
			"PP_minlen")	PP_minlen=${config[1]};;
			"PP_mindiff")	PP_mindiff=${config[1]};;
			"PP_minCpG")	PP_minCpG=${config[1]};;
		esac
	done < $path_to_config_file

# Start to create Snakemake file
rm -f snakefile
touch snakefile
#echo -e "configfile: \"config.yaml\"\n\n" >> snakefile

# Add dummy rule
echo -e "rule finalize:\n\tinput:$t1$t2$t3$t4$t5\n\tshell:\n\t\t\"rm -rf BSmooth_dummy.txt\"\n\n" >> snakefile

# Add pre-processing rule
echo -e "rule run_pre_processing:\n\toutput:\n\t\t\"data_prep.csv\"\n\tinput:\n\t\t\"data.csv\",\n\t\t\"config.tsv\"\n\tshell:\n\t\t\"01_scripts/Pre-Processing/pre-process.sh\"\n\n" >> snakefile

echo -e "The following tools will be executed:\n"
if [ $bSmooth = "TRUE" ]; then
	echo " - BSmooth"
	echo -e "rule run_BSmooth:\n\tinput:\n\t\t\"data_prep.csv\"\n\toutput:\n\t\t\"02_data/BSmooth/group_A.txt\",\n\t\t\"02_data/BSmooth/group_B.txt\",\n\t\t\"03_results/BSmooth/BSmooth_DMRs_raw.tsv\",\n\t\t\"03_results/BSmooth/BSmooth_DMRs_std.tsv\",\n\t\t\"BSmooth_dummy.txt\"\n\tconda:\n\t\t\"envs/BSmooth.yml\"\n\tshell:\n\t\t\"01_scripts/BSmooth/run.sh -i -d -o\"\n\n" >> snakefile
fi
	mkdir -p 02_data/BSmooth
	mkdir -p 03_results/BSmooth
if [ $dss = "TRUE" ]; then
	echo " - DSS"
	echo -e "rule run_DSS:\n\tinput:\n\t\t\"data_prep.csv\"\n\toutput:\n\t\t\"02_data/DSS/group_A.txt\",\n\t\t\"02_data/DSS/group_B.txt\",\n\t\t\"03_results/DSS/DSS_DMRs_raw.tsv\",\n\t\t\"03_results/DSS/DSS_DMRs_std.tsv\",\n\t\t\"DSS_dummy.txt\"\n\tconda:\n\t\t\"envs/DSS.yml\"\n\tshell:\n\t\t\"01_scripts/DSS/run.sh -i -d -o\"\n\n" >> snakefile
	mkdir -p 02_data/DSS
	mkdir -p 03_results/DSS
fi
if [ $methylKit = "TRUE" ]; then
	echo " - MethylKit"
	echo -e "rule run_MethylKit:\n\tinput:\n\t\t\"data_prep.csv\"\n\toutput:\n\t\t\"03_results/MethylKit/MethylKit_DMRs_raw.tsv\",\n\t\t\"03_results/MethylKit/MethylKit_DMRs_std.tsv\",\n\t\t\"MethylKit_dummy.txt\"\n\tconda:\n\t\t\"envs/MethylKit.yml\"\n\tshell:\n\t\t\"01_scripts/MethylKit/run.sh -i -d -o\"\n\n" >> snakefile
	mkdir -p 02_data/MethylKit
	mkdir -p 03_results/MethylKit
fi
if [ $metilene = "TRUE" ]; then
	echo " - Metilene"
	echo -e "rule run_Metilene:\n\tinput:\n\t\t\"data_prep.csv\"\n\toutput:\n\t\t\"03_results/Metilene/Metilene_DMRs_raw.tsv\",\n\t\t\"03_results/Metilene/Metilene_DMRs_std.tsv\",\n\t\t\"Metilene_dummy.txt\"\n\tconda:\n\t\t\"envs/Metilene.yml\"\n\tshell:\n\t\t\"01_scripts/Metilene/run.sh -i -d -o\"\n\n" >> snakefile
	mkdir -p 02_data/Metilene
	mkdir -p 03_results/Metilene
fi
if [ $rnBeads = "TRUE" ]; then
	echo " - RnBeads"
	echo -e "rule run_RnBeads:\n\tinput:\n\t\t\"data.csv\"\n\toutput: \n\t\t\"02_data/RnBeads/sample_annotation.csv\",\n\t\t\"03_results/RnBeads/RnBeads_DMRs_raw.tsv\",\n\t\t\"03_results/RnBeads/RnBeads_DMRs_std.tsv\",\n\t\t\"RnBeads_dummy.txt\"\n\tconda:\n\t\t\"envs/RnBeads.yml\"\n\tshell:\n\t\t\"01_scripts/RnBeads/run.sh -i -d -o\"\n\n" >> snakefile
	mkdir -p 02_data/RnBeads
	mkdir -p 03_results/RnBeads
fi

snakemake --cores 12 --use-conda

wait

awk -vOFS='\t' -vmincpg="$PP_minCpG" -vminlen="$PP_minlen" -vmindiff="$PP_mindiff" '$4 >= mincpg && ($7 >= mindiff || $7 <= mindiff*(-1)) && ($3-$2+1) >= minlen {print}' '03_results/BSmooth/BSmooth_DMRs_std.tsv' > '03_results/BSmooth/BSmooth_DMRs_std_PP.tsv'

awk -vOFS='\t' -vmincpg="$PP_minCpG" -vminlen="$PP_minlen" -vmindiff="$PP_mindiff" '$4 >= mincpg && ($7 >= mindiff || $7 <= mindiff*(-1)) && ($3-$2+1) >= minlen {print}' '03_results/DSS/DSS_DMRs_std.tsv' > '03_results/DSS/DSS_DMRs_std_PP.tsv'

awk -vOFS='\t' -vmincpg="$PP_mindiff" -vminlen="$PP_minlen" '($7*100 >= mindiff || $7*100 <= mindiff*(-1)) && ($3-$2+1) >= minlen {print}' '03_results/MethylKit/MethylKit_DMRs_std.tsv' > '03_results/MethylKit/MethylKit_DMRs_std_PP.tsv'

awk -vOFS='\t' -vmincpg="$PP_minCpG" -vminlen="$PP_minlen" -vmindiff="$PP_mindiff" '$4 >= mincpg && ($7 >= mindiff || $7 <= mindiff*(-1)) && ($3-$2+1) >= minlen {print}' '03_results/Metilene/Metilene_DMRs_std.tsv' > '03_results/Metilene/Metilene_DMRs_std_PP.tsv'

awk -vOFS='\t' -vmincpg="$PP_minCpG" -vminlen="$PP_minlen" -vmindiff="$PP_mindiff" '$4 >= mincpg && ($7 >= mindiff || $7 <= mindiff*(-1)) && ($3-$2+1) >= minlen {print}' '03_results/RnBeads/RnBeads_DMRs_std.tsv' > '03_results/RnBeads/RnBeads_DMRs_std_PP.tsv'


<<<<<<< HEAD

=======
>>>>>>> df63897037eb1c759196430e3b9c977da0260d9d
