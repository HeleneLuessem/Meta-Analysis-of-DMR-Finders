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

# Start to create directory structure
mkdir -p 02_data
mkdir -p 02_data/prep
mkdir -p 03_results

# Start to create Snakemake file
rm -f snakefile
echo -e "configfile: \"config.yaml\"\n\n" >> snakefile

# Add Pre-Processing to Snakemake file


# Read in relevant parameters from config file
while IFS=$'\t', read -r -a config
	do
		case "${config[0]}" in
			"BSmooth")	bSmooth=${config[1]};;
			"DSS")		dss=${config[1]};;
			"MethylKit")	methylKit=${config[1]};;
			"Metilene")	metilene=${config[1]};;
			"RnBeads")	rnBeads=${config[1]};;
		esac
	done < $path_to_config_file

echo -e "The following tools will be executed:\n"
if [ $bSmooth = "TRUE" ]; then
	echo " - BSmooth"
	echo -e "rule run_BSmooth:\n\toutput: \n\t\t\"02_data/BSmooth/group_A.txt\",\n\t\t\"02_data/BSmooth/group_B.txt\",\n\t\t\"03_results/BSmooth/BSmooth_DMRs_raw.tsv\"\n\t\t\"03_results/BSmooth/BSmooth_DMRs_std.tsv\"\n\tconda:\n\t\t\"envs/BSmooth.yml\"\n\tshell:\n\t\t\"01_scripts/BSmooth/run.sh -i -d -o\"\n\n" >> snakefile
fi
	mkdir -p 02_data/BSmooth
	mkdir -p 03_results/BSmooth
if [ $dss = "TRUE" ]; then
	echo " - DSS"
	echo -e "rule run_DSS:\n\toutput: \n\t\t\"02_data/DSS/group_A.txt\",\n\t\t\"02_data/DSS/group_B.txt\",\n\t\t\"03_results/DSS/DSS_DMRs_raw.tsv\"\n\t\t\"03_results/DSS/DSS_DMRs_std.tsv\"\n\tconda:\n\t\t\"envs/DSS.yml\"\n\tshell:\n\t\t\"01_scripts/DSS/run.sh -i -d -o\"\n\n" >> snakefile
	mkdir -p 02_data/DSS
	mkdir -p 03_results/DSS
fi
if [ $methylKit = "TRUE" ]; then
	echo " - MethylKit"
	echo -e "rule run_MethylKit:\n\toutput: \n\t\t\"03_results/MethylKit/MethylKit_DMRs_raw.tsv\"\n\t\t\"03_results/MethylKit/MethylKit_DMRs_std.tsv\"\n\tconda:\n\t\t\"envs/MethylKit.yml\"\n\tshell:\n\t\t\"01_scripts/MethylKit/run.sh -i -d -o\"\n\n" >> snakefile
	mkdir -p 02_data/MethylKit
	mkdir -p 03_results/MethylKit
fi
if [ $metilene = "TRUE" ]; then
	echo " - Metilene"
	echo -e "rule run_Metilene:\n\toutput: \n\t\t\"03_results/Metilene/Metielne_DMRs_raw.tsv\"\n\t\t\"03_results/Metilene/Metilene_DMRs_std.tsv\"\n\tconda:\n\t\t\"envs/Metilene.yml\"\n\tshell:\n\t\t\"01_scripts/Metilene/run.sh -i -d -o\"\n\n" >> snakefile
	mkdir -p 02_data/Metilene
	mkdir -p 03_results/Metilene
fi
if [ $rnBeads = "TRUE" ]; then
	echo " - RnBeads"
	echo -e "rule run_RnBeads:\n\toutput: \n\t\t\"02_data/RnBeads/sample_annotation.csv\"\n\t\t\"03_results/RnBeads/RnBeads_DMRs_raw.tsv\"\n\t\t\"03_results/RnBeads/RnBeads_DMRs_std.tsv\"\n\tconda:\n\t\t\"envs/RnBeads.yml\"\n\tshell:\n\t\t\"01_scripts/RnBeads/run.sh -i -d -o\"\n\n" >> snakefile
	mkdir -p 02_data/RnBeads
	mkdir -p 03_results/RnBeads
fi

