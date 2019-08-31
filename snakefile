configfile: "config.yaml"


rule dummy:
	input:
		"BSmooth_dummy.txt",
		"DSS_dummy.txt",
		"MethylKit_dummy.txt",
		"Metilene_dummy.txt",
		"RnBeads_dummy.txt"
	shell:
		"rm -rf BSmooth_dummy.txt"


rule run_pre_processing:
	output:
		"data_prep.csv"
	input:
		"data.csv",
		"config.tsv"
	shell:
		"01_scripts/Pre-Processing/pre-process.sh"


rule run_BSmooth:
	input:
		"data_prep.csv"
	output:
		"02_data/BSmooth/group_A.txt",
		"02_data/BSmooth/group_B.txt",
		"03_results/BSmooth/BSmooth_DMRs_raw.tsv",
		"03_results/BSmooth/BSmooth_DMRs_std.tsv",
		"BSmooth_dummy.txt"
	conda:
		"envs/BSmooth.yml"
	shell:
		"01_scripts/BSmooth/run.sh -i -d -o"


rule run_DSS:
	input:
		"data_prep.csv"
	output:
		"02_data/DSS/group_A.txt",
		"02_data/DSS/group_B.txt",
		"03_results/DSS/DSS_DMRs_raw.tsv",
		"03_results/DSS/DSS_DMRs_std.tsv",
		"DSS_dummy.txt"
	conda:
		"envs/DSS.yml"
	shell:
		"01_scripts/DSS/run.sh -i -d -o"


rule run_MethylKit:
	input:
		"data_prep.csv"
	output:
		"03_results/MethylKit/MethylKit_DMRs_raw.tsv",
		"03_results/MethylKit/MethylKit_DMRs_std.tsv",
		"MethylKit_dummy.txt"
	conda:
		"envs/MethylKit.yml"
	shell:
		"01_scripts/MethylKit/run.sh -i -d -o"


rule run_Metilene:
	input:
		"data_prep.csv"
	output:
		"03_results/Metilene/Metielne_DMRs_raw.tsv",
		"03_results/Metilene/Metilene_DMRs_std.tsv",
		"Metilene_dummy.txt"
	conda:
		"envs/Metilene.yml"
	shell:
		"01_scripts/Metilene/run.sh -i -d -o"


rule run_RnBeads:
	input:
		"data.csv"
	output: 
		"02_data/RnBeads/sample_annotation.csv",
		"03_results/RnBeads/RnBeads_DMRs_raw.tsv",
		"03_results/RnBeads/RnBeads_DMRs_std.tsv",
		"RnBeads_dummy.txt"
	conda:
		"envs/RnBeads.yml"
	shell:
		"01_scripts/RnBeads/run.sh -i -d -o"


