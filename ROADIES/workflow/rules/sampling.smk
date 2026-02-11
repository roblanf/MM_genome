import glob
from collections import OrderedDict
import random,os
from pathlib import Path
import subprocess

num_species = len(os.listdir(config["GENOMES"]))
num_genomes = len(SAMPLES)

# Determine how many loci each genome contributes
loci_counts = OrderedDict([(key, 0) for key in SAMPLES])

focal_genome = config.get("FOCAL_GENOME", "")
if focal_genome:
	# Strip extension to get sample name (handle .fa.gz and .fa)
	focal_name = focal_genome
	if focal_name.endswith(".fa.gz"):
		focal_name = focal_name[:-6]
	elif focal_name.endswith(".fa"):
		focal_name = focal_name[:-3]

	if focal_name not in loci_counts:
		raise ValueError(
			f"FOCAL_GENOME '{focal_genome}' not found in input genomes. "
			f"Available genomes: {list(loci_counts.keys())}"
		)
	loci_counts[focal_name] = num
else:
	for i in range(num):
		index = random.randint(0, num_genomes - 1)
		loci_counts[SAMPLES[index]] = loci_counts[SAMPLES[index]] + 1

# Calculate gene ID ranges (1-indexed, inclusive [start, end])
# Genomes with 0 loci get start=0, end=0 as a sentinel for "skip"
gene_starts = OrderedDict([(key, 0) for key in SAMPLES])
gene_ends = OrderedDict([(key, 0) for key in SAMPLES])
next_id = 1
for sample in SAMPLES:
	count = loci_counts[sample]
	if count > 0:
		gene_starts[sample] = next_id
		gene_ends[sample] = next_id + count - 1
		next_id += count

def get_index_s(wildcards):
	return gene_starts[wildcards]

def get_index_e(wildcards):
	return gene_ends[wildcards]

def should_skip(wildcards):
	return loci_counts[wildcards] == 0

rule sequence_select:
	input:
		config["GENOMES"] + "/{sample}." + ("fa.gz" if EXTENSION[0]=="gz" else "fa")
	params:
		LENGTH=config["LENGTH"],
		KFAC=lambda wildcards: get_index_s(wildcards.sample),
		KFAC_e=lambda wildcards: get_index_e(wildcards.sample),
		THRES=config["UPPER_CASE"],
		SKIP=lambda wildcards: should_skip(wildcards.sample)
	benchmark:
		config["OUT_DIR"]+"/benchmarks/{sample}.sample.txt"
	threads: lambda wildcards: int(config['num_threads'])
	output:
        	config["OUT_DIR"]+"/samples/{sample}_temp.fa"
	shell:
			'''
			if [ "{params.SKIP}" = "True" ]; then
				echo "Skipping {input} (no loci assigned)"
				touch {output}
			else
				echo "We are starting to sample {input}"
				echo "./workflow/scripts/sampling/build/sampling -i {input} -o {output} -l {params.LENGTH} -s {params.KFAC} -e {params.KFAC_e} -t {params.THRES}"
				time ./workflow/scripts/sampling/build/sampling -i {input} -o {output} -l {params.LENGTH} -s {params.KFAC} -e {params.KFAC_e} -t {params.THRES}
			fi
			'''

rule sequence_merge:
	input:
		expand(config["OUT_DIR"]+"/samples/{sample}_temp.fa", sample=SAMPLES),
	params:
		gene_dir = config["OUT_DIR"]+"/samples",
		plotdir = config["OUT_DIR"]+"/plots",
		statdir = config["OUT_DIR"]+"/statistics"
	output:
        	config["OUT_DIR"]+"/samples/out.fa",
			report(config["OUT_DIR"]+"/plots/sampling.png",caption="../report/sampling.rst",category='Sampling Report')
	shell:
		"python3 workflow/scripts/sequence_merge.py {params.gene_dir} {output} {params.plotdir} {params.statdir}"
