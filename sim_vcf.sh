#!/bin/bash

python make_vcf_simulation.py \
	--input data/snv_BRCA_more \
	--feature data/rloops_hg19_sorted.bed \
	--chrom_lengths data/hg19.genome \
	--perm 5 \
	--hg ../../hg19.fa \
	--hgver hg19 \
	--cosmicver 3.0 \
	--activities ../valid_pcawg_published_activities.txt \
	--out out_BRCA_SBS5_sim \
	--temp_folder temp3
