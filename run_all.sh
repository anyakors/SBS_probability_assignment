#!/bin/bash

stage=0
stop_stage=2

. ./parse_options.sh || exit 1;

if [ ${stage} -le 0 ] && [ ${stop_stage} -ge 0 ]; then

	python calc_occupancy.py \
		--input data/snv_BRCA \
		--feature data/ctcf_sites_mcf7.bed \
		--chrom_lengths data/hg19.genome \
		--glob \
		--hg ../../hg19.fa \
		--activities data/valid_pcawg_published_activities_BRCA.txt \
		--out out_proba_CTCF_mcf7_global

	rm temp/*.bed

fi

if [ ${stage} -le 1 ] && [ ${stop_stage} -ge 1 ]; then

	python run_stats_single_bootstrap.py \
		--folder out_proba_CTCF_mcf7_global \
		--activities data/valid_pcawg_published_activities_BRCA.txt \
		--p_cutoff 0.05 \
		--N_perm 1000 \
		--N_sampled 10000 

fi

if [ ${stage} -le 2 ] && [ ${stop_stage} -ge 2 ]; then

	python run_stats_multisample.py \
		--folder out_proba_CTCF_mcf7_global \
		--activities data/valid_pcawg_published_activities_BRCA.txt \
		--p_cutoff 0.05 \

fi