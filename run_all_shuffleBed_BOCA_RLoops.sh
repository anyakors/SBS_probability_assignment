#!/bin/bash

stage=0
stop_stage=2

. ./parse_options.sh || exit 1;

if [ ${stage} -le 0 ] && [ ${stop_stage} -ge 0 ]; then

	python calc_occupancy_random_shuffleBed_pad.py \
		--input ../snv_named \
		--feature data/rloops_hg19_sorted.bed \
		--chrom_lengths data/hg19.genome \
		--perm 5 \
		--hg ../../hg19.fa \
		--hgver hg19 \
		--cosmicver 3.0 \
		--activities ../valid_pcawg_published_activities.txt \
		--out out_proba_BOCA_RLoops_bootstrap_trinuclcorr \
		--temp_folder temp3 \
		--pad 100

	rm temp3/*.bed

fi

if [ ${stage} -le 1 ] && [ ${stop_stage} -ge 1 ]; then

	python run_stats_single_bootstrap_shuffleBed.py \
		--folder out_proba_BOCA_RLoops_bootstrap_trinuclcorr \
		--activities ../valid_pcawg_published_activities.txt \
		--p_cutoff 0.05 \
		--N_perm 1000 \
		--N_sampled 10000 

fi

if [ ${stage} -le 2 ] && [ ${stop_stage} -ge 2 ]; then

	python run_stats_multisample_fromBootstrap.py \
		--folder out_proba_BOCA_RLoops_bootstrap_trinuclcorr \
		--activities ../valid_pcawg_published_activities.txt \
		--p_cutoff 0.05 \

fi