#!/bin/bash

stage=1
stop_stage=2

. ./parse_options.sh || exit 1;

if [ ${stage} -le 0 ] && [ ${stop_stage} -ge 0 ]; then

	python calc_occupancy_random_shuffleBed.py \
		--input data/snv_LIRI \
		--feature data/cpg_islands_hg19.bed \
		--chrom_lengths data/hg19.genome \
		--perm 5 \
		--hg ../../hg19.fa \
		--hgver hg19 \
		--cosmicver 3.0 \
		--activities data/valid_pcawg_published_activities_LIRI.txt \
		--out out_proba_LIRI_CpG_bootstrap \
		--temp_folder temp3

	rm temp3/*.bed

fi

if [ ${stage} -le 1 ] && [ ${stop_stage} -ge 1 ]; then

	python run_stats_single_bootstrap_shuffleBed.py \
		--folder out_proba_LIRI_CpG_bootstrap \
		--activities data/valid_pcawg_published_activities_LIRI.txt \
		--p_cutoff 0.05 \
		--N_perm 1000 \
		--N_sampled 10000 

fi

if [ ${stage} -le 2 ] && [ ${stop_stage} -ge 2 ]; then

	python run_stats_multisample_fromBootstrap.py \
		--folder out_proba_LIRI_CpG_bootstrap \
		--activities data/valid_pcawg_published_activities_LIRI.txt \
		--p_cutoff 0.05 \

fi