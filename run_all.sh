#!/bin/bash

stage=0
stop_stage=0

. ./parse_options.sh || exit 1;

if [ ${stage} -le 0 ] && [ ${stop_stage} -ge 0 ]; then

	python calc_occupancy_new.py \
		--input data/snv_sbs8 \
		--feature data/GSM5669226_H9ESC_LADs_hg19.bed \
		--chrom_lengths data/hg19.genome \
		--glob \
		--hg ../../hg19.fa \
		--hgver hg19 \
		--cosmicver 3.0 \
		--activities data/valid_pcawg_published_activities_sbs8.txt \
		--out out_proba_LAD_H9ESC_sbs8_

	#rm temp/*.bed

fi

if [ ${stage} -le 1 ] && [ ${stop_stage} -ge 1 ]; then

	python run_stats_single_bootstrap.py \
		--folder out_proba_LAD_H9ESC_sbs8_ \
		--activities data/valid_pcawg_published_activities_sbs8.txt \
		--p_cutoff 0.05 \
		--N_perm 1000 \
		--N_sampled 10000 

fi

if [ ${stage} -le 2 ] && [ ${stop_stage} -ge 2 ]; then

	python run_stats_multisample.py \
		--folder out_proba_LAD_H9ESC_sbs8_ \
		--activities data/valid_pcawg_published_activities_sbs8.txt \
		--p_cutoff 0.05 \

fi