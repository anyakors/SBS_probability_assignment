#!/bin/bash

stage=0
stop_stage=2

. ./parse_options.sh || exit 1;

if [ ${stage} -le 0 ] && [ ${stop_stage} -ge 0 ]; then

	python calc_occupancy_random_shuffleBed_cells.py \
		--input data/snv_cell_lines \
		--feature data/methylation_1kb_upstream_TSS \
		--chrom_lengths data/hg19.genome \
		--perm 5 \
		--hg ../../hg19.fa \
		--hgver hg19 \
		--cosmicver 3.0 \
		--activities data/cancer_cell_lines_sigs.txt \
		--out out_proba_cell_lines_meth_trinucl_corr \
		--temp_folder temp2 \
		--correct_trinucl

	rm temp2/*.bed

fi

if [ ${stage} -le 1 ] && [ ${stop_stage} -ge 1 ]; then

	python run_stats_single_bootstrap_shuffleBed.py \
		--folder out_proba_cell_lines_meth_trinucl_corr \
		--activities data/cancer_cell_lines_sigs.txt \
		--p_cutoff 0.05 \
		--N_perm 1000 \
		--N_sampled 10000 

fi

if [ ${stage} -le 2 ] && [ ${stop_stage} -ge 2 ]; then

	python run_stats_multisample_fromBootstrap.py \
		--folder out_proba_cell_lines_meth_trinucl_corr \
		--activities data/cancer_cell_lines_sigs.txt \
		--p_cutoff 0.05 \

fi