python simulate_signatures.py \
	--signatures SBS3 SBS5 \
	--activities 0.8 0.2 \
	--hg19 ../../hg19.fa \
	--hgver hg19 \
	--N_samples 3 \
	--N_mut 5000 \
	--N_noise 1000 \
	--FCs 3 2 \
	--cosmicver 3.0 \
	--feature data/rloops_hg19_sorted.bed \
	--output out_simulated \