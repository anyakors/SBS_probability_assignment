# SBS_probability_assignment

SBS_probability_assignment is a library for the assignment of mutational signature probabilities at genomic end epigenomic features of interest. 

## Usage

bedtools should be installed in the environment where the scripts are executed.
Input folder should contain .vcf files; feature is a .bed file with genomic features of interest. 
Reference genome should be provided. Activities of COSMIC signatures should be pre-assigned.

```bash
	python calc_occupancy.py \
		--input data/snv_BRCA \
		--feature data/ctcf_sites_mcf7.bed \
		--chrom_lengths data/hg19.genome \
		--glob \
		--hg ../../hg19.fa \
		--activities data/valid_pcawg_published_activities_BRCA.txt \
		--out out_proba_CTCF_mcf7_global
```

## License

[MIT](https://choosealicense.com/licenses/mit/)
