# var_utils

A collection of small command line tools for handling genetic variants.

## Installation

```bash
pip install git+https://github.com/adamyhe/var_utils.git
```

## ld_fetch

This tool uses the Ensembl REST API to fetch variants in LD with a provided
list of variants based on 1K Genomes populations. It can launch multiple API calls at once to reduce runtime. NOTE: you should probably use a relatively small number still (e.g., <=8) to avoid getting in trouble with Ensembl sysadmins.

```bash
ld_fetch -h
ld_fetch --input rsids.txt --output rsids_ld.txt [--pop CEU --threshold 0.8 --metric d_prime/r2 --wsize 200 --nthreads 1 --verbose]
```

## lookup_dbsnp

This tool looks up coordinates, alleles, and allele frequency data from the UCSC dbSNP tables (e.g., https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp151.txt.gz) for a given list of rsIDs.

```bash
lookup_dbsnp -h
lookup_dbsnp --input rsids.txt --dbsnp snp151.txt.gz --output snp_data.tsv.gz [--tmp tmp_file.txt]
```

## cleanup_snps

This tool cleans up a tsv containing variants of interest (excludes all variants that contain non-ACGT characters nearby and will cause edge run-over).

```bash
cleanup_snps -h
cleanup_snps --input snps.tsv --fasta genome.fna --output snps_ACGT.tsv [--in_window 2114 --allowed_chars ACGT --verbose]
```

## match_allele_frequency

This tool samples a background variant set to match the allele frequency distribution of a target set of variants. Specifically, it bins the variants by MAF, then randomly samples matching counts of background variants.

```bash
match_allele_frequency -h
match_allele_frequency --df_target target_snps.tsv --df_background background_snps.tsv --output matched_snps.tsv [--frequency_bin_size 0.01 --random_seed 47 --verbose]
```
