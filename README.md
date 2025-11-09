# var_utils

A collection of small command line tools for handling genetic variants.

## Installation

```bash
pip install git+https://github.com/adamyhe/var_utils.git
```

## ld_fetch

This tool uses the Ensembl REST API to fetch variants in LD with a provided
list of variants based on 1K Genomes populations. It can launch multiple API calls at once to reduce runtime. NOTE: you should probably use a relatively small number still (e.g., <0.8) to avoid getting in trouble with Ensembl sysadmins.

```bash
ld_fetch --input rsids.txt --output rsids_ld.txt [--pop CEU --threshold 0.8 --metric d_prime/r2 --wsize 200 --nthreads 1 --verbose]
```

## lookup_dbsnp

This tool looks up coordinates, alleles, and allele frequency data from the UCSC dbSNP tables (e.g., https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp151.txt.gz) for a given list of rsIDs.

```bash
lookup_dbsnp --input rsids.txt --dbsnp snp151.txt.gz --output snp_data.tsv.gz [--tmp tmp_file.txt]
```
