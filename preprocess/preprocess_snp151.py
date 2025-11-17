# Preprocess dbSNP file from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp151.txt.gz
# For use with match_allele_frequency

import subprocess

import pandas as pd
import tqdm

# Filter dbSNP text file using rsIDs and keep only SNPs on autosomes
subprocess.run(
    "zcat snp151.txt.gz | awk '$2 ~ /^chr[0-9]+$/' | awk -F'\t' '$22' | pigz > snp151.wfreqs.txt.gz",
    shell=True,
)

# Parse filtered dbSNP file and extract relevant columns
df = pd.read_csv("snp151.wfreqs.txt.gz", sep="\t", header=None)
df = df[df[11] == "single"]
has_two_alleles = df[22].str.len() == 4
df = df[has_two_alleles]

freqs = df[24].str.split(",", expand=True)
out_df = pd.DataFrame(
    {
        "chrom": df[1],
        "pos": df[2],
        "rsid": df[4],
        "strand": df[6],
        "ref": df[7],
    }
)
snps = df[22].str.split(",", expand=True)

complement = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
}

alts = []
which = []
for i in tqdm.trange(snps.shape[0]):
    if out_df["strand"].iloc[i] == "+":
        if out_df["ref"].iloc[i] != snps[0].iloc[i]:
            alts.append(snps[0].iloc[i])
            which.append(0)
        else:
            alts.append(snps[1].iloc[i])
            which.append(1)
    elif out_df["strand"].iloc[i] == "-":
        if out_df["ref"].iloc[i] != complement[snps[0].iloc[i]]:
            alts.append(complement[snps[0].iloc[i]])
            which.append(0)
        else:
            alts.append(complement[snps[1].iloc[i]])
            which.append(1)
    else:
        raise ValueError(f"Invalid strand: {out_df['strand'].iloc[i]}")

out_df["alt"] = alts
out_df["alt_freq"] = pd.to_numeric(
    pd.Series(
        [
            (freqs[0].iloc[i] if which[i] == 0 else freqs[1].iloc[i])
            for i in tqdm.trange(snps.shape[0])
        ]
    )
)
out_df["ref_freq"] = pd.to_numeric(
    pd.Series(
        [
            (freqs[1].iloc[i] if which[i] == 0 else freqs[0].iloc[i])
            for i in tqdm.trange(snps.shape[0])
        ]
    )
)
is_biallelic = out_df["ref_freq"] + out_df["alt_freq"] == 1
out_df = out_df[is_biallelic]

out_df.to_csv("snp151.wfreqs.biallelic.tsv.gz", sep="\t", index=False)
