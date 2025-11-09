# lookup_dbsnp.py
# Author: Adam Youlin He <adamyhe@gmail.com>

"""
Script to fetch coordinates and allele information from a UCSC dbSNP file for
a given list of rsIDs. For hg38, download from
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp151.txt.gz
"""

import argparse
import os
import random
import string
import subprocess

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        help="Text file with rsIDs to look up.",
        required=True,
    )
    parser.add_argument(
        "-d",
        "--dbsnp",
        type=str,
        help="Path to gzipped UCSC dbSNP file.",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Output file path.",
        required=True,
    )
    parser.add_argument(
        "-t",
        "--tmp",
        type=str,
        default=None,
        help="Temporary file path to store intermediate results.",
    )
    args = parser.parse_args()

    # Determine temporary file path
    rand_string = "".join(random.choices(string.ascii_uppercase + string.digits, k=8))
    if args.tmp is None:
        tmp = f"./lookup_dbsnp_alleles_tmp.{rand_string}.txt"
    elif os.path.isfile(args.tmp):
        tmp = args.tmp
    elif os.path.isdir(args.tmp):
        tmp = os.path.join(args.tmp, f"lookup_dbsnp_alleles_tmp.{rand_string}.txt")
    else:
        raise ValueError(f"Invalid tmp argument: {args.tmp}")

    # Filter dbSNP text file using rsIDs and keep only SNPs on autosomes
    subprocess.run(
        f"zgrep -wf {args.input} {args.dbsnp} | awk '$2 ~ /^chr[0-9]+$/' > {tmp}",
        shell=True,
    )

    # Parse filtered dbSNP file and extract relevant columns
    df = pd.read_csv(tmp, sep="\t", header=None)
    df = df[df[11] == "single"]
    freqs = df[24].str.split(",", expand=True)
    out_df = pd.DataFrame(
        {
            "chrom": df[1],
            "pos": df[2],
            "rsid": df[4],
            "ref": df[8],
            "ref_freq": pd.to_numeric(freqs[0]),
            "alt_freq": pd.to_numeric(freqs[1]),
        }
    )
    snps = df[22].str.split(",", expand=True)
    out_df["alt"] = [
        (snps[0].iloc[i] if snps[0].iloc[i] != out_df.ref.iloc[i] else snps[1].iloc[i])
        for i in range(snps.shape[0])
    ]

    # Clean nans
    out_df.dropna(inplace=True)
    # Write output
    out_df.to_csv(args.output, sep="\t", index=False)
    # Clean up tmp file
    subprocess.run(f"rm {tmp}", shell=True)


if __name__ == "__main__":
    main()
