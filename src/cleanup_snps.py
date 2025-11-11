# cleanup_snps.py
# Adam He <adamyhe@gmail.com>

"""
Tool to clean up a tsv containing SNPs of interest (such as one produced by
lookup_dbsnp). Excludes all SNPs that contain non-ACGT characters nearby
(default = 2114 bases).
"""

import argparse

import pandas as pd
import pyfastx
import tqdm


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input SNPs tsv. Required columns: chrom, pos",
    )
    parser.add_argument("-f", "--fasta", type=str, required=True)
    parser.add_argument("-o", "--output", type=str, required=True)
    parser.add_argument(
        "-w", "--in_window", type=int, default=2114, help="Window size to check"
    )
    parser.add_argument(
        "-a", "--allowed_chars", type=str, default="ACGT", help="Allowed characters"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Whether to print progress bar"
    )
    args = parser.parse_args()

    # Read in sequences around SNPs
    snps = pd.read_csv(args.snps_tsv, sep="\t")
    coords = snps[["chrom", "pos"]].copy()
    coords["end"] = coords["pos"] + 1
    fa = pyfastx.Fasta(args.fa_fname)
    wholesome = []
    for row in tqdm.tqdm(
        snps.itertuples(), total=snps.shape[0], disable=not args.verbose
    ):
        chrom = row[1]
        center = row[2]
        start = max(0, center - args.in_window // 2)
        end = center + args.in_window // 2 - 1  # pyfastx includes the end
        seq = fa.fetch(chrom, (start, end)).upper()
        is_wholesome = (
            all([c in args.allowed_chars for c in seq]) and len(seq) == args.in_window
        )
        wholesome.append(is_wholesome)

    # Dump cleaned SNPs to tsv
    print(
        "Filtered out {} due to disallowed characters or length != {}".format(
            snps.shape[0] - sum(wholesome), args.in_window
        )
    )
    snps = snps[wholesome]
    snps.reset_index(drop=True, inplace=True)
    snps.to_csv(args.out_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
