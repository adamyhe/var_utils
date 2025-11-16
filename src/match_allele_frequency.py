# match_allele_frequency.py
# Adam He <adamyhe@gmail.com>

"""
Tool to sample a background variant set to match the allele frequency distribution
of a target set of SNPs.
"""

import argparse
import logging

import numpy as np
import pandas as pd
import tqdm


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-t",
        "--df_target",
        type=str,
        help="Tab separated file with target SNPs to match (must have column 'ref_freq').",
        required=True,
    )
    parser.add_argument(
        "-b",
        "--df_background",
        type=str,
        help="Tab separated file with background SNPs that will be sampled (must have column 'ref_freq').",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Output file path to save matched SNPs.",
        required=True,
    )
    parser.add_argument(
        "-f",
        "--frequency_bin_size",
        type=float,
        default=0.01,
        help="Size of bins to group SNP frequencies by.",
    )
    parser.add_argument(
        "-s",
        "--random_seed",
        type=int,
        default=None,
        help="Seed for random number generator.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Whether to print progress bar.",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    if args.frequency_bin_size >= 1 or args.frequency_bin_size <= 0:
        raise ValueError("Bin size must be between 0 and 1.")
    if args.random_seed is not None:
        np.random.seed(args.random_seed)

    # Load data
    df_target = pd.read_csv(args.df_target, sep="\t")
    if args.verbose:
        logging.info(f"Loaded {len(df_target)} target SNPs")
    df_bg = pd.read_csv(args.df_background, sep="\t")
    if args.verbose:
        logging.info(f"Loaded {len(df_bg)} background SNPs")

    # Bin SNPs by frequency
    bins = np.arange(0, 1 + args.frequency_bin_size / 10, args.frequency_bin_size)
    df_target["bin"] = pd.cut(df_target["ref_freq"], bins=bins, include_lowest=True)
    df_bg["bin"] = pd.cut(df_bg["alt_freq"], bins=bins, include_lowest=True)

    # Sample SNPs by frequency
    matched = []
    for b in tqdm.tqdm(
        df_target["bin"].unique(),
        disable=not args.verbose,
        desc=f"Sampling in bins of {args.frequency_bin_size}",
    ):
        n = sum(df_target["bin"] == b)
        pool = df_bg[df_bg["bin"] == b]
        sampled = pool.sample(n, replace=False)
        matched.append(sampled)

    # Concatenate sampled SNPs
    df_matched = pd.concat(matched)
    df_matched.drop(columns=["bin"], inplace=True)
    df_matched.sort_values(["chrom", "pos"], inplace=True)

    # Save matched SNPs
    df_matched.to_csv(args.output, sep="\t", index=False)
    if args.verbose:
        logging.info(f"Saved {df_matched.shape[0]} matched SNPs to {args.output}")


if __name__ == "__main__":
    main()
