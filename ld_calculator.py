import argparse
import sys

import requests
import tqdm


def main(snps, pop="CEU", threshold=0.8, metric="d_prime", wsize=200, verbose=False):
    server = "https://rest.ensembl.org"
    return_list = []

    for snp in tqdm.tqdm(
        snps, desc="Fetching LD data", total=len(snps), verbose=verbose
    ):
        ext = f"/ld/human/{snp}/1000GENOMES:phase_3:{pop}?{metric}={threshold};window_size={wsize}"
        r = requests.get(server + ext, headers={"Content-Type": "application/json"})
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        decoded = r.json()
        return_list.extend([d["variation2"] for d in decoded])

    unique_snps = list(set(return_list))
    return unique_snps


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A short script that uses the REST API to retrieve variants in LD with a provided list.")
    parser.add_argument(
        "-s",
        "--snps",
        type=str,
        help="Text file with rsIDs to calculate LD",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Output file to save LD SNPs (if not provided, prints to stdout)",
        default=None,
    )
    parser.add_argument(
        "-p", "--pop", type=str, help="Population to calculate LD", default="CEU"
    )
    parser.add_argument(
        "-t", "--threshold", type=float, help="Threshold for LD", default=0.8
    )
    parser.add_argument(
        "-m", "--metric", type=str, help="Metric for LD", default="d_prime"
    )
    parser.add_argument(
        "-w", "--wsize", type=int, help="Window size for LD (kb)", default=200
    )
    parser.add_argument(
        "-v",
        "--verbose",
        help="print progress bar (disabled if -o provided)",
        action="store_true",
    )
    args = parser.parse_args()

    with open(args.snps, "r") as f:
        snps = f.read().splitlines()

    ld_snps = main(
        snps,
        args.pop,
        args.threshold,
        args.metric,
        args.wsize,
        args.verbose and args.output is None,
    )

    if args.output:
        with open(args.output, "w") as f:
            for snp in ld_snps:
                f.write(snp + "\n")

    else:
        for snp in ld_snps:
            print(snp)
