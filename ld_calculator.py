import argparse
import concurrent.futures
import itertools
import sys

import requests
import tqdm


def rest_api_call(snp, pop="CEU", threshold=0.8, metric="d_prime", wsize=200):
    server = "https://rest.ensembl.org"
    ext = f"/ld/human/{snp}/1000GENOMES:phase_3:{pop}?{metric}={threshold};window_size={wsize}"

    r = requests.get(server + ext, headers={"Content-Type": "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()
    return [d["variation2"] for d in decoded]


def main(
    snps,
    pop="CEU",
    threshold=0.8,
    metric="d_prime",
    wsize=200,
    nthreads=8,
    verbose=False,
):
    if nthreads > 1:
        params = zip(
            [snp for snp in snps],
            itertools.repeat(pop),
            itertools.repeat(threshold),
            itertools.repeat(metric),
            itertools.repeat(wsize),
        )
        with concurrent.futures.ThreadPoolExecutor(max_workers=nthreads) as executor:
            results = list(
                tqdm.tqdm(
                    executor.map(lambda p: rest_api_call(*p), params),
                    total=len(snps),
                    disable=not verbose,
                    desc="Calculating LD",
                )
            )
    else:
        results = []
        for snp in tqdm.tqdm(
            snps, disable=not verbose, desc="Calculating LD", total=len(snps)
        ):
            results.append(rest_api_call(snp, pop, threshold, metric, wsize))

    flatten_dedup_results = list(set(itertools.chain.from_iterable(results)))
    return flatten_dedup_results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A short script that uses the REST API to retrieve variants in LD with a provided list."
    )
    parser.add_argument(
        "-s",
        "--snps",
        type=str,
        help="Text file with rsIDs to calculate LD",
        required=True,
    )
    parser.add_argument(
        "-o", "--output", type=str, help="Output file to save LD SNPs", required=True
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
        "-n", "--nthreads", type=int, help="Number of threads", default=1
    )
    parser.add_argument(
        "-v", "--verbose", help="print progress bar", action="store_true"
    )
    args = parser.parse_args()

    with open(args.snps, "r") as f:
        snps = list(set(f.read().splitlines()))

    print(f"Read {len(snps)} unique rsIDs from {args.snps}.")

    ld_snps = main(
        snps,
        args.pop,
        args.threshold,
        args.metric,
        args.wsize,
        args.nthreads,
        args.verbose,
    )
    with open(args.output, "w") as f:
        for snp in ld_snps:
            f.write(snp + "\n")
