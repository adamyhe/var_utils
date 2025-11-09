# ld_fetch.py
# Author: Adam Youlin He <adamyhe@gmail.com>

"""
A short script that uses the Ensembl REST API to fetch variants in LD with a provided
list of variants.
"""

import argparse
import concurrent.futures
import itertools
import sys

import requests
import tqdm


def rest_api_call(var, pop="CEU", threshold=0.8, metric="d_prime", wsize=200):
    """
    Function to make API call to Ensembl REST API

    Parameters
    ----------
    var : str
        rsID to calculate LD on
    pop : str, optional
        Population to calculate LD, by default "CEU"
    threshold : float, optional
        Threshold for LD, by default 0.8
    metric : str, optional
        Metric for LD, by default "d_prime"
    wsize : int, optional
        Window size for LD (kb), by default 200

    Returns
    -------
    list
        List of variant IDs
    """
    # generate request
    server = "https://rest.ensembl.org"
    ext = f"/ld/human/{var}/1000GENOMES:phase_3:{pop}?{metric}={threshold};window_size={wsize}"

    # send request
    r = requests.get(server + ext, headers={"Content-Type": "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    # parse results for variant IDs and return
    decoded = r.json()
    return [d["variation2"] for d in decoded]


def caller(
    in_var,
    pop="CEU",
    threshold=0.8,
    metric="d_prime",
    wsize=200,
    nthreads=8,
    verbose=False,
):
    """
    Wrapper function for rest_api_call that handles multithreading of API calls
    over a list of variants

    Parameters
    ----------
    in_var : list
        List of rsIDs to calculate LD on
    pop : str, optional
        Population to calculate LD, by default "CEU"
    threshold : float, optional
        Threshold for LD, by default 0.8
    metric : str, optional
        Metric for LD, by default "d_prime"
    wsize : int, optional
        Window size for LD (kb), by default 200
    nthreads : int, optional
        Number of threads, by default 8
    verbose : bool, optional
        print progress bar, by default False

    Returns
    -------
    list
        List of deduplicated variant IDs from API calls
    """
    # async multithreading of API calls
    if nthreads > 1:
        params = zip(
            [var for var in in_var],
            itertools.repeat(pop),
            itertools.repeat(threshold),
            itertools.repeat(metric),
            itertools.repeat(wsize),
        )
        with concurrent.futures.ThreadPoolExecutor(max_workers=nthreads) as executor:
            results = list(
                tqdm.tqdm(
                    executor.map(lambda p: rest_api_call(*p), params),
                    total=len(in_var),
                    disable=not verbose,
                    desc="Calculating LD",
                )
            )
    # or just run in for loop if nthreads == 1
    else:
        results = []
        for var in tqdm.tqdm(
            in_var, disable=not verbose, desc="Calculating LD", total=len(in_var)
        ):
            results.append(rest_api_call(var, pop, threshold, metric, wsize))

    # consolidate results and return
    flatten_dedup_results = list(set(itertools.chain.from_iterable(results)))
    return flatten_dedup_results


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        help="Text file with rsIDs to calculate LD on",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Output file to save LD variants",
        required=True,
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

    with open(args.input, "r") as f:
        in_var = list(set(f.read().splitlines()))

    print(f"Read {len(in_var)} unique rsIDs from {args.input}.")
    print(
        f"Filtering for variants in LD "
        f"({args.metric} > {args.threshold}, distance < {args.wsize}, population={args.pop})."
    )

    out_var = caller(
        in_var,
        args.pop,
        args.threshold,
        args.metric,
        args.wsize,
        args.nthreads,
        args.verbose,
    )

    print(f"Retrived {len(out_var)} unique rsIDs in LD.")
    with open(args.output, "w") as f:
        for var in out_var:
            f.write(var + "\n")


if __name__ == "__main__":
    main()
