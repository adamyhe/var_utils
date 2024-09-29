import random
import string
import subprocess
import sys

import pandas as pd

rsids = sys.argv[1]
dbsnp = sys.argv[2]
output = sys.argv[3]
rand_string = "".join(random.choices(string.ascii_uppercase + string.digits, k=8))
tmp = f"./lookup_dbsnp_alleles_tmp.{rand_string}.txt"

# Filter dbSNP text file using rsIDs and keep only SNPs on autosomes
subprocess.run(
    f"zgrep -wf {rsids} {dbsnp} | awk '$2 ~ /^chr[0-9]+$/' > {tmp}", shell=True
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

# Filter out SNPs with more than 2 alleles
out_df = out_df[pd.to_numeric(freqs[2]).isnull()]
# Clean nans
out_df.dropna(inplace=True)
# Write output
out_df.to_csv(output, sep="\t", index=False)
# Clean up tmp file
subprocess.run(f"rm {tmp}", shell=True)
