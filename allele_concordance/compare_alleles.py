#!/usr/bin/env python
# can be run as python compare.py <illumina.txt> <nanopore.txt> <detail_output> <genotype_summary>
# or as an executable ./compare.py <illumina.txt> <nanopore.txt> <detail_output> <genotype_summary>

import csv
from collections import defaultdict

import sys
if len(sys.argv) != 5:
	print("Usage: python compare.py <illumina_file> <nanopore_file> <detail_output> <genotype_summary>")
	sys.exit(1)

ILLUMINA_FILE = sys.argv[1]
NANOPORE_FILE = sys.argv[2]
OUT_DETAIL = sys.argv[3]
OUT_SUMMARY = sys.argv[4]

def parse_value(val):
    """Return (gt, dp), normalizing missing values."""
    if val is None or val == "" or val == ".":
        return "./.", "NA"
    if ":" not in val:
        return val, "NA"
    gt, dp = val.split(":", 1)
    return gt, dp


def is_het(gt):
    return gt in ("0/1", "1/0")


def classify(gtI, gtN):
    """Return classification type from the AWK logic."""
    # Missingness logic
    if gtI == "./." and gtN == "./.":
        return "missing_both"
    if gtI == "./.":
        return "missing_illumina"
    if gtN == "./.":
        return "missing_nanopore"

    # Exact match (including "0/1" ↔ "1/0")
    if gtI == gtN or (gtI, gtN) in (("0/1", "1/0"), ("1/0", "0/1")):
        return "match"

    # Partial match logic (het ↔ homo transitions)
    partial_pairs = {
        ("0/1", "0/0"), ("0/1", "1/1"),
        ("1/0", "0/0"), ("1/0", "1/1")
    }

    if (gtI, gtN) in partial_pairs or (gtN, gtI) in partial_pairs:
        return "partial_match"

    # Heterozygous ↔ homozygous transitions for 1/1 and 0/0
    if gtI == "1/1" and gtN in ("0/1", "1/0"):
        return "partial_match"
    if gtI == "0/0" and gtN in ("0/1", "1/0"):
        return "partial_match"

    return "no_match"


# ------------------------------------------------------------------------------
# Read Illumina file
# ------------------------------------------------------------------------------
illum = defaultdict(dict)
seen_illum = set()
all_samples = set()

with open(ILLUMINA_FILE) as f:
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)
    illum_samples = header[2:]
    all_samples.update(illum_samples)

    for row in reader:
        chrom, pos = row[0], row[1]
        key = f"{chrom}:{pos}"
        seen_illum.add(key)

        for s, val in zip(illum_samples, row[2:]):
            illum[key][s] = val


# ------------------------------------------------------------------------------
# Read Nanopore file + compare
# ------------------------------------------------------------------------------
count = defaultdict(lambda: defaultdict(int))
seen = set()

# open detail output now
detail = open(OUT_DETAIL, "w")
detail.write("CHROM\tPOS\tSAMPLE\tTYPE\tDEPTH_illumina\tDEPTH_nanopore\tHET_illumina\tHET_nanopore\n")

with open(NANOPORE_FILE) as f:
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)
    nano_samples = header[2:]
    all_samples.update(nano_samples)

    for row in reader:
        chrom, pos = row[0], row[1]
        key = f"{chrom}:{pos}"
        seen.add(key)

        # Build nanopore sample dict for this row
        nano_vals = dict(zip(nano_samples, row[2:]))

        for s in all_samples:
            # Illumina value (may not exist)
            illum_val = illum[key].get(s, "./.:NA")
            gtI, dpI = parse_value(illum_val)

            # Nanopore value (may not exist)
            nano_val = nano_vals.get(s, "./.:NA")
            gtN, dpN = parse_value(nano_val)

            hetI = "yes" if is_het(gtI) else "no"
            hetN = "yes" if is_het(gtN) else "no"

            type_ = classify(gtI, gtN)
            count[s][type_] += 1

            detail.write(
                f"{chrom}\t{pos}\t{s}\t{type_}\t{dpI}\t{dpN}\t{hetI}\t{hetN}\n"
            )

detail.close()


# ------------------------------------------------------------------------------
# Handle Illumina-only SNPs
# ------------------------------------------------------------------------------
detail = open(OUT_DETAIL, "a")

for key in seen_illum:
    if key in seen:
        continue

    chrom, pos = key.split(":")

    for s in all_samples:
        illum_val = illum[key].get(s, "./.:NA")
        gtI, dpI = parse_value(illum_val)
        gtN, dpN = "./.", "NA"

        hetI = "yes" if is_het(gtI) else "no"
        hetN = "no"

        if gtI == "./.":
            type_ = "missing_both"
        else:
            type_ = "missing_nanopore"

        count[s][type_] += 1

        detail.write(
            f"{chrom}\t{pos}\t{s}\t{type_}\t{dpI}\t{dpN}\t{hetI}\t{hetN}\n"
        )

detail.close()


# ------------------------------------------------------------------------------
# Write summary file
# ------------------------------------------------------------------------------
with open(OUT_SUMMARY, "w") as out:
    out.write(
        "SAMPLE\tmatch\tpartial_match\tno_match\tmissing_illumina\tmissing_nanopore\tmissing_both\n"
    )
    for s in sorted(all_samples):
        out.write(
            f"{s}\t"
            f"{count[s]['match']}\t"
            f"{count[s]['partial_match']}\t"
            f"{count[s]['no_match']}\t"
            f"{count[s]['missing_illumina']}\t"
            f"{count[s]['missing_nanopore']}\t"
            f"{count[s]['missing_both']}\n"
        )
