"""
Step 1 - Homozygosity mapping.

Real-world equivalent: PLINK's `--homozyg-group` (paper section "Homozygosity mapping").

The idea: if a disease is autosomal recessive, an affected dog inherited the SAME broken
copy of a gene from both parents. Near that gene, the affected dog's two chromosome
copies will look identical for a long stretch - a "run of homozygosity" (ROH). If several
affected dogs are related, they'll often share the SAME stretch, because they inherited
it from the same ancestor. Unaffected dogs usually won't share that exact stretch.

So the rule we're looking for, at each SNP marker, is:
    1. All AFFECTED dogs ("cases") are homozygous (e.g. "CC", not "CG") at this marker.
    2. All affected dogs have the SAME homozygous genotype (all "CC", not one "CC" and
       one "GG" - that would mean they got homozygous by chance for different alleles).
    3. The CONTROL dogs do NOT all share that same homozygous genotype - i.e. this
       stretch is specific to the disease, not just a common homozygous region in the
       breed (real quote from the paper: ROH must be "allelically shared by all three
       cases and either allelically different or absent in the controls").

A single marker matching that rule could be a coincidence. A long, unbroken STRETCH of
consecutive markers all matching it is what PLINK calls a "run of homozygosity", and
that's what we report.

Input:  data/snp_genotypes.csv   (toy SNP-array genotypes for 3 cases + 3 controls)
Output: results/step1_roh_regions.csv  (chrom, start_bp, end_bp, n_markers)
"""

import csv
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
INPUT_FILE = REPO_ROOT / "data" / "snp_genotypes.csv"
OUTPUT_FILE = REPO_ROOT / "results" / "step1_roh_regions.csv"

CASE_COLUMNS = ["case1", "case2", "case3"]
CONTROL_COLUMNS = ["control1", "control2", "control3"]

# A stretch shorter than this many markers is too short to count as a real ROH -
# it's the toy-scale stand-in for PLINK's --homozyg-snp minimum-marker-count setting.
MIN_MARKERS_PER_ROH = 4


def is_homozygous(genotype: str) -> bool:
    """'CC' -> True (both letters the same), 'CG' -> False."""
    return genotype[0] == genotype[1]


def marker_qualifies(row: dict) -> bool:
    """
    Does this one SNP marker look like it belongs to a case-specific ROH?
    Implements the three-part rule described in the module docstring.
    """
    case_genotypes = [row[col] for col in CASE_COLUMNS]
    control_genotypes = [row[col] for col in CONTROL_COLUMNS]

    all_cases_homozygous = all(is_homozygous(gt) for gt in case_genotypes)
    all_cases_match_each_other = len(set(case_genotypes)) == 1
    if not (all_cases_homozygous and all_cases_match_each_other):
        return False

    shared_case_genotype = case_genotypes[0]
    all_controls_match_cases = all(gt == shared_case_genotype for gt in control_genotypes)
    # We WANT this to be False: controls must differ from (or be absent from) the
    # cases' shared genotype, otherwise this is just a common homozygous spot in
    # the breed, not something specific to the affected dogs.
    return not all_controls_match_cases


def load_markers():
    with open(INPUT_FILE, newline="") as f:
        rows = list(csv.DictReader(f))
    # Group by chromosome, keep each chromosome's markers sorted by position -
    # a run of homozygosity only makes sense along one chromosome at a time.
    by_chrom = {}
    for row in rows:
        row["position_bp"] = int(row["position_bp"])
        by_chrom.setdefault(row["chrom"], []).append(row)
    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda r: r["position_bp"])
    return by_chrom


def find_roh_runs(markers_by_chrom):
    """
    Walk each chromosome's markers in position order and collect consecutive
    runs of qualifying markers. Returns a list of dicts describing each run
    long enough to count as a ROH.
    """
    roh_regions = []
    for chrom, markers in markers_by_chrom.items():
        current_run = []
        for marker in markers:
            if marker_qualifies(marker):
                current_run.append(marker)
            else:
                if current_run:
                    roh_regions.append(summarize_run(chrom, current_run))
                current_run = []
        if current_run:  # chromosome ended while still inside a run
            roh_regions.append(summarize_run(chrom, current_run))
    # Drop runs that are too short to be a believable shared ROH
    return [r for r in roh_regions if r["n_markers"] >= MIN_MARKERS_PER_ROH]


def summarize_run(chrom, run_markers):
    return {
        "chrom": chrom,
        "start_bp": run_markers[0]["position_bp"],
        "end_bp": run_markers[-1]["position_bp"],
        "n_markers": len(run_markers),
    }


def main():
    print("STEP 1: Homozygosity mapping")
    print(f"  reading {INPUT_FILE.relative_to(REPO_ROOT)}")

    markers_by_chrom = load_markers()
    all_runs_before_length_filter = []
    for chrom, markers in markers_by_chrom.items():
        n_qualifying = sum(1 for m in markers if marker_qualifies(m))
        print(f"  {chrom}: {n_qualifying}/{len(markers)} markers match the case-shared, "
              f"control-different rule")

    roh_regions = find_roh_runs(markers_by_chrom)

    print(f"  found {len(roh_regions)} case-specific ROH region(s) with >= "
          f"{MIN_MARKERS_PER_ROH} consecutive markers:")
    for region in roh_regions:
        span_mb = (region["end_bp"] - region["start_bp"]) / 1_000_000
        print(f"    {region['chrom']}:{region['start_bp']:,}-{region['end_bp']:,} "
              f"({region['n_markers']} markers, {span_mb:.1f} Mb)")

    OUTPUT_FILE.parent.mkdir(exist_ok=True)
    with open(OUTPUT_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["chrom", "start_bp", "end_bp", "n_markers"])
        writer.writeheader()
        writer.writerows(roh_regions)
    print(f"  wrote {OUTPUT_FILE.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
