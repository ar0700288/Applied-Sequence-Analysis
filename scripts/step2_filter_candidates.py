"""
Step 2 - Recessive-model candidate filtering.

Real-world equivalent: the webGQT filtering described in the paper's "Variant analysis"
section (control-panel filtering) combined with the "case-specific ROH" intersection
from the Results section.

Whole-genome/exome sequencing of the two affected dogs turns up many DNA variants versus
the reference genome - most of them are normal breed variation, not the disease cause.
We narrow that list down using four independent checks. A variant only survives if it
passes ALL four:

  1. RECESSIVE INHERITANCE  - both sequenced affected dogs must be homozygous for the
     variant (hom_alt/hom_alt). A dominant or randomly-inherited variant wouldn't show
     this clean pattern.
  2. RARE IN CONTROLS        - a real disease-causing recessive variant should be rare.
     The paper allowed at most 2 heterozygous carriers, and 0 homozygous-affected
     individuals, in its control-genome panel. A variant that's common in unaffected
     dogs can't be the cause of a rare disease.
  3. INSIDE A CASE-SPECIFIC ROH - the variant's position must fall inside one of the
     shared homozygous regions found in Step 1. This is what actually ties the sequencing
     result back to the homozygosity-mapping result.
  4. PROTEIN-CHANGING (missense) - only variants predicted to change the protein
     sequence can plausibly break the protein's function, so only "missense" survives
     out of {missense, intronic, UTR, synonymous, ...}.

Input:  data/candidate_variants.csv     (variants found in the two sequenced case dogs)
        data/control_panel.csv          (control-cohort genotypes at those sites)
        results/step1_roh_regions.csv   (output of step 1)
Output: results/step2_filtered_candidates.csv
"""

import csv
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
CANDIDATES_FILE = REPO_ROOT / "data" / "candidate_variants.csv"
CONTROL_PANEL_FILE = REPO_ROOT / "data" / "control_panel.csv"
ROH_FILE = REPO_ROOT / "results" / "step1_roh_regions.csv"
OUTPUT_FILE = REPO_ROOT / "results" / "step2_filtered_candidates.csv"

MAX_HET_ALLOWED_IN_CONTROLS = 2
MAX_HOM_ALT_ALLOWED_IN_CONTROLS = 0


def load_roh_regions():
    with open(ROH_FILE, newline="") as f:
        regions = list(csv.DictReader(f))
    for region in regions:
        region["start_bp"] = int(region["start_bp"])
        region["end_bp"] = int(region["end_bp"])
    return regions


def load_control_genotypes():
    """variant_id -> list of control genotype calls, e.g. ['het', 'hom_ref', ...]"""
    with open(CONTROL_PANEL_FILE, newline="") as f:
        rows = list(csv.DictReader(f))
    control_columns = [c for c in rows[0].keys() if c != "variant_id"]
    return {row["variant_id"]: [row[c] for c in control_columns] for row in rows}


def passes_case_genotype_check(variant: dict) -> bool:
    """Both sequenced affected dogs must be homozygous for the variant."""
    return variant["seq_case_wes_gt"] == "hom_alt" and variant["seq_case_wgs_gt"] == "hom_alt"


def passes_control_frequency_check(control_genotypes: list) -> bool:
    """Rare-in-controls check: <=2 het carriers, 0 homozygous-affected controls."""
    n_het = control_genotypes.count("het")
    n_hom_alt = control_genotypes.count("hom_alt")
    return n_het <= MAX_HET_ALLOWED_IN_CONTROLS and n_hom_alt <= MAX_HOM_ALT_ALLOWED_IN_CONTROLS


def passes_roh_check(variant: dict, roh_regions: list) -> bool:
    """Variant position must fall inside one of Step 1's shared homozygous regions."""
    position = int(variant["position_bp"])
    for region in roh_regions:
        if variant["chrom"] == region["chrom"] and region["start_bp"] <= position <= region["end_bp"]:
            return True
    return False


def passes_consequence_check(variant: dict) -> bool:
    """Only variants predicted to change the protein sequence are plausible candidates."""
    return variant["consequence"] == "missense"


def main():
    print("STEP 2: Recessive-model candidate filtering")
    print(f"  reading {CANDIDATES_FILE.relative_to(REPO_ROOT)}, "
          f"{CONTROL_PANEL_FILE.relative_to(REPO_ROOT)}, "
          f"{ROH_FILE.relative_to(REPO_ROOT)}")

    roh_regions = load_roh_regions()
    control_genotypes_by_variant = load_control_genotypes()
    with open(CANDIDATES_FILE, newline="") as f:
        candidates = list(csv.DictReader(f))

    survivors = []
    for variant in candidates:
        checks = {
            "recessive (both cases hom_alt)": passes_case_genotype_check(variant),
            "rare in controls": passes_control_frequency_check(
                control_genotypes_by_variant[variant["variant_id"]]
            ),
            "inside case-specific ROH": passes_roh_check(variant, roh_regions),
            "missense (protein-changing)": passes_consequence_check(variant),
        }
        passed_all = all(checks.values())
        status = "KEEP" if passed_all else "drop"
        failed = [name for name, ok in checks.items() if not ok]
        reason = "" if passed_all else f" (failed: {', '.join(failed)})"
        print(f"  [{status}] {variant['variant_id']} ({variant['gene']}, "
              f"{variant['consequence']}){reason}")
        if passed_all:
            survivors.append(variant)

    OUTPUT_FILE.parent.mkdir(exist_ok=True)
    with open(OUTPUT_FILE, "w", newline="") as f:
        fieldnames = list(candidates[0].keys())
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(survivors)

    print(f"  {len(survivors)}/{len(candidates)} candidates survived all four filters")
    print(f"  wrote {OUTPUT_FILE.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
