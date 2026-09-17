"""
Stage 2 - MIXED (real VCF + 2 small fake files). Check each real candidate variant
against three rules from the paper's "Variant analysis" section:

  1. RECESSIVE  - both sequenced dogs must be homozygous (1/1) in the real VCF.
  2. RARE IN CONTROLS - <=2 heterozygous carriers, 0 homozygous, in the (fake) control
     panel. No real control-genome panel is public, so this input is fabricated - see
     ../fake_control_panel.csv for why.
  3. INSIDE A ROH REGION - the variant's position must fall inside a (fake) shared
     homozygous region. No real SNP-array data is public, so this input is also
     fabricated, and was never computed by any algorithm - see ../fake_roh_regions.csv.

A fourth fact - whether the variant is missense (protein-changing) - is not computed
either. It's copied directly from the paper's own reported annotation (Ensembl/NCBI, run
by the original authors), not re-derived here.

Because both candidate variants in the real VCF are already known finalists from the
paper, this stage does not discover anything - it confirms both pass every check, for
specific, checkable reasons. It cannot demonstrate rejecting a bad candidate, because
there is no bad candidate in the real input.

Input:  ../stage1_real_evidence/real_variants.vcf
        ../fake_control_panel.csv
        ../fake_roh_regions.csv
Output: filtered_candidates.csv
"""

import csv
from pathlib import Path

STAGE_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = STAGE_DIR.parent
VCF_FILE = PIPELINE_ROOT / "stage1_real_evidence" / "real_variants.vcf"
CONTROL_PANEL_FILE = PIPELINE_ROOT / "fake_control_panel.csv"
ROH_FILE = PIPELINE_ROOT / "fake_roh_regions.csv"
OUTPUT_FILE = STAGE_DIR / "filtered_candidates.csv"

MAX_HET_ALLOWED = 2
MAX_HOM_ALT_ALLOWED = 0

# Copied directly from the paper's Results section - not computed by any annotation
# tool here. Real annotation (Ensembl release 100 + NCBI Annotation Release 105) was
# run by the original authors; this project does not re-run it.
KNOWN_CONSEQUENCE = {
    "LOXHD1": "missense",
    "MROH8": "missense",
}


def load_vcf(path: Path):
    """Minimal VCF reader - just enough for this project's 2-row files."""
    variants = []
    sample_names = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                sample_names = line.split("\t")[9:]
                continue
            fields = line.split("\t")
            chrom, pos, _id, ref, alt, _qual, _filter, info = fields[:8]
            gene = dict(item.split("=") for item in info.split(";"))["GENE"]
            genotypes = {}
            for name, sample_field in zip(sample_names, fields[9:]):
                gt = sample_field.split(":")[0]
                genotypes[name] = gt
            variants.append({
                "chrom": chrom, "position": int(pos), "ref": ref, "alt": alt,
                "gene": gene, "genotypes": genotypes,
            })
    return variants


def load_csv_skipping_comments(path: Path):
    lines = [l for l in path.read_text().splitlines() if not l.startswith("#") and l.strip()]
    return list(csv.DictReader(lines))


def passes_recessive_check(variant: dict) -> bool:
    return all(gt == "1/1" for gt in variant["genotypes"].values())


def passes_control_frequency_check(variant: dict, control_rows: list) -> bool:
    row = next((r for r in control_rows if r["gene"] == variant["gene"]), None)
    if row is None:
        return False
    return (int(row["n_het"]) <= MAX_HET_ALLOWED
            and int(row["n_hom_alt"]) <= MAX_HOM_ALT_ALLOWED)


def passes_roh_check(variant: dict, roh_rows: list) -> bool:
    for region in roh_rows:
        if (variant["chrom"] == region["chrom"]
                and int(region["start_bp"]) <= variant["position"] <= int(region["end_bp"])):
            return True
    return False


def passes_consequence_check(variant: dict) -> bool:
    return KNOWN_CONSEQUENCE.get(variant["gene"]) == "missense"


def main():
    print("STAGE 2 (real VCF + fake control panel + fake ROH regions): filtering")
    variants = load_vcf(VCF_FILE)
    control_rows = load_csv_skipping_comments(CONTROL_PANEL_FILE)
    roh_rows = load_csv_skipping_comments(ROH_FILE)

    survivors = []
    for variant in variants:
        checks = {
            "recessive (both dogs 1/1)": passes_recessive_check(variant),
            "rare in controls (fake)": passes_control_frequency_check(variant, control_rows),
            "inside ROH region (fake)": passes_roh_check(variant, roh_rows),
            "missense (from paper's annotation)": passes_consequence_check(variant),
        }
        passed_all = all(checks.values())
        status = "KEEP" if passed_all else "drop"
        failed = [name for name, ok in checks.items() if not ok]
        reason = "" if passed_all else f" (failed: {', '.join(failed)})"
        print(f"  [{status}] {variant['gene']} {variant['chrom']}:{variant['position']}"
              f"{reason}")
        if passed_all:
            survivors.append(variant)

    OUTPUT_FILE.parent.mkdir(exist_ok=True)
    with open(OUTPUT_FILE, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["gene", "chrom", "position", "ref", "alt"])
        for v in survivors:
            writer.writerow([v["gene"], v["chrom"], v["position"], v["ref"], v["alt"]])

    print(f"  {len(survivors)}/{len(variants)} candidates survived")
    print(f"  wrote {OUTPUT_FILE.relative_to(PIPELINE_ROOT)}")


if __name__ == "__main__":
    main()
