"""
Step 4 - Validate the result against the paper.

Real-world equivalent: everything after the computational analysis in the paper - PCR +
Sanger sequencing of 585 Rottweilers, conservation analysis across 99 species, and
population screening of >800,000 dogs - all of which independently confirmed that
chr7:44,806,821 G>C in LOXHD1 is the real, causal variant (see WORKFLOW.md steps 15-17).

This toy pipeline obviously can't run PCR or Sanger-sequence real dogs. What it CAN do is
the one thing that's actually meaningful for a from-scratch reproduction: check whether
the pipeline's own reasoning (Steps 1-3, run independently on synthetic data built from
the paper's numbers) lands on the same answer the paper reports. That's the real
pass/fail signal for "did I understand the method well enough to rebuild its logic."

Input:  results/step3_ranked_candidates.csv
Output: results/validation_report.txt
"""

import csv
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
INPUT_FILE = REPO_ROOT / "results" / "step3_ranked_candidates.csv"
OUTPUT_FILE = REPO_ROOT / "results" / "validation_report.txt"

# The variant actually reported in the paper (Results section):
# "a G>C missense variant at chr7:44,806,821 in lipoxygenase homology domains 1
#  (LOXHD1) ... predicted to result in a glycine-to-alanine substitution p.(G1914A)"
TARGET_VARIANT = {
    "gene": "LOXHD1",
    "chrom": "chr7",
    "position_bp": "44806821",
    "ref": "G",
    "alt": "C",
    "protein_change": "p.(G1914A)",
}


def main():
    print("STEP 4: Validate against the paper's reported result")
    print(f"  reading {INPUT_FILE.relative_to(REPO_ROOT)}")

    with open(INPUT_FILE, newline="") as f:
        ranked_candidates = list(csv.DictReader(f))
    top_candidate = ranked_candidates[0]

    checks = {
        "gene": top_candidate["gene"] == TARGET_VARIANT["gene"],
        "chromosome": top_candidate["chrom"] == TARGET_VARIANT["chrom"],
        "position": top_candidate["position_bp"] == TARGET_VARIANT["position_bp"],
        "reference allele": top_candidate["ref"] == TARGET_VARIANT["ref"],
        "alternate allele": top_candidate["alt"] == TARGET_VARIANT["alt"],
    }
    passed = all(checks.values())

    lines = []
    lines.append("LOXHD1 reproduction pipeline - validation report")
    lines.append("=" * 50)
    lines.append(f"Paper's reported causal variant: {TARGET_VARIANT['gene']} "
                 f"{TARGET_VARIANT['chrom']}:{TARGET_VARIANT['position_bp']} "
                 f"{TARGET_VARIANT['ref']}>{TARGET_VARIANT['alt']} "
                 f"({TARGET_VARIANT['protein_change']})")
    lines.append(f"Pipeline's top candidate:        {top_candidate['gene']} "
                 f"{top_candidate['chrom']}:{top_candidate['position_bp']} "
                 f"{top_candidate['ref']}>{top_candidate['alt']}")
    lines.append("")
    for name, ok in checks.items():
        lines.append(f"  [{'match' if ok else 'MISMATCH'}] {name}")
    lines.append("")
    lines.append("RESULT: " + ("PASS - pipeline recovered the paper's reported variant"
                                if passed else
                                "FAIL - pipeline's top candidate does not match the paper"))

    report = "\n".join(lines)
    print()
    print(report)

    OUTPUT_FILE.parent.mkdir(exist_ok=True)
    OUTPUT_FILE.write_text(report + "\n")
    print(f"\n  wrote {OUTPUT_FILE.relative_to(REPO_ROOT)}")

    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
