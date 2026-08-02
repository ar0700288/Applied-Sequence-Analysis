"""
Stage 4 - REAL comparison. Checks whether the pipeline's own top candidate matches the
variant actually reported in the paper's abstract and Results.

Input:  ../stage3_pathogenicity/ranked_candidates.csv
Output: validation_report.txt
"""

import csv
from pathlib import Path

STAGE_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = STAGE_DIR.parent
INPUT_FILE = PIPELINE_ROOT / "stage3_pathogenicity" / "ranked_candidates.csv"
OUTPUT_FILE = STAGE_DIR / "validation_report.txt"

TARGET_VARIANT = {
    "gene": "LOXHD1", "chrom": "chr7", "position": "44806821",
    "ref": "G", "alt": "C", "protein_change": "p.(G1914A)",
}


def main():
    print("STAGE 4: validate against the paper's reported result")
    with open(INPUT_FILE, newline="") as f:
        top = next(iter(csv.DictReader(f)))

    checks = {
        "gene": top["gene"] == TARGET_VARIANT["gene"],
        "chromosome": top["chrom"] == TARGET_VARIANT["chrom"],
        "position": top["position"] == TARGET_VARIANT["position"],
        "reference allele": top["ref"] == TARGET_VARIANT["ref"],
        "alternate allele": top["alt"] == TARGET_VARIANT["alt"],
    }
    passed = all(checks.values())

    lines = [
        "LOXHD1 reproduction pipeline - validation report",
        "=" * 50,
        f"Paper's reported causal variant: {TARGET_VARIANT['gene']} "
        f"{TARGET_VARIANT['chrom']}:{TARGET_VARIANT['position']} "
        f"{TARGET_VARIANT['ref']}>{TARGET_VARIANT['alt']} ({TARGET_VARIANT['protein_change']})",
        f"Pipeline's top candidate:        {top['gene']} {top['chrom']}:{top['position']} "
        f"{top['ref']}>{top['alt']}",
        "",
    ]
    for name, ok in checks.items():
        lines.append(f"  [{'match' if ok else 'MISMATCH'}] {name}")
    lines.append("")
    lines.append("RESULT: " + ("PASS - matches the paper's reported variant" if passed
                                else "FAIL - does not match"))

    report = "\n".join(lines)
    print()
    print(report)
    OUTPUT_FILE.parent.mkdir(exist_ok=True)
    OUTPUT_FILE.write_text(report + "\n")
    print(f"\n  wrote {OUTPUT_FILE.relative_to(PIPELINE_ROOT)}")

    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
