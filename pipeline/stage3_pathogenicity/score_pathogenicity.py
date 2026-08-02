"""
Stage 3 - REAL (values copied from the paper, not recomputed). PROVEAN and PolyPhen-2 are
trained statistical models over large protein-sequence databases - not something to
reimplement in a script. So this stage does what happens in practice: it looks up the
value a scientist would get back from those tools. The numbers below are the actual
PROVEAN/PolyPhen-2 scores reported in the paper's Results for these exact two variants,
plus a fact from the paper's Discussion: LOXHD1 was already a known hearing-loss gene in
humans/mice, MROH8 was not.

Input:  ../stage2_filtering/filtered_candidates.csv
Output: ranked_candidates.csv
"""

import csv
from pathlib import Path

STAGE_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = STAGE_DIR.parent
INPUT_FILE = PIPELINE_ROOT / "stage2_filtering" / "filtered_candidates.csv"
OUTPUT_FILE = STAGE_DIR / "ranked_candidates.csv"

# Real scores reported in the paper's Results section.
PATHOGENICITY_SCORES = {
    "LOXHD1": {"provean_score": -4.517, "provean_call": "deleterious",
               "polyphen_score": 0.992, "polyphen_call": "probably damaging"},
    "MROH8": {"provean_score": -3.336, "provean_call": "deleterious",
              "polyphen_score": 0.550, "polyphen_call": "possibly damaging"},
}

# From the paper's Discussion: is this gene already known to cause hearing loss,
# independent of this study?
KNOWN_HEARING_GENE = {"LOXHD1": True, "MROH8": False}


def score_candidate(variant: dict) -> dict:
    scores = PATHOGENICITY_SCORES[variant["gene"]]
    known_hearing_gene = KNOWN_HEARING_GENE[variant["gene"]]
    points = 0
    points += 1 if scores["provean_call"] == "deleterious" else 0
    points += 1 if scores["polyphen_call"] in ("probably damaging", "possibly damaging") else 0
    points += 2 if known_hearing_gene else 0
    return {**variant, **scores, "known_hearing_gene": known_hearing_gene, "score": points}


def main():
    print("STAGE 3 (real published scores): pathogenicity ranking")
    with open(INPUT_FILE, newline="") as f:
        candidates = list(csv.DictReader(f))

    scored = [score_candidate(c) for c in candidates]
    scored.sort(key=lambda c: c["score"], reverse=True)

    for rank, c in enumerate(scored, start=1):
        print(f"  #{rank} {c['gene']}: PROVEAN {c['provean_score']} ({c['provean_call']}), "
              f"PolyPhen-2 {c['polyphen_score']} ({c['polyphen_call']}), "
              f"known hearing gene: {c['known_hearing_gene']}, score={c['score']}")

    OUTPUT_FILE.parent.mkdir(exist_ok=True)
    with open(OUTPUT_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(scored[0].keys()))
        writer.writeheader()
        writer.writerows(scored)

    print(f"  top candidate: {scored[0]['gene']}")
    print(f"  wrote {OUTPUT_FILE.relative_to(PIPELINE_ROOT)}")


if __name__ == "__main__":
    main()
