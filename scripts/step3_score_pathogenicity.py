"""
Step 3 - Pathogenicity scoring and ranking.

Real-world equivalent: PROVEAN + PolyPhen-2 scoring, plus the paper's own reasoning in
the Results section for why LOXHD1 was chosen over MROH8 once both had two variants left.

After Step 2 there are usually still a handful of missense variants left (in the paper:
exactly 2 - LOXHD1 and MROH8). Both COULD be damaging. To pick a final candidate, the
paper used two more pieces of evidence per variant:

  1. PREDICTED PROTEIN DAMAGE - PROVEAN and PolyPhen-2 are tools that estimate, from
     protein sequence alone, whether a substitution is likely to break the protein's
     function. These are not something you can meaningfully reimplement in a few lines
     of Python (they use trained statistical models over large sequence databases), so
     this script does what a scientist actually does with them: it looks up the value
     that tool reported and uses it. The numbers below are the REAL PROVEAN/PolyPhen-2
     scores reported in the paper for these two exact variants.
  2. PHENOTYPE PLAUSIBILITY - does the gene make biological sense as a cause of hearing
     loss? The paper notes LOXHD1 was already a known hearing-loss gene in humans and
     mice, whereas MROH8 had only been linked to unrelated traits (red blood cell
     volume, BMI, telomere length, hippocampal atrophy) - i.e. no prior connection to
     hearing at all. This isn't a number you compute; it's a literature fact you look up
     per candidate gene.

We combine both into a simple score and rank the candidates - exactly the reasoning
the paper's Discussion walks through in prose.

Input:  results/step2_filtered_candidates.csv
Output: results/step3_ranked_candidates.csv
"""

import csv
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
INPUT_FILE = REPO_ROOT / "results" / "step2_filtered_candidates.csv"
OUTPUT_FILE = REPO_ROOT / "results" / "step3_ranked_candidates.csv"

# Real PROVEAN / PolyPhen-2 results reported in the paper (Results section) for the
# two variants that survived Step 2. A variant not in this table would need to actually
# be submitted to those tools - there's no formula standing in for them here.
PATHOGENICITY_SCORES = {
    "chr7_44806821": {  # LOXHD1
        "provean_score": -4.517,
        "provean_call": "deleterious",       # PROVEAN calls anything <= -2.5 "deleterious"
        "polyphen_score": 0.992,
        "polyphen_call": "probably damaging",  # PolyPhen-2 HumVar: >0.85 = probably damaging
    },
    "chr24_25785932": {  # MROH8
        "provean_score": -3.336,
        "provean_call": "deleterious",
        "polyphen_score": 0.550,
        "polyphen_call": "possibly damaging",  # 0.15-0.85 = possibly damaging
    },
}

# Is this gene already known to matter for hearing, independent of this study?
# (From the paper's Discussion: LOXHD1 causes hearing loss in humans and mice;
# MROH8 has GWAS hits for unrelated traits and no known hearing role.)
KNOWN_HEARING_GENE = {
    "LOXHD1": True,
    "MROH8": False,
}


def score_candidate(variant: dict) -> dict:
    scores = PATHOGENICITY_SCORES[variant["variant_id"]]
    is_known_hearing_gene = KNOWN_HEARING_GENE.get(variant["gene"], False)

    # Simple, transparent point system - not a statistical model, just a readable
    # stand-in for "how many independent pieces of evidence point the same way".
    points = 0
    points += 1 if scores["provean_call"] == "deleterious" else 0
    points += 1 if scores["polyphen_call"] in ("probably damaging", "possibly damaging") else 0
    points += 2 if is_known_hearing_gene else 0  # prior biological plausibility counts double

    return {
        "variant_id": variant["variant_id"],
        "gene": variant["gene"],
        "chrom": variant["chrom"],
        "position_bp": variant["position_bp"],
        "ref": variant["ref"],
        "alt": variant["alt"],
        "provean_score": scores["provean_score"],
        "provean_call": scores["provean_call"],
        "polyphen_score": scores["polyphen_score"],
        "polyphen_call": scores["polyphen_call"],
        "known_hearing_gene": is_known_hearing_gene,
        "score": points,
    }


def main():
    print("STEP 3: Pathogenicity scoring and ranking")
    print(f"  reading {INPUT_FILE.relative_to(REPO_ROOT)}")

    with open(INPUT_FILE, newline="") as f:
        candidates = list(csv.DictReader(f))

    scored = [score_candidate(v) for v in candidates]
    scored.sort(key=lambda v: v["score"], reverse=True)

    for rank, candidate in enumerate(scored, start=1):
        print(f"  #{rank} {candidate['gene']} ({candidate['variant_id']}): "
              f"PROVEAN {candidate['provean_score']} ({candidate['provean_call']}), "
              f"PolyPhen-2 {candidate['polyphen_score']} ({candidate['polyphen_call']}), "
              f"known hearing gene: {candidate['known_hearing_gene']}, "
              f"score={candidate['score']}")

    OUTPUT_FILE.parent.mkdir(exist_ok=True)
    with open(OUTPUT_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(scored[0].keys()))
        writer.writeheader()
        writer.writerows(scored)

    top = scored[0]
    print(f"  top candidate: {top['gene']} ({top['variant_id']})")
    print(f"  wrote {OUTPUT_FILE.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
