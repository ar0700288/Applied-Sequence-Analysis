# LOXHD1 Rottweiler deafness — pipeline reproduction (learning project)

This is a small, dependency-free reimplementation of the **analytical logic** behind:

> Hytönen MK, Niskanen JE, Arumilli M, Brookhart-Knox CA, Donner J, Lohi H (2021).
> *Missense variant in LOXHD1 is associated with canine nonsyndromic hearing loss.*
> Human Genetics 140:1611-1618. https://doi.org/10.1007/s00439-021-02286-z

The study found that a Rottweiler litter's congenital deafness was caused by a single
missense variant, **chr7:44,806,821 G>C**, in the gene **LOXHD1**. They found it by
narrowing a whole genome down to one variant in four stages:

1. **Homozygosity mapping** — find genome regions where all affected dogs share the same
   homozygous DNA, and unaffected dogs don't (this is what "autosomal recessive" looks
   like in genotype data).
2. **Variant filtering** — from all the DNA variants found by sequencing, keep only the
   ones that (a) fall inside those shared regions, (b) are rare/absent in a healthy
   control population, and (c) actually change a protein.
3. **Pathogenicity scoring** — of the handful of variants left, rank them by how likely
   they are to break the protein and how plausible they are as a cause of *this* disease.
4. **Validation** — check that the winner matches what was actually reported.

**This repo reproduces exactly that four-stage funnel**, in plain Python, on a small
synthetic dataset built to mirror the real numbers in the paper. Given the toy inputs,
it independently re-derives the paper's answer: LOXHD1, chr7:44,806,821 G>C.

## What this is (and isn't)

This is a **learning / interview-prep project**, not a production bioinformatics
pipeline. Being upfront about that is more useful than pretending otherwise:

- The full study also involves whole-exome/genome sequencing, read alignment (BWA),
  variant calling (GATK4), structural-variant/mobile-element calling (DELLY/MELT),
  functional annotation (Ensembl/VEP), a private 637-genome control panel, PROVEAN/
  PolyPhen-2 web tools, BLASTp/COBALT conservation analysis, and Sanger sequencing across
  585 dogs plus a 771,864-dog commercial screen. Those need real multi-GB sequencing data,
  licensed software, and wet-lab work — see [`WORKFLOW.md`](WORKFLOW.md) for the complete
  17-step table of what real tool does what, with no shortcuts taken in the description.
- What's genuinely reproducible with light code — homozygosity mapping, recessive-model
  filtering, and pathogenicity-based ranking — is implemented here for real, on synthetic
  data sized so a first read-through takes minutes, not a GPU cluster.
- The synthetic data is built by hand (not randomly generated) so every input value is
  traceable to a specific line in the paper. See [`docs/CODE_WALKTHROUGH.md`](docs/CODE_WALKTHROUGH.md)
  for exactly how.

## Directory layout

```
loxhd1-repro/
├── README.md                       this file
├── WORKFLOW.md                     the full study workflow as a Step/Input/Output/Tools table
├── docs/
│   └── CODE_WALKTHROUGH.md         file-by-file, section-by-section explanation of the code
├── data/                           hand-built synthetic inputs (see walkthrough for how these map to the paper)
│   ├── snp_genotypes.csv           toy SNP-array genotypes for 6 dogs (3 affected, 3 unaffected)
│   ├── candidate_variants.csv      toy sequencing variants found in the two sequenced affected dogs
│   └── control_panel.csv           toy control-cohort genotypes at those same variant sites
├── scripts/                        one script per pipeline stage, each runnable on its own
│   ├── step1_homozygosity_mapping.py
│   ├── step2_filter_candidates.py
│   ├── step3_score_pathogenicity.py
│   └── step4_validate_result.py
├── results/                        generated when you run the pipeline (empty until then)
└── run_pipeline.py                 runs all four steps in order
```

## How to run it

No installation needed — everything uses only the Python standard library.

```bash
python run_pipeline.py
```

This runs the four steps in sequence, printing what each one is doing and why, and
writes its output into `results/`. The last line tells you whether the pipeline's own
top candidate matches the variant actually reported in the paper.

You can also run any single stage on its own, e.g. just the homozygosity mapping:

```bash
python scripts/step1_homozygosity_mapping.py
```

## Where to read next

- Read [`WORKFLOW.md`](WORKFLOW.md) first — it's the paper's whole methodology as a table,
  independent of this code, and it's the main artifact for explaining the study's
  workflow in an interview.
- Read [`docs/CODE_WALKTHROUGH.md`](docs/CODE_WALKTHROUGH.md) for what each script and
  each data file actually does, and why it's written the way it is.
