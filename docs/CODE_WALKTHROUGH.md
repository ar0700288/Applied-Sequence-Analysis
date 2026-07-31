# Code walkthrough

This explains every file in the project: what it contains, why it's built the way it is,
and how it maps back to the paper. Read [`README.md`](../README.md) first for the
overall picture; this doc is the detailed follow-up.

## The big picture: how data flows through the pipeline

```
data/snp_genotypes.csv ──────► step1_homozygosity_mapping.py ──► results/step1_roh_regions.csv
                                                                          │
data/candidate_variants.csv ─┐                                          │
data/control_panel.csv ──────┼──► step2_filter_candidates.py ◄──────────┘
                              │            │
                              │            ▼
                              │    results/step2_filtered_candidates.csv
                              │            │
                              │            ▼
                              │    step3_score_pathogenicity.py
                              │            │
                              │            ▼
                              │    results/step3_ranked_candidates.csv
                              │            │
                              │            ▼
                              └──►  step4_validate_result.py
                                           │
                                           ▼
                                results/validation_report.txt
```

Every arrow is a CSV file, not a shared Python object or database. This is deliberate:
it's how real bioinformatics pipelines work (each tool reads files, writes files, and
knows nothing about the tool before or after it), and it means you can open any
`results/*.csv` in a spreadsheet and see exactly what one stage handed to the next.

---

## `data/snp_genotypes.csv`

This is the toy stand-in for the CanineHD SNP-array output from the paper's Homozygosity
mapping section: **6 dogs (3 affected "cases", 3 unaffected "controls") genotyped at a
set of markers along the genome.** A real array has 173,662 markers genome-wide; this
toy file has 48, placed on 4 fake chromosome regions, each one built to demonstrate a
specific rule the real algorithm has to get right:

| Region | What it contains | What it's there to prove |
|---|---|---|
| `chr7:41.2M-45.4M` | all 3 cases homozygous `CC`, all 3 controls het `CG` | this is the region that really matters - it's built to sit exactly where the paper says the real ROH was (chr7:41.2-45.5 Mb, Fig. 2b), and it's what step 1 should detect |
| `chr24:25.0M-26.5M` | same pattern, cases `AA` / controls `AG` | a second, smaller true ROH (stands in for one of the paper's other 21 ROH regions) |
| `chr1:10.0M-11.75M` | cases genotypes deliberately don't match each other | proves the algorithm requires cases to share the *same* homozygous genotype, not just *any* homozygous genotype each |
| `chr5:5.0M-5.75M` | cases **and** controls all homozygous `TT` | proves the algorithm correctly rejects a region just because it's homozygous - it also has to be *different from the controls*, otherwise it's just a common fixed spot in the breed, not disease-specific |

Every value was written by hand (not randomly generated), so each one is traceable back
to a rule in the paper rather than being an opaque fixture.

## `data/candidate_variants.csv`

Stand-in for the variant list that comes out of WES/WGS + read alignment + variant
calling (paper steps 6-9 in `WORKFLOW.md`) - the DNA differences found by sequencing the
two affected dogs, before any filtering. Real column meanings:

- `seq_case_wes_gt` / `seq_case_wgs_gt` - genotype call in the exome-sequenced dog and
  the genome-sequenced dog (`hom_alt` = homozygous for the variant, matching a real
  VCF's `1/1` call).
- `consequence` - what kind of change this is (`missense` = changes an amino acid;
  `intronic`, `UTR`, `synonymous` = doesn't change the protein sequence).

10 rows, each earning its place for a specific reason:

1. **`chr7_44806821` (LOXHD1, missense)** - the real variant from the paper. Should
   survive every filter.
2. **`chr24_25785932` (MROH8, missense)** - the real second candidate from the paper
   (the one ultimately ruled out). Should also survive Step 2's filters, then lose in
   Step 3's ranking.
3-7. **Non-missense variants inside the real ROH blocks** (`SLC26A9`, `ANKRD11`, `TMC1`,
   a second `LOXHD1` variant, a second `MROH8` variant) - these represent the other 5
   of the paper's "seven case-specific variants that resided in ROH" that weren't
   exonic/missense. They're here to prove Step 2's consequence filter actually does
   something, rather than everything conveniently already being a missense variant.
8. **`chr7_43700000` (COL4A3, missense)** - inside the real ROH block, and missense, but
   deliberately **common in the control panel** (see below). Proves the control-frequency
   filter, not just the ROH filter, is doing real work.
9. **`chr2_20000000` (RANDOMX, missense)** - on a chromosome that has no detected ROH at
   all. Proves variants outside any shared region get rejected regardless of anything
   else about them.
10. **`chr7_41000000` (RANDOMY, missense)** - on chr7, but 200 kb *before* the ROH block
    starts. A boundary case: proves the ROH-membership check is a real position
    comparison, not "same chromosome as something interesting."

## `data/control_panel.csv`

Stand-in for genotyping the same 10 variant sites in a panel of unaffected control dogs
(the paper used 637 real control genomes; this uses 10 toy ones - same rule, smaller
number). For every variant except one, only 1-2 controls are heterozygous carriers and
none are homozygous - i.e. the variant is rare, consistent with a recessive disease
allele. The exception is `chr7_43700000` (row 8 above), which has 3 heterozygous
controls **and** 1 homozygous-affected control - i.e. it's common enough in healthy dogs
that it can't be the cause of a rare disease. That's what step 2's control-frequency
filter is designed to catch.

---

## `scripts/step1_homozygosity_mapping.py`

Reimplements PLINK's `--homozyg-group`, simplified to the essential rule.

- `is_homozygous(genotype)` - one-line helper: `"CC"[0] == "CC"[1]` → `True`.
- `marker_qualifies(row)` - the heart of the algorithm. For one SNP marker, checks (a)
  all 3 cases are homozygous, (b) all 3 cases have the *same* genotype, and (c) the
  controls do **not** all share that same genotype. Only markers passing all three can
  be part of a shared, disease-specific ROH.
- `load_markers()` - reads the CSV and groups rows by chromosome (a run of homozygosity
  is a concept that only makes sense within one chromosome), sorting each chromosome's
  markers by position.
- `find_roh_runs(markers_by_chrom)` - walks each chromosome's markers in order, and
  every time it sees a `marker_qualifies() == True` immediately followed by another one,
  they're part of the same run; the first non-qualifying marker (or the end of the
  chromosome) closes the run. `MIN_MARKERS_PER_ROH = 4` filters out short runs that could
  just be coincidence - this is the toy-scale version of PLINK's `--homozyg-snp`
  minimum-marker-count parameter.
- `main()` - orchestrates the above, prints a per-chromosome summary of how many markers
  qualified (so you can *see* chr1 and chr5 getting correctly rejected), and writes the
  surviving regions to `results/step1_roh_regions.csv`.

## `scripts/step2_filter_candidates.py`

Reimplements the recessive-inheritance filtering described in the paper's "Variant
analysis" section, plus the ROH intersection from the Results section. Four independent
checks, each its own function so each one can be tested/read/explained on its own:

- `passes_case_genotype_check(variant)` - both sequenced dogs must be `hom_alt`.
- `passes_control_frequency_check(control_genotypes)` - counts `"het"` and `"hom_alt"`
  in the control list and compares against `MAX_HET_ALLOWED_IN_CONTROLS = 2` and
  `MAX_HOM_ALT_ALLOWED_IN_CONTROLS = 0` - literally the paper's own stated rule
  ("allowing a maximum of two heterozygous calls for each variant in the controls").
- `passes_roh_check(variant, roh_regions)` - reads `results/step1_roh_regions.csv` (the
  *output* of step 1 becomes step 2's *input* - this is the file-based hand-off shown in
  the diagram above) and checks whether the variant's position falls inside any region's
  `[start_bp, end_bp]` interval.
- `passes_consequence_check(variant)` - `consequence == "missense"`.

`main()` runs all four checks per variant, prints a `[KEEP]`/`[drop]` line explaining
*which* check failed when a variant is dropped (this is the most useful debugging output
in the whole pipeline - if a variant you expect to survive doesn't, this line tells you
exactly why), and writes survivors to `results/step2_filtered_candidates.csv`.

## `scripts/step3_score_pathogenicity.py`

This is the one stage that *can't* be meaningfully computed from scratch - PROVEAN and
PolyPhen-2 are trained models over large sequence databases, not something to
reimplement in a script. So this script does what happens in practice: it looks up the
value a scientist would get back from those tools. `PATHOGENICITY_SCORES` and
`KNOWN_HEARING_GENE` are literal facts taken from the paper's Results and Discussion
sections (real PROVEAN/PolyPhen-2 scores, and whether each gene had a prior known link to
hearing).

- `score_candidate(variant)` - looks up the two lookups above for one variant, and
  combines them into a simple, readable point total (`points`): 1 point for a
  "deleterious" PROVEAN call, 1 for a damaging PolyPhen-2 call, 2 for being a gene
  already known to cause hearing loss (weighted higher because it's the strongest,
  least-noisy signal - a plausible mechanism, not just a statistical prediction).
- `main()` - scores every survivor from step 2, sorts by score (highest first), prints
  the ranking, and writes it to `results/step3_ranked_candidates.csv`. The top row is the
  pipeline's final answer.

## `scripts/step4_validate_result.py`

The sanity check. `TARGET_VARIANT` is the exact variant reported in the paper's
abstract and Results (`chr7:44,806,821 G>C`, `LOXHD1`, `p.(G1914A)`). `main()` reads the
step-3 ranking, takes the #1 row, and compares each field (gene, chromosome, position,
ref allele, alt allele) against `TARGET_VARIANT`, printing `match`/`MISMATCH` per field
and an overall `PASS`/`FAIL`. It also `raise SystemExit(1)` on failure, so if you wire
this into a CI check later, a wrong answer makes the whole pipeline fail loudly instead
of silently writing a report nobody reads.

## `run_pipeline.py`

The only file that isn't itself one pipeline stage. It just calls the four
`scripts/step*.py` files in order using `subprocess.run`, in a fresh Python process each
time (mirroring how real workflow tools like Snakemake or Nextflow launch each stage as
its own process), and stops immediately if any stage exits with a non-zero return code -
which `step4` does on purpose if validation fails, so a broken pipeline can't silently
report success.

---

## If you wanted to make this "real"

The honest next step, in order of effort, would be:

1. Swap `data/snp_genotypes.csv` for a real PLINK `.bed/.bim/.fam` triplet and call
   actual `plink --homozyg-group` instead of `step1`'s hand-rolled version.
2. Download the real public reads (SRA BioProject `PRJNA702911`) and run BWA → GATK4 →
   DELLY as described in `WORKFLOW.md` steps 6-11, instead of hand-written
   `candidate_variants.csv`.
3. Get access to a real control-genome panel (the paper's is private/DBVDC) to replace
   `control_panel.csv`.
4. Submit the real candidate sequences to PROVEAN/PolyPhen-2/BLASTp instead of the
   hardcoded lookup tables in `step3`.

None of that changes the *shape* of the pipeline - it's still homozygosity mapping →
recessive filtering → pathogenicity ranking → validation. That shape is the actual thing
worth understanding and being able to explain.
