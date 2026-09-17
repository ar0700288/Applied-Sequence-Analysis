# LOXHD1 Rottweiler deafness, pipeline reproduction

This project tries to redo the reasoning behind one paper:

> Hytönen MK, Niskanen JE, Arumilli M, Brookhart-Knox CA, Donner J, Lohi H (2021).
> Missense variant in LOXHD1 is associated with canine nonsyndromic hearing loss.
> Human Genetics 140:1611-1618. https://doi.org/10.1007/s00439-021-02286-z

The paper found that a litter of deaf Rottweiler puppies all shared one broken gene,
LOXHD1, at the exact spot chr7:44,806,821, where a normal G was changed to a C. This
project checks that claim against real data where real data exists, and states plainly,
up front, where it doesn't.

## What the paper's full method actually looked like

Before getting into what this project does, here is the whole method the paper used, in
order, so it's clear what's being compared against.

| # | Step | What happened | Input | Output | Tools |
|---|------|----------------|-------|--------|-------|
| 1 | Sample collection & DNA extraction | Collect blood from 585 Rottweilers, including 4 affected littermates; extract DNA | Blood samples | Genomic DNA | Chemagen robot; NanoDrop; Qubit |
| 2 | SNP array genotyping | Genotype 4 affected + 2 unaffected dogs on a genome-wide SNP chip | Genomic DNA (6 dogs) | Raw genotype calls | Illumina CanineHD BeadChip |
| 3 | Genotype QC | Filter samples and markers by call rate and Hardy-Weinberg equilibrium | Raw genotype calls | QC'd genotypes, 3 cases + 3 controls | PLINK |
| 4 | ROH parameter optimization | Simulate a fully homozygous genome to tune detection settings | QC'd genotypes | Optimal ROH parameters | PLINK; in-house Python |
| 5 | Homozygosity mapping | Find regions shared by all affected dogs and different in controls | QC'd genotypes + parameters | 22 case-specific ROH regions | PLINK |
| 6 | Whole-exome sequencing | Sequence the exome of one affected dog | Genomic DNA | Raw FASTQ reads | NimbleGen kit; Illumina NextSeq500 |
| 7 | Whole-genome sequencing | Sequence the whole genome of a second affected dog | Genomic DNA | Raw FASTQ reads | Illumina HiSeq2000 |
| 8 | Read alignment | Align reads to the CanFam3 reference genome | FASTQ reads | Sorted BAM files | BWA; SAMBLASTER; Sambamba |
| 9 | Small-variant calling | Call SNVs and indels from the aligned reads | BAM files | Joint-genotyped VCF | GATK4 |
| 10 | Structural variant calling | Call larger deletions, duplications, inversions | WGS BAM | SV calls | DELLY |
| 11 | Mobile element calling | Call transposon insertions | WGS BAM | MEI calls | MELT |
| 12 | Functional annotation | Label each variant with its gene and consequence type | VCFs | Annotated variant list | Ensembl; NCBI annotation |
| 13 | Candidate filtering | Keep variants that are homozygous in both sequenced dogs, rare in a 637-genome control panel, and inside a shared ROH region | Annotated variants + control panel + ROH regions | 2 finalists: LOXHD1, MROH8 | webGQT |
| 14 | Pathogenicity prediction | Score how damaging each finalist's mutation looks | 2 candidate variants | PROVEAN and PolyPhen-2 scores | PROVEAN; PolyPhen-2 |
| 15 | Conservation analysis | Check whether the mutated position is unchanged across ~99 other species | LOXHD1 protein sequence | Confirms the position is conserved | NCBI BLASTp; COBALT |
| 16 | Sanger validation | Confirm the variant by direct sequencing across 585 dogs | Genomic DNA | Confirmed segregation with deafness | PCR; Sanger sequencing |
| 17 | Population screening | Check how common the variant is across breeds | DNA from 800,000+ dogs | Variant is Rottweiler-specific and rare | Commercial genetic testing |

## The first thing to know: two pieces of that method are private, so two files here are fake

Steps 2 through 5 and step 13's control panel both rely on data that was never made
public. The SNP-array genotypes for the study's own 6-dog cohort were never deposited
anywhere, and the 637-genome control panel belongs to a private consortium (DBVDC) that
never released it either. There is nothing to download for either one, no matter how
much time or compute is available.

Because of that, this project fills both gaps with small, clearly labeled fake files
instead of skipping them silently or faking them and pretending otherwise.

`pipeline/fake_roh_regions.csv` states, as plain numbers, the region boundaries the paper
itself reports. It does not run any algorithm to derive them, because there is no real
genotype data to run an algorithm on. The file's own header says this.

`pipeline/fake_control_panel.csv` states, as invented counts, how many "healthy" dogs
carry each candidate variant. There is no real control panel to draw real numbers from,
so these numbers were made up, specifically to exercise the "is this rare" check, not to
represent anything real. The file's own header says this too.

Steps 1, 16, and 17 need physical dogs and a lab, so nothing here reproduces those
either. Step 15 was left as a manual web-tool step in the real study and is not automated
here.

Everything else, the sequencing data itself, and the pathogenicity scores, uses real
public data or real values published in the paper, not recomputed. The rest of this
document explains which is which.

## Why it's split into separate stage folders instead of one script

Each stage of this project is its own small folder, with its own script, instead of one
long program. That mirrors how the table above is laid out: sequencing, then variant
calling, then filtering, then scoring, each one a distinct step that takes a specific
input and hands off a specific output to the next step. Splitting the code up the same
way means each stage can be read on its own, run on its own, and its output checked on
its own.

## Stage 1, real evidence (covers steps 6 through 9 above)

This is the one stage that touches real, live, external data. The paper's raw
sequencing reads for this exact study, from the two real affected dogs it sequenced, were
deposited publicly on NCBI's Sequence Read Archive, under BioProject PRJNA702911.

Downloading the full reads was not practical here: the exome run is about 3.4 GB, the
whole-genome run about 12.9 GB, and turning either into a usable aligned file on ordinary
hardware takes real time, on the order of a day for the exome data alone and several days
for the genome data, mainly because BWA and GATK have no native Windows build and this
machine is a two-core laptop, not a compute server. So rather than downloading everything
and aligning it from scratch, this stage asks NCBI directly for just the handful of reads
that overlap two exact positions: chr7:44,806,821, the real LOXHD1 variant, and
chr24:25,785,932, MROH8, the second candidate the paper's own text says it ruled out.
Those two positions were not discovered by this project. They were read out of the paper
and typed directly into the code, because scanning the whole genome blind for unknown
variants would require the full alignment and variant-calling run this stage specifically
avoids, for the resource reasons just described.

Concretely, the script asks NCBI's servers for only the reads in a roughly 250 base pair
window around each of those two positions, for both sequenced dogs, using the fact that
the original submission was already aligned to the CanFam3.1 reference genome when the
paper's authors uploaded it years ago. That keeps each download to a few kilobytes
instead of gigabytes. For every one of those four downloads, two dogs times two variants,
the script then reads through every individual DNA read, walks its alignment string to
work out exactly which letter that read has at the target position, and tallies up the
results. From those tallies it writes out a real VCF file, the standard format a variant
caller would produce, stating each dog's genotype at each position based only on what the
real reads actually show.

The result, right now: all four checks come back homozygous for the alternate allele, 34
out of 34 real reads and 16 out of 16 real reads showing the LOXHD1 C allele in the two
dogs, and 34 out of 35 and 14 out of 14 showing the MROH8 T allele, the one non-matching
read being an ordinary sequencing error rather than a second true allele. The raw
downloaded reads themselves are kept in `evidence/*.sam`, unedited, as the literal proof
of what was pulled.

## Stage 2, filtering (covers step 13 above)

This stage takes the real VCF from stage 1 and checks each of its two variants against
four rules taken from the paper's own filtering logic. Both dogs must be homozygous for
the variant, checked against the real VCF genotypes. The variant must be rare in healthy
dogs, checked against `fake_control_panel.csv`, since no real control panel exists. The
variant's position must fall inside a shared region, checked against
`fake_roh_regions.csv`, for the same reason. The variant must change the protein, which
is checked against a small list copied directly from the paper's own reported annotation
(step 12 above), not computed by any annotation tool here.

Because both variants entering this stage are already the paper's known final answers,
this stage cannot demonstrate rejecting a bad candidate. There is no bad candidate in the
input to reject. It only confirms that two already-correct answers pass four specific,
individually checkable rules.

## Stage 3, pathogenicity scoring (covers step 14 above)

PROVEAN and PolyPhen-2 are statistical tools trained on large protein databases, not
something that can be meaningfully rebuilt from scratch in a script. So this stage uses
the actual scores the paper reports for these two exact variants: PROVEAN calls both
deleterious, PolyPhen-2 calls LOXHD1 probably damaging and MROH8 only possibly damaging.
It also uses one more fact straight from the paper's discussion, that LOXHD1 was already
a known hearing-loss gene in humans and mice before this study, while MROH8 had no such
prior link. Combining those, LOXHD1 outranks MROH8, matching the paper's own conclusion.

## Stage 4, validation

The simplest stage. It takes the pipeline's own top answer from stage 3 and compares it,
field by field, gene name, chromosome, position, both DNA letters, against the exact
variant the paper reports in its abstract and results. This is not scientific evidence of
anything. It is a check that the code, run start to finish, produces the answer it is
supposed to, the same way a programmer writes a test that asserts an expected result.

Run `python pipeline/run_pipeline.py` and all four stages run in order, ending in a
report that currently says PASS, LOXHD1, chr7:44,806,821 G>C.

## Directory layout

```
loxhd1-repro-Project-01/
├── README.md                          this file
└── pipeline/
    ├── run_pipeline.py                 runs all 4 stages in order
    ├── fake_control_panel.csv          fake, no real control panel is public
    ├── fake_roh_regions.csv            fake, just states the paper's reported region
    │
    ├── stage1_real_evidence/
    │   ├── pull_and_build_vcf.py       pulls real reads, builds a real VCF
    │   ├── evidence/*.sam              the raw downloaded reads, kept as proof
    │   └── real_variants.vcf           generated output
    │
    ├── stage2_filtering/
    │   ├── filter_candidates.py
    │   └── filtered_candidates.csv     generated output
    │
    ├── stage3_pathogenicity/
    │   ├── score_pathogenicity.py
    │   └── ranked_candidates.csv       generated output
    │
    └── stage4_validation/
        ├── validate_result.py
        └── validation_report.txt      generated output
```

## How to run it

Stages 2 through 4 need nothing but plain Python, no installation. Stage 1 needs NCBI's
SRA Toolkit only if you want to re-download the reads yourself. The `.sam` files it needs
are already saved in `stage1_real_evidence/evidence/`, so by default it reuses those
instead of downloading anything.

```bash
python pipeline/run_pipeline.py
```

Or run one stage at a time, for example just the filtering step:

```bash
python pipeline/stage2_filtering/filter_candidates.py
```

To force stage 1 to actually re-fetch from NCBI, delete the `.sam` files in
`stage1_real_evidence/evidence/`, download the SRA Toolkit from
https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-win64.zip, and
point the SRATOOLS_BIN environment variable at its `bin` folder.

## License

The code in this project is released under the [MIT License](LICENSE). The SRA sequencing data and reference files keep their own terms of use.
