# Workflow to regenerate the results of Hytönen et al. 2021 (LOXHD1 / Rottweiler deafness)

Paper: Hytönen MK, Niskanen JE, Arumilli M, Brookhart-Knox CA, Donner J, Lohi H (2021).
*Missense variant in LOXHD1 is associated with canine nonsyndromic hearing loss.*
Human Genetics 140:1611-1618. https://doi.org/10.1007/s00439-021-02286-z

This table is built directly from the paper's **Materials and methods** section, in the
order the authors actually performed the work: sample collection → array genotyping →
homozygosity mapping → sequencing → variant calling → filtering → in silico prediction →
wet-lab validation → population screening.

| # | Step | Process details | Input | Output | Tools |
|---|------|------------------|-------|--------|-------|
| 1 | Sample collection & DNA extraction | Collect EDTA blood from 585 Rottweilers (incl. 4 affected littermates); extract genomic DNA; QC concentration/purity | Blood samples | Genomic DNA aliquots | Chemagen extraction robot; NanoDrop ND-1000; Qubit 3.0 |
| 2 | SNP array genotyping | Genotype 4 affected + 2 unaffected relatives on a genome-wide SNP chip (173,662 markers) | Genomic DNA (6 dogs) | Raw genotype calls (PLINK bed/bim/fam) | Illumina CanineHD BeadChip |
| 3 | Genotype QC | Filter samples/markers: call rate >95%, HWE p < 1e-8; drop 1 poor-call-rate case | Raw genotype calls | QC'd genotype set: 3 cases, 3 controls, 154,235 markers | PLINK v1.9 |
| 4 | ROH parameter optimization | Simulate a fully homozygous genome; sweep `--homozyg-density` (10-125 kb/SNP) and `--homozyg-gap` (20-1000 kb) to find the values that maximize genome coverage | QC'd genotypes + simulated genome | Optimal ROH parameters (density = 30 kb/SNP, gap = 200 kb, min 70 SNPs/run) | PLINK 1.9; in-house Python (heterozygosity calc) |
| 5 | Homozygosity mapping | Detect runs of homozygosity (`--homozyg-group`) with optimized params; keep only ROH allelically shared by all 3 cases and different/absent in all 3 controls | QC'd genotypes + optimal params | 22 case-specific ROH regions, 62.3 Mb total (largest: chr1, 21.3 Mb) | PLINK 1.9 |
| 6 | Whole-exome sequencing (WES) | Capture + sequence exome of 1 affected dog, 2×150 bp, 38× coverage | Genomic DNA (1 case) | Raw FASTQ reads | NimbleGen SeqCap EZ kit; Illumina NextSeq500 |
| 7 | Whole-genome sequencing (WGS) | Sequence whole genome of a second affected littermate, 2×100 bp, 15× coverage | Genomic DNA (1 case) | Raw FASTQ reads | Illumina HiSeq2000 |
| 8 | Read alignment & processing | Align reads to CanFam3 reference; mark duplicates; sort/compress | FASTQ reads (WES + WGS) | Sorted, dedup BAM files | BWA v0.7.15; SpeedSeq; SAMBLASTER; Sambamba |
| 9 | Small-variant calling (SNV/indel) | Per-sample gVCF calling (HaplotypeCaller) → CombineGVCFs → joint genotyping | BAM files | Joint-genotyped VCF (SNVs + indels) | GATK 4.1 (HaplotypeCaller, GenotypeGVCFs) |
| 10 | Structural variant (SV) calling | Detect deletions/duplications/inversions/insertions in the WGS sample | WGS BAM | SV calls (VCF) | DELLY |
| 11 | Mobile element insertion (MEI) calling | Detect transposon insertions in the WGS sample against a Repbase reference | WGS BAM + Repbase library | MEI calls | MELT |
| 12 | Functional annotation | Annotate all SNV/indel/SV/MEI calls (exonic, intronic, splicing, UTR, other) | VCFs from steps 9-11 | Annotated variant tables | Ensembl release 100 + NCBI Annotation Release 105 (VEP-equivalent) |
| 13 | Recessive-model candidate filtering | Filter case variants against a control-genome panel (637 genomes for SNV/indel, 290 for SV/MEI); require homozygous in both sequenced cases, ≤2 heterozygous calls allowed in controls; intersect with case-specific ROH from step 5 | Annotated variants + control panel + ROH regions | 32 candidate SNV/indel (6 in ROH) + 63 SV (0 in ROH) + 32 MEI (1 in ROH) → 7 variants total in ROH → 2 exonic missense finalists (LOXHD1, MROH8) | webGQT variant server |
| 14 | In silico pathogenicity prediction | Score the two missense finalists for predicted functional impact | 2 candidate missense variants | PROVEAN scores (LOXHD1 −4.517, MROH8 −3.336, both "deleterious"); PolyPhen-2 scores (LOXHD1 0.992 "probably damaging", MROH8 0.550 "possibly damaging") | PROVEAN; PolyPhen-2 |
| 15 | Conservation / homology analysis | Retrieve LOXHD1 orthologs across 99 Eutherian species; align and check conservation of the affected residue (G1914) | LOXHD1 protein sequence (XP_022277134.1) | Multiple sequence alignment showing the residue is fully conserved | NCBI BLASTp; COBALT |
| 16 | Sanger validation & segregation | PCR-amplify the LOXHD1 site with flanking primers; Sanger-sequence 585 Rottweilers (4 cases + 581 relatives); check segregation with hearing status | Genomic DNA (585 dogs) | Confirmed genotypes; complete segregation with deafness; allele freq 2.6%, carrier freq 5.3% (excluding the affected family) | PCR (Taq polymerase); ABI 3730 capillary sequencer; Sequencher / UGENE |
| 17 | Population / breed screening | Screen the variant in a large commercial-testing cohort to assess breed specificity and general-population frequency | DNA from 28,116 dogs (374 breeds) + 771,864 mixed-breed dogs | Variant found only in Rottweiler-associated dogs; allele freq 0.04%, carrier freq 0.08% in the large cohort; 63.4% of carriers show Rottweiler ancestry | Genoscoper Laboratories genetic test; Wisdom Panel™ |

## Where this reproduction project stops

Steps 1-2, 6-7, 16-17 are wet-lab / commercial-database steps — no code reproduces them.
Step 15 is a one-off web-tool lookup. Steps 4, 5, 13, 14 are the purely **computational,
algorithmic** core of the paper — the part that determines *which variant gets picked* —
and that's what [`scripts/`](scripts/) reimplements, in simplified form, on a small
synthetic dataset so it runs instantly with no installs and no multi-GB downloads. See
[README.md](README.md) for exactly what is simplified and why, and
[docs/CODE_WALKTHROUGH.md](docs/CODE_WALKTHROUGH.md) for a line-by-line explanation of
the code.
