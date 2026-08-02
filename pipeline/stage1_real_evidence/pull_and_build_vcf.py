"""
Stage 1 - REAL. Pull real reads from the real public sequencing data for this study, and
turn what they show into a real VCF file.

Data source: SRA BioProject PRJNA702911 (the actual dogs from this actual study).
    WES run (affected dog 1): SRR13743383
    WGS run (affected dog 2): SRR13743384

This does NOT scan the whole genome for unknown variants - it checks the 2 exact
positions the paper reports (LOXHD1, MROH8), the same way steps 2-4 in this doc's
docstrings assume as a starting point. Finding those 2 positions blind, from nothing,
would require running a real aligner + variant caller across the whole exome/genome -
a much bigger job (see the project's own notes on why that wasn't done here).

Dependency: NCBI's native-Windows SRA Toolkit (sam-dump.exe) - not available via conda on
Windows. Download: https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-win64.zip
Point the SRATOOLS_BIN environment variable at its bin/ folder before running this script.
If evidence/*.sam already exist (already pulled once), this script reuses them instead of
re-fetching, so it can be re-run offline.

Output: real_variants.vcf - a real, standard-format VCF with 2 rows (LOXHD1, MROH8),
one genotype column per sequenced dog, built directly from real read counts.
"""

import os
import re
import subprocess
from pathlib import Path

STAGE_DIR = Path(__file__).resolve().parent
EVIDENCE_DIR = STAGE_DIR / "evidence"
OUTPUT_VCF = STAGE_DIR / "real_variants.vcf"

SRATOOLS_BIN_DEFAULT = r"C:\sratoolkit\bin"
SRATOOLS_BIN = Path(os.environ.get("SRATOOLS_BIN", SRATOOLS_BIN_DEFAULT))

# The two real candidate variants reported in the paper's Results section.
TARGET_VARIANTS = [
    {"gene": "LOXHD1", "chrom": "chr7", "position": 44_806_821, "ref": "G", "alt": "C",
     "window_start": 44_806_700, "window_end": 44_806_950},
    {"gene": "MROH8", "chrom": "chr24", "position": 25_785_932, "ref": "C", "alt": "T",
     "window_start": 25_785_800, "window_end": 25_786_050},
]

# The two real sequenced affected dogs from this exact study.
SEQUENCED_DOGS = [
    {"label": "wes_case", "sra_run": "SRR13743383"},
    {"label": "wgs_case", "sra_run": "SRR13743384"},
]

EXCLUDE_FLAG_MASK = 0x904  # unmapped / secondary / supplementary
MIN_MAPPING_QUALITY = 20


def fetch_region_sam(sra_run: str, chrom: str, start: int, end: int, out_path: Path) -> None:
    """Ask NCBI's remote SRA index for just the reads overlapping one small window."""
    if out_path.exists():
        return  # already pulled once - reuse it, don't re-hit the network
    sam_dump_exe = SRATOOLS_BIN / "sam-dump.exe"
    if not sam_dump_exe.exists():
        raise FileNotFoundError(
            f"sam-dump.exe not found at {sam_dump_exe} and {out_path.name} doesn't exist "
            f"yet. Download NCBI's SRA Toolkit and set SRATOOLS_BIN, or supply the .sam "
            f"files directly in {EVIDENCE_DIR}."
        )
    region = f"{chrom}:{start}-{end}"
    result = subprocess.run(
        [str(sam_dump_exe), "--aligned-region", region, sra_run],
        capture_output=True, text=True, timeout=300,
    )
    if result.returncode != 0:
        raise RuntimeError(f"sam-dump failed for {sra_run} {region}: {result.stderr}")
    out_path.write_text(result.stdout)


def base_at_reference_position(pos: int, cigar: str, seq: str, target: int):
    """Walk one SAM record's CIGAR to find the read's letter at one reference position."""
    ref_pos = pos
    seq_pos = 0
    for length_str, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        length = int(length_str)
        if op in ("M", "=", "X"):
            if ref_pos <= target < ref_pos + length:
                return seq[seq_pos + (target - ref_pos)]
            ref_pos += length
            seq_pos += length
        elif op in ("I", "S"):
            seq_pos += length
        elif op in ("D", "N"):
            if ref_pos <= target < ref_pos + length:
                return "-"
            ref_pos += length
    return None


def pileup(sam_path: Path, target_position: int):
    allele_counts = {}
    n_reads_used = 0
    with open(sam_path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            flag, pos, mapq, cigar, seq = (
                int(fields[1]), int(fields[3]), int(fields[4]), fields[5], fields[9]
            )
            if flag & EXCLUDE_FLAG_MASK or mapq < MIN_MAPPING_QUALITY:
                continue
            base = base_at_reference_position(pos, cigar, seq, target_position)
            if base is None:
                continue
            n_reads_used += 1
            allele_counts[base] = allele_counts.get(base, 0) + 1
    return n_reads_used, allele_counts


def call_genotype(allele_counts: dict, ref: str, alt: str) -> str:
    """Returns a VCF-style genotype string: 0/0, 0/1, or 1/1."""
    total = sum(allele_counts.values())
    if total == 0:
        return "./."
    ref_frac = allele_counts.get(ref, 0) / total
    alt_frac = allele_counts.get(alt, 0) / total
    if alt_frac >= 0.85:
        return "1/1"
    if ref_frac >= 0.85:
        return "0/0"
    if alt_frac >= 0.2 and ref_frac >= 0.2:
        return "0/1"
    return "./."


def write_vcf(variant_calls: list) -> None:
    """variant_calls: list of dicts, one per variant, each with genotype calls per dog."""
    lines = [
        "##fileformat=VCFv4.2",
        "##source=pull_and_build_vcf.py (real read counts, SRA BioProject PRJNA702911)",
        '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene containing this variant (from paper annotation)">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for ref,alt (real read counts)">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(dog["label"] for dog in SEQUENCED_DOGS),
    ]
    for call in variant_calls:
        sample_fields = []
        for dog in SEQUENCED_DOGS:
            gt = call["genotypes"][dog["label"]]["gt"]
            ref_count = call["genotypes"][dog["label"]]["ref_count"]
            alt_count = call["genotypes"][dog["label"]]["alt_count"]
            sample_fields.append(f"{gt}:{ref_count},{alt_count}")
        lines.append(
            f"{call['chrom']}\t{call['position']}\t.\t{call['ref']}\t{call['alt']}\t.\t"
            f"PASS\tGENE={call['gene']}\tGT:AD\t" + "\t".join(sample_fields)
        )
    OUTPUT_VCF.write_text("\n".join(lines) + "\n")


def main():
    print("STAGE 1 (REAL): pulling real reads and building a real VCF")
    variant_calls = []
    for variant in TARGET_VARIANTS:
        print(f"  {variant['gene']} {variant['chrom']}:{variant['position']} "
              f"{variant['ref']}>{variant['alt']}")
        genotypes = {}
        for dog in SEQUENCED_DOGS:
            sam_path = EVIDENCE_DIR / f"{dog['label']}_{variant['gene'].lower()}.sam"
            fetch_region_sam(dog["sra_run"], variant["chrom"],
                              variant["window_start"], variant["window_end"], sam_path)
            n_reads, counts = pileup(sam_path, variant["position"])
            gt = call_genotype(counts, variant["ref"], variant["alt"])
            print(f"    {dog['label']} ({dog['sra_run']}): {n_reads} reads, {counts} -> {gt}")
            genotypes[dog["label"]] = {
                "gt": gt,
                "ref_count": counts.get(variant["ref"], 0),
                "alt_count": counts.get(variant["alt"], 0),
            }
        variant_calls.append({**variant, "genotypes": genotypes})

    write_vcf(variant_calls)
    print(f"  wrote {OUTPUT_VCF.relative_to(STAGE_DIR.parent.parent)}")


if __name__ == "__main__":
    main()
