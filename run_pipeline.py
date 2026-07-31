"""
Single entry point: runs all four pipeline stages in order.

Each stage is also a normal, standalone script (see scripts/step*.py) that reads its
input file(s) and writes its output file(s) into results/. This script just calls them
in the right order and stops immediately if one of them fails - the same "files are the
interface between stages" pattern real pipeline tools (Snakemake, Nextflow, WDL) use,
just written out by hand instead of hidden behind a framework.

Usage:
    python run_pipeline.py
"""

import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent

PIPELINE_STEPS = [
    "scripts/step1_homozygosity_mapping.py",
    "scripts/step2_filter_candidates.py",
    "scripts/step3_score_pathogenicity.py",
    "scripts/step4_validate_result.py",
]


def main():
    for step_script in PIPELINE_STEPS:
        print("=" * 70, flush=True)
        result = subprocess.run([sys.executable, step_script], cwd=REPO_ROOT)
        if result.returncode != 0:
            print(f"\nPipeline stopped: {step_script} exited with code {result.returncode}")
            sys.exit(result.returncode)
        print(flush=True)

    print("=" * 70)
    print("Pipeline finished. See results/ for all intermediate and final outputs.")


if __name__ == "__main__":
    main()
