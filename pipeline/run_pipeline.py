"""
Runs all 4 stages in order. Each stage is also independently runnable on its own -
see the stageN_*/*.py file directly.

Usage:
    python run_pipeline.py
"""

import subprocess
import sys
from pathlib import Path

PIPELINE_ROOT = Path(__file__).resolve().parent

STAGES = [
    "stage1_real_evidence/pull_and_build_vcf.py",
    "stage2_filtering/filter_candidates.py",
    "stage3_pathogenicity/score_pathogenicity.py",
    "stage4_validation/validate_result.py",
]


def main():
    for stage_script in STAGES:
        print("=" * 70, flush=True)
        result = subprocess.run([sys.executable, stage_script], cwd=PIPELINE_ROOT)
        if result.returncode != 0:
            print(f"\nPipeline stopped: {stage_script} exited with code {result.returncode}")
            sys.exit(result.returncode)
        print(flush=True)

    print("=" * 70)
    print("Pipeline finished.")


if __name__ == "__main__":
    main()
