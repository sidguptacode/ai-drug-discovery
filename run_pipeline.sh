#!/usr/bin/env bash
# Do NOT run all pipeline steps in a row. Steps must be run one at a time and evaluated after each.
# Use: bash run_step_1.sh, bash run_step_2.sh, ... bash run_step_10.sh
# After each step: check outputs, decide if results are reasonable, then proceed or re-run.

echo "Do not run the pipeline as a single batch." >&2
echo "Run one step at a time with run_step_N.sh (N=1..10), then evaluate outputs before the next step." >&2
echo "" >&2
echo "  bash run_step_1.sh   # QC" >&2
echo "  bash run_step_2.sh   # Integration" >&2
echo "  bash run_step_3.sh   # Clustering" >&2
echo "  bash run_step_4.sh   # Annotation" >&2
echo "  bash run_step_5.sh   # Export" >&2
echo "  bash run_step_6.sh   # Load samples (Python)" >&2
echo "  bash run_step_7.sh   # Preprocess" >&2
echo "  bash run_step_8.sh   # LR scoring" >&2
echo "  bash run_step_9.sh   # CCI" >&2
echo "  bash run_step_10.sh  # Aggregate and rank" >&2
echo "" >&2
echo "Python step scripts (6–10) activate conda env 'ai-drug-discovery'; R steps use default R. Config: config.yml or PIPELINE_STEP_CONFIG." >&2
exit 1
