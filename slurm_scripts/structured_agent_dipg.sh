#!/bin/bash
#SBATCH --job-name=cluster_var_score_06
#SBATCH --nodelist=cpunode4
#SBATCH --partition=cpunodes
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/w/20251/cchin/ai-drug-discovery/logs/slurm_%j.out
#SBATCH --error=/w/20251/cchin/ai-drug-discovery/logs/slurm_%j.err

set -euo pipefail

PROJECT_DIR="/w/20251/cchin/ai-drug-discovery"
cd "$PROJECT_DIR"

# Load environment variables (OPENAI_API_KEY etc.)
if [ -f "$PROJECT_DIR/.env" ]; then
    set -a
    source "$PROJECT_DIR/.env"
    set +a
fi

source "$PROJECT_DIR/.venv/bin/activate"

python agent_controller.py \
    --message "Initialize a run named cluster_var_score_06 from config.yml, then run through the pipeline, evaluating and potentially iterating based on the output of each step. \n\n
    Below is some additional context about the task and disease. 
    Background:
    Diffuse Intrinsic Pontine Glioma (DIPG) in humans is an aggressive paediatric brain cancer arising in the brainstem, predominantly the pons. It carries no curative treatment options and a median survival of under one year. DIPG tumours are characterised by significant intratumoural cellular heterogeneity and a highly invasive phenotype. Tumour cells are known to actively interact with their surrounding tumour microenvironment (TME), which includes endothelial, neuronal, myeloid, and glial cell populations. These interactions are thought to drive disease progression and invasion. Understanding how DIPG tumour cells communicate with the TME through ligand-receptor (LR) signalling is a key step toward identifying novel therapeutic targets. \n\n
        Task: Ligand-Receptor Pair Identification from DIPG Spatial Transcriptomics Data" 
    \ --model "gpt-4o-2024-08-06"
