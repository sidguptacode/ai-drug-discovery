#!/bin/bash
#SBATCH --job-name=ov1_run01
#SBATCH --nodelist=cpunode4
#SBATCH --partition=cpunodes
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=5:00:00
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
    --message "Initialize a run named ov1_run01 from configs/config_ovarian1.yml, then run through the pipeline, evaluating and iterating based on the output of each step. \n\n
    Below is some additional context about the task and disease as well as guidelines to follow when reflecting and iterating on certain steps. \n 
    **Data Background: **\n
    High-grade serous ovarian carcinoma (HGSOC) is the most common and lethal subtype of ovarian cancer, characterised by near-universal TP53 mutations, widespread copy number alterations (CNAs), and a high degree of intratumoural heterogeneity. HGSOC tumours are polyclonal, harbouring spatially segregated subclones with distinct CNA profiles associated with chemotherapy resistance and poor prognosis. Tumour cells actively interact with a complex tumour microenvironment (TME) comprising immune cells (T cells, macrophages, plasma cells), fibroblasts, endothelial cells, and epithelial populations. These interactions, mediated through ligand-receptor (LR) signalling, are thought to shape local immune infiltration patterns and drive disease progression. Understanding which LR pairs mediate communication between HGSOC cell populations is a key step toward identifying novel therapeutic targets. \n\n
    **Reflection guidlines: **: \n
    **Annotation reflection: **\n
Do the labels make sense in the context of the **disease**, **data sample**, and **user’s target** (e.g. cell types expected in the tissue)? \n
Are EnrichR scores good enough (e.g. many clusters with adjusted P < 0.05)? If labels are generic or scores are poor, consider: \n
  - Changing 'annotation.primary_db' or 'annotation.enrichr_dbs' \n
  - Using 'label_filter_preset', 'label_prefer_patterns', or 'label_disqualify_patterns' in config to steer or filter terms \n 
Cross-check assigned labels against what is known from config: 'species' (labels should match the dataset species when the database mixes species), and tissue/region(labels should be plausible for the disease or sample context—e.g. brain disease → brain-relevant cell types; lung sample -> lung/airway plausible). \n
- If any label contradicts dataset identity (wrong species, or tissue/body region that does not fit the disease or sample), do not treat the run as acceptable. Remediate (e.g. label_prefer_patterns / label_disqualify_patterns to align with species and context, or change enrichr_dbs), re-run from step 4, and document in 'runs/<run_id>/run_info.json'.
    "
    \ --model "gpt-4o-2024-08-06"
    \ --instruction "ov1"
