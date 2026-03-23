dipg_prompt="""
   Task: Ligand-Receptor Pair Identification from DIPG Spatial Transcriptomics Data
Background:
Diffuse Intrinsic Pontine Glioma (DIPG) is an aggressive paediatric brain cancer arising in the brainstem, predominantly the pons. It carries no curative treatment options and a median survival of under one year. DIPG tumours are characterised by significant intratumoural cellular heterogeneity and a highly invasive phenotype. Tumour cells are known to actively interact with their surrounding tumour microenvironment (TME), which includes endothelial, neuronal, myeloid, and glial cell populations. These interactions are thought to drive disease progression and invasion. Understanding how DIPG tumour cells communicate with the TME through ligand-receptor (LR) signalling is a key step toward identifying novel therapeutic targets.
Data:
You are provided with processed spatial transcriptomics count data from a DIPG-infiltrated human brainstem, spanning multiple sequential tissue regions from midbrain to spinal cord. The data is located at ./data/evals/dipg/spatial_transcriptomics_samples/ and contains 10 samples (sample1–sample10). Each sample directory contains three files: a filtered feature barcode matrix (*_filtered_feature_bc_matrix.h5), a tissue positions list (*_tissue_positions_list.csv) containing spatial (x, y) coordinates for each spot, and a scalefactors JSON (*_scalefactors_json.json) containing 10x Visium scale factors needed to interpret spatial coordinates. Load and integrate all 10 samples for analysis, retaining sample identity as a metadata field. 
Using the provided data, identify candidate ligand-receptor (LR) pairs that mediate intercellular communication between cell populations present in the tissue. Your final output should be a ranked list of the top candidate LR pairs, scored by a quantitative metric of your choice that reflects the strength or significance of the interaction. The output should be organised by the directional cell-type pairing involved (i.e. which cell population is the sender expressing the ligand, and which is the receiver expressing the receptor). Focus particularly on interactions between non-tumour TME populations and tumour cell populations.
Output Format:
For each directional cell-type pairing, return a ranked list of LR pairs in the following format:
Sender cell type → Receiver cell type

Ligand_Gene : Receptor_Gene | Score | Brief rationale
Ligand_Gene : Receptor_Gene | Score | Brief rationale
...

Cover all major sender-receiver pairings identified in your analysis, not only those you consider most biologically interesting. Document any analytical decisions you make, including how you identified cell types, how you scored LR pairs, and any filtering or quality control steps applied.
    """

ovarian_1_prompt="""
Task: Ligand-Receptor Pair Identification from HGSOC Spatial Transcriptomics Data

Background:
High-grade serous ovarian carcinoma (HGSOC) is the most common and lethal subtype of ovarian cancer, characterised by near-universal TP53 mutations, widespread copy number alterations (CNAs), and a high degree of intratumoural heterogeneity. HGSOC tumours are polyclonal, harbouring spatially segregated subclones with distinct CNA profiles associated with chemotherapy resistance and poor prognosis. Tumour cells actively interact with a complex tumour microenvironment (TME) comprising immune cells (T cells, macrophages, plasma cells), fibroblasts, endothelial cells, and epithelial populations. These interactions, mediated through ligand-receptor (LR) signalling, are thought to shape local immune infiltration patterns and drive disease progression. Understanding which LR pairs mediate communication between HGSOC cell populations is a key step toward identifying novel therapeutic targets.

Data:
You are provided with processed spatial transcriptomics count data from HGSOC tumour tissue sections profiled using 10x Visium. The data is located at ./data/evals/ovarian_1/ and contains 8 samples. All files are in the same directory and can be identified by the prefix SP1 through SP8. For each sample, the three files in 10x Visium sparse format are: GSM*_SP[N]_barcodes.tsv (spot barcodes), GSM*_SP[N]_features.tsv (gene identifiers), and GSM*_SP[N]_matrix.mtx (sparse count matrix), where [N] runs from 1 to 8. Load and integrate all 8 samples for analysis, retaining sample identity as a metadata field.

Task:
Using the provided data, identify candidate ligand-receptor (LR) pairs that mediate intercellular communication between all major cell populations present in the tissue. Your final output should be a ranked list of the top candidate LR pairs per directional cell-type pairing, scored by a quantitative metric of your choice that reflects the strength or significance of the interaction. Cover all major sender-receiver pairings identified in your analysis. Relevant non-malignant TME cell types to consider include, but are not limited to: macrophage subpopulations, T cells, fibroblasts, endothelial cells, and plasma cells. Pay particular attention to interactions between tumour and TME populations in both sending and receiving directions, as well as any autocrine loops within a single cell population.

Output Format:
For each directional cell-type pairing, return a ranked list of LR pairs in the following format:

Sender cell type -> Receiver cell type
Ligand_Gene : Receptor_Gene | Score | Brief rationale
Ligand_Gene : Receptor_Gene | Score | Brief rationale
...

Cover all major sender-receiver pairings identified in your analysis, not only those you consider most biologically interesting. Where an interaction involves an autocrine loop (sender and receiver are the same population), label it explicitly as [AUTOCRINE]. Document any analytical decisions you make, including how you identified cell types, how you scored LR pairs, and any filtering or quality control steps applied.
"""