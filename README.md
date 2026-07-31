# VitisSeedOmics

This repository contains the principal computational workflows used in the study:

> **Multi-omics and explainable machine learning identify key enzymatic regulators of flavonoid biosynthesis in grape seeds**

The scripts cover graph-pangenome construction, pan-transcriptome indexing, RNA-seq processing, differential expression analysis, orthogroup classification, gene-structure comparison, and explainable machine learning with SHAP interaction analysis.
<img width="3482" height="2479" alt="Graphical abstract" src="https://github.com/user-attachments/assets/88b43b4e-53c0-4504-8a97-09702f62a6bb" />

## Repository contents

| Script | Purpose |
| --- | --- |
| `01_graphpangenome.sh` | Construct a grape graph pangenome with Cactus |
| `02_pantranscriptome.sh` | Build and index a VG pan-transcriptome |
| `03_rnaseq.sh` | Perform reference-based RNA-seq quality control, alignment, quantification, normalization, and QC |
| `04_DESeq_analysis.r` | Conduct two-group differential expression analysis with DESeq2 |
| `05_orthoFinder.sh` | Identify orthogroups and infer phylogenetic relationships with OrthoFinder |
| `06_classify_orthogroups.py` | Classify orthogroups as core, soft-core, shell, or cloud |
| `07_plot_genestructure.r` | Compare gene-structure characteristics among pan-genome categories |
| `08_SHAP_analysis.py` | Perform ensemble feature selection and XGBoost–SHAP interpretation |

## Workflow overview

The scripts form four related analysis branches:

1. **Graph-pangenome branch:** `01_graphpangenome.sh` → `02_pantranscriptome.sh`
2. **RNA-seq branch:** `03_rnaseq.sh` → `04_DESeq_analysis.r`
3. **Orthology branch:** `05_orthoFinder.sh` → `06_classify_orthogroups.py` → `07_plot_genestructure.r`
4. **Explainable machine-learning branch:** `08_SHAP_analysis.py`

## Software requirements

### Graph-pangenome and pan-transcriptome

- Cactus `v3.1.4`
- Singularity or Apptainer
- VG toolkit
- Toil, provided through the Cactus container
- SLURM-compatible high-performance computing environment

### RNA-seq analysis

- fastp
- HISAT2
- SAMtools
- Subread/featureCounts
- gffread
- MultiQC, optional
- R and Rscript
- R packages: `ggplot2`, `pheatmap`, `getopt`, `DESeq2`, `RColorBrewer`, and `grid`

### Orthogroup and gene-structure analysis

- OrthoFinder
- DIAMOND
- MAFFT
- IQ-TREE 3
- Python 3 with `pandas`
- R packages: `ggplot2`, `dplyr`, and `gridExtra`

### Explainable machine learning

- Python 3
- `numpy`
- `pandas`
- `matplotlib`
- `seaborn`
- `shap`
- `networkx`
- `statsmodels`
- `scikit-learn`
- `xgboost`
- `lightgbm`

The exact software versions used for a publication-quality reproduction should be recorded with the corresponding analysis results.

## Usage

### 1. Graph-pangenome construction

`01_graphpangenome.sh` runs `cactus-pangenome` in a Singularity container.

#### Required inputs

- `cactus_v3.1.4.sif`: Cactus container image
- `02_Config/vitis_seqfile.txt`: Cactus sequence configuration file
- Genome assemblies referenced by the sequence configuration file
- A reference genome named `V148` in the sequence configuration

The script uses:

- `03_JobStore/` as the persistent Toil job store;
- node-local storage under `/tmp` for temporary files;
- automatic `--restart` when `03_JobStore/` already exists;
- `V148` as the graph reference;
- 64 maximum CPU cores and 480 GB maximum memory.

Run the script inside an appropriate SLURM allocation:

```bash
bash 01_graphpangenome.sh
```

Principal outputs are written to:

```text
Output/
├── VitisPan.gbz
├── VitisPan.*
└── additional graph, index, and VCF outputs
```

The workflow requests clipped GBZ, ODGI, chromosome-level VG/OG, Giraffe, and VCF outputs. Interrupted runs can be resumed from `03_JobStore/`.

### 2. Pan-transcriptome construction

`02_pantranscriptome.sh` projects the reference annotation onto the graph and builds indexes for downstream graph-based RNA-seq mapping.

#### Required inputs

- `VitisPan.gbz`
- `Vv.graphnames.gff3`

The first column of the GFF3 file must use graph-compatible path names, for example:

```text
V148#0#Chr01
```

Before running, edit the `WORKDIR`, `GRAPH`, and `GFF` variables if necessary:

```bash
bash 02_pantranscriptome.sh
```

Principal outputs:

```text
RNAseq_graph/
├── graph.refpaths.txt
├── pantranscript.pg
├── pantranscript.xg
├── pantranscript.gcsa
├── pantranscript.gcsa.lcp
├── pantranscript.dist
├── pantranscript.transcripts.gbwt
├── pantranscript.info.tsv
└── pantranscript.files.txt
```

The XG, GCSA, distance, and transcript GBWT indexes can be used in subsequent VG RNA-seq mapping workflows such as `vg mpmap`.

### 3. Reference-based RNA-seq workflow

`03_rnaseq.sh` performs:

1. GFF3-to-GTF conversion and annotation validation;
2. extraction of known splice sites;
3. HISAT2 index construction;
4. paired-end read quality control with fastp;
5. splice-aware alignment with HISAT2;
6. BAM sorting, indexing, and validation;
7. gene-level fragment counting with featureCounts;
8. calculation of count, TPM, and FPKM matrices;
9. PCA, sample-correlation analysis, and MultiQC reporting.

#### Required project structure

```text
project/
├── Vv.fa
├── Vv.gff3
├── rawdata/
│   ├── GroupA-1_R1.fq.gz
│   ├── GroupA-1_R2.fq.gz
│   ├── GroupA-2_R1.fq.gz
│   └── GroupA-2_R2.fq.gz
└── 03_rnaseq.sh
```

Read files must follow:

```text
<sample>_R1.fq.gz
<sample>_R2.fq.gz
```

For automatic group and replicate parsing, the recommended sample name is:

```text
<group>-<replicate>
```

Run:

```bash
PROJECT_DIR=/path/to/project FC_STRAND=0 bash 03_rnaseq.sh
```

`FC_STRAND` accepts:

| Value | Library type |
| --- | --- |
| `0` | Unstranded |
| `1` | Stranded |
| `2` | Reversely stranded |

Principal expression outputs:

```text
06.expression/
├── gene_counts_matrix.tsv
├── gene_TPM_matrix.tsv
├── gene_FPKM_matrix.tsv
└── sample_expression_statistics.tsv
```

Quality-control outputs are written to `05.logs/` and `07.qc/`. Existing intermediate files are reused only after integrity and timestamp checks.

### 4. DESeq2 differential expression analysis

`04_DESeq_analysis.r` is designed for a two-group comparison.

#### Required inputs

- Raw gene-count matrix with gene IDs in the first column
- FPKM or TPM matrix with gene IDs in the first column
- Two-column, tab-delimited group file with a header
- Reference/control group name

Example group file:

```text
ID	Group
Vv-1	Vv
Vv-2	Vv
Vh-1	Vh
Vh-2	Vh
```

Run:

```bash
Rscript 04_DESeq_analysis.r \
  -i 06.expression/gene_counts_matrix.tsv \
  -g comparison_groups.tsv \
  -k 06.expression/gene_TPM_matrix.tsv \
  -r Vv \
  -f 0.05 \
  -c 2 \
  -o DESeq2_results \
  -p Vh_vs_Vv
```

Parameters:

| Option | Description |
| --- | --- |
| `-i`, `--input` | Raw count matrix |
| `-g`, `--group` | Sample group file |
| `-k`, `--fpkm` | FPKM/TPM matrix |
| `-r`, `--ref` | Reference group |
| `-f`, `--fdr` | Adjusted P-value threshold; default `0.05` |
| `-c`, `--fc` | Fold-change threshold; default `2` |
| `-o`, `--outdir` | Output directory |
| `-p`, `--prefix` | Output prefix |

Outputs include:

- complete differential expression results;
- significant DEG table;
- volcano plot;
- MA plot;
- group-level expression correlation plot.

### 5. OrthoFinder analysis

Place one protein FASTA file per genome in the `data/` directory:

```text
data/
├── Genome1.fa
├── Genome2.fa
└── Genome3.fa
```

Submit the SLURM job:

```bash
sbatch 05_orthoFinder.sh
```

The workflow uses:

- DIAMOND for sequence-similarity searches;
- MAFFT for multiple sequence alignment;
- IQ-TREE 3 for tree inference;
- 64 search threads and 32 analysis threads.

Results are written to `results/`.

### 6. Orthogroup pan-genome classification

`06_classify_orthogroups.py` uses two standard OrthoFinder output files:

- `Orthogroups.GeneCount.tsv`
- `Orthogroups.tsv`

Run:

```bash
python 06_classify_orthogroups.py \
  --count results/Orthogroups/Orthogroups.GeneCount.tsv \
  --tsv results/Orthogroups/Orthogroups.tsv \
  --output grape_pan
```

Classification rules:

| Category | Definition |
| --- | --- |
| Core | Present in all genomes |
| Soft-core | Present in more than 90% but not all genomes |
| Shell | Present in at least 10% and no more than 90% of genomes |
| Cloud | Present in fewer than 10% of genomes |

Outputs:

```text
grape_pan.orthogroup_category.tsv
grape_pan.summary.tsv
grape_pan.detail.tsv
```

### 7. Gene-structure comparison

`07_plot_genestructure.r` compares gene, CDS, exon, and intron lengths among core, soft-core, shell, and cloud genes.

#### Required inputs

`gene_info.txt`, with a header containing at least:

```text
geneID	GeneLen	CDSLen	ExonLen	IntronLen
```

`gene_family.txt`, without a header:

```text
Family	geneID
```

`family.freq_class.txt`, without a header:

```text
Family	Size	Class
```

Run:

```bash
Rscript 07_plot_genestructure.r
```

The script applies `log10(length + 1)` transformation and performs Wilcoxon rank-sum tests for:

- Core versus Soft-core;
- Core versus Shell;
- Core versus Cloud.

The four panels are saved as:

```text
Comparison_4panels.pdf
```

### 8. Ensemble feature selection and SHAP analysis

`08_SHAP_analysis.py` compares the Vv and Vh metabolomic profiles using weighted feature selection and an interpretable XGBoost model.

#### Required input

The script reads `data.txt` as a tab-delimited table:

- the first column contains compound or feature names;
- sample columns containing `Vv_` are assigned to class 0;
- sample columns containing `Vh_` are assigned to class 1;
- other annotation columns are ignored.

Run:

```bash
python 08_SHAP_analysis.py
```

The workflow:

1. removes zero-variance features;
2. calculates group means and `log2FC(Vh/Vv)`;
3. ranks features using a weighted consensus score;
4. retains features with `log2FC(Vh/Vv) > 1`;
5. selects up to 50 top-ranked features;
6. trains a regularized XGBoost classifier;
7. calculates SHAP values and SHAP interaction values;
8. generates global, feature-level, sample-level, and pairwise interaction visualizations.

The consensus score is calculated as:

```text
0.75 × ANOVA F-score
+ 0.05 × logistic regression importance
+ 0.05 × linear SVM importance
+ 0.05 × random forest importance
+ 0.05 × XGBoost importance
+ 0.05 × LightGBM importance
```

Principal outputs are written to `ML_Results/`:

```text
ML_Results/
├── MultiModel_Feature_Ranking_ALL.csv
├── MultiModel_Feature_Ranking_Vh_Upregulated.csv
├── Fig1_Classification_Performance.*
├── Fig2_Global_Contribution.*
├── Fig3_Dependence_Plots/
├── Fig4_Main_vs_Interaction.*
├── Fig5_Interaction_Matrix.*
├── Fig7_Interaction_Network.*
├── Fig8_SHAP_Heatmap.*
└── Fig10_2D_PDP/
```

Figures are exported as 300-dpi PNG and editable PDF files. The random seed is fixed at `42`.

> **Important:** the current script fits and explains the XGBoost model using all supplied samples. The classification-probability plot therefore describes the fitted dataset and should not be interpreted as independent test-set performance.

## Reproducibility notes

- Review and update hard-coded input paths before running each script.
- Keep sample names identical across count, expression, and group files.
- Use raw integer counts, not TPM or FPKM values, as DESeq2 input.
- Confirm RNA-seq library strandedness before setting `FC_STRAND`.
- Preserve the Cactus `03_JobStore/` directory when resuming an interrupted graph-pangenome run.
- Record software versions, parameters, reference files, and random seeds with each analysis.
- Large genome graphs, indexes, sequencing reads, and intermediate BAM files are not intended for Git version control.

## Citation

If you use these scripts, please cite the associated manuscript:

> *Multi-omics and explainable machine learning identify key enzymatic regulators of flavonoid biosynthesis in grape seeds.*

Author list, journal information, and DOI will be added after publication.

## Contact

For questions about the workflows or input formats, please open an issue in this repository.
