#for drought transcriptomics (sample: root)

# Cacao Drought Transcriptomics Pipeline

This repository contains a Snakemake-based pipeline for analyzing transcriptomics data from cacao under drought stress. It includes:

- **DESeq2** for differential expression
- **GO enrichment** using clusterProfiler
- **GSEA** for pathway analysis
- **WGCNA** for module detection
workflow/
├── Snakefile
├── scripts/
│   ├── run_deseq2.R
│   ├── run_go_enrichment.R
│   ├── run_gsea.R
│   └── run_wgcna.R
├── envs/
│   ├── deseq2.yaml
│   ├── go_enrichment.yaml
│   ├── gsea.yaml
│   └── wgcna.yaml
├── config/
│   └── sample_metadata.csv
├── results/
data/
│   └── expression_matrix.csv

## Usage

```bash
bash run_pipeline.sh
docker build -t cacao_pipeline .
docker run --rm -v $(pwd):/pipeline cacao_pipeline

```bash
mv "Additional file S1. DESeq2_Top100_Genes_root.xlsx" data/root_DEGs.xlsx
git add data/root_DEGs.xlsx
git commit -m "Add root DEGs Excel file"
git push
