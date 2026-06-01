# Hepatitis C Virus Can Induce Gene Expression Changes Associated with Hepatocarcinogenesis

This repository contains analysis code and supporting result files for the paper:

**Hepatitis C Virus Can Induce Gene Expression Changes Associated with Hepatocarcinogenesis**
Tuyana Boldanova, Fredrik Trulsson, Fahim Ebrahimi, Andrej Benjak, et al.
*JHEP Reports*, 2026.
DOI: [10.1016/j.jhepr.2026.101897](https://doi.org/10.1016/j.jhepr.2026.101897)

## Overview

The analyses in this repository investigate whether hepatitis C virus (HCV) infection can induce persistent transcriptomic changes in chronic hepatitis C (CHC) liver tissue that are associated with hepatocellular carcinoma (HCC) development.

The repository includes notebooks and scripts used for differential expression analysis, comparison of baseline and follow-up transcriptomic profiles, identification of persistent HCV-associated genes, and evaluation of transcriptomic predictors of HCC-related phenotypes.

## Repository contents

* `Baseline_vs_follow_ups.ipynb`
  Analysis comparing baseline and follow-up samples.

* `pyDEseq_transcriptome_baseline.ipynb`
  Differential expression analysis of baseline transcriptomic profiles.

* `predict_hcc_HCV_genes_20112025.ipynb`
  Analysis of HCV-associated genes in relation to HCC-associated expression changes.

* `predict_hcc_elastic_net.ipynb`
  Elastic-net modelling for HCC-related transcriptomic prediction.

* `predict_hcc_linear_regression.ipynb`
  Linear-regression-based analysis of transcriptomic predictors.

* `imsig_analysis.R`
  Immune signature analysis.

* `base_hcv_induced_persistant_22122025.csv`
  Processed table of persistent HCV-induced genes.

* `elasticnet_parallel_bootstrap_results_22122025.csv`
  Processed elastic-net bootstrap results.

* `imsig_input_metadata.csv`
  Metadata used for immune signature analysis.

## Usage

The notebooks are provided as an archival record of the analyses performed for the manuscript. They can be opened and run using Jupyter Notebook or JupyterLab.

Required Python and R package versions may depend on the local analysis environment. Users wishing to reproduce the analyses should inspect the import statements and package loading commands in each notebook or script.

## Data availability

This repository contains code and processed analysis outputs associated with the manuscript. Raw sequencing data and additional study data are described in the data availability statement of the paper.

## Citation

If you use this code, please cite the associated paper:

Boldanova T, Trulsson F, Ebrahimi F, Benjak A, et al. **Hepatitis C Virus Can Induce Gene Expression Changes Associated with Hepatocarcinogenesis.** *JHEP Reports*. 2026. https://doi.org/10.1016/j.jhepr.2026.101897

## License

This repository is released under the MIT License. See [`LICENSE`](LICENSE) for details.
