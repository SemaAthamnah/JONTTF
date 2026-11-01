# JONTTF — Signed Brain Network Community Detection

**Author:** Sema Athamnah  
**Email:** athamnah@msu.edu  
**Affiliation:** Department of Biomedical Engineering, Michigan State University

## Method overview

![JONTTF overview](images/JONTTF.jpg)

*Overview of JONTTF for signed multiplex networks, jointly modeling positive and negative edges through a shared node membership tensor.*

## Contents
- `main.m`: Run one SSBM sample from one noise setting and perform a JONTTF hyperparameter grid search.
- `funs/`: JONTTF solver, grid search, clustering metrics, and required numerical helpers.
- `data/`: Place the SSBM `.mat` files here (not included yet).

## Usage
1. Place `NTF_SSBM_NB_1s.mat` and/or `NTF_SSBM_NW_1s.mat` in `data/`.
2. Open MATLAB in this directory and run `main`.
3. Change `dataset = 'NB'` to `'NW'` to select the other dataset.
4. Adjust parameter grids in `main.m` if needed.

The script selects the first available `SSBM_i`/`GT_i` setting and first sample, infers `k` from the ground-truth labels, and reports clustering metrics. The ground truth is used for evaluation and to infer `k`, not for factor fitting. This is an illustrative benchmark run, not a full multi-sample evaluation. No output files are written by default.

**Dependencies:** MATLAB with the functions in `funs/`. The provided clustering-measure implementation may require additional MATLAB functionality depending on your installation.
