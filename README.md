# OUD-survival

This repository contains **example code** for preparing time-interval data and running survival analysis for opioid use disorder (OUD) onset, remission, and relapse.

The scripts are designed for use on the **All of Us Researcher Workbench**. Because All of Us data are controlled-access, **no individual-level data are included in this repository**. Only generic, de-identified *example* scripts are provided.

---

## Overview

Our project:

> **Integrating Survival Models and a Therapy Knowledge Graph to Characterize Onset, Remission, and Relapse in Opioid Use Disorder**  
> *Mengman Wei\*, Stanislav Listopad, and Qian Peng\**

We use interval-based survival models to study OUD onset, remission, and relapse, and integrate the results with a therapy knowledge graph for drug repurposing and treatment characterization.

This repo focuses on the **methodological side**: how to  
1. Construct interval-split survival data from longitudinal EHR-like sources,  
2. Fit Cox / time-varying survival models on those intervals, and  
3. Build a knowledge graph for drug repurposing.

---

## What’s in this repository?

> **Note:** Adjust file paths/names below to match your actual repo structure.

- `prepare_intervals.py`  
  Example R script showing how to:
  - Read visit-level and condition/exposure tables  
  - Define OUD onset/remission/relapse events  
  - Create episode/interval-level data (start time, end time, event indicator)  
  - Add time-varying covariates (e.g., diagnoses, medications, encounters)

- `survival_analysis.R`  
  Example R script showing how to:
  - Load the prepared interval dataset  
  - Fit Cox models (including time-varying covariates) using the **survival** package  
  - Extract hazard ratios, confidence intervals, and model diagnostics  

- `knowledge_graph.py`  
  Example Python script for building a drug repurposing knowledge graph using publicly available gene information.

No real All of Us table IDs, SQL queries, or participant-level data are included.

---

## Working with All of Us data

The **All of Us Research Program** data are **controlled-access**. To use these scripts with real data, you must obtain Researcher Workbench access and complete all required training.

---

## Reference and citation

If you use these scripts or ideas in your own work, please cite:

> **Integrating Survival Models and a Therapy Knowledge Graph to Characterize Onset, Remission, and Relapse in Opioid Use Disorder**  
> Mengman Wei\*, Stanislav Listopad, and Qian Peng\*

You may also cite this repository directly (e.g., as a software / code resource) with its GitHub URL and commit/tag.


