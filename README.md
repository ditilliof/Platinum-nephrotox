# Platinum-nephrotox – Analysis Workflow

This repository contains the complete workflow for processing data, performing key event (KE) mapping, calibrating qAOP and PK models in CmdStan, and performing QIVIVE and benchmark dose (BMD) analysis for the study *Quantitative Adverse Outcome Pathway Modeling of Cisplatin-Induced Nephrotoxicity: Developing In Vitro and In Vivo Models for Predictive Extrapolation*.

---

## 1. Data Analysis (normalisation)

### In vivo data

In `data_analysis/`, use the **BioSpyder_analysis_pipeline_splitReplicate.R** script to compute the log2FC dataset from the raw rat kidney expression data.

- **Input:** raw in vivo count/expression tables  
- **Output:** log2FC in vivo datasets ready for upload to TXG-MAPr  

### In vitro data (RPTEC/TERT1)

Still in `data_analysis/`, run:

- `tempOseq_count_to_DEG.Rmd`

This script computes the log2FC dataset from the raw Tempo-Seq counts.

- **Input:** raw RPTEC/TERT1 Tempo-Seq counts  
- **Output:** log2FC in vitro data ready for upload to TXG-MAPr  

---

## 2. Eigengene Scores via TXG-MAPr

Upload the normalised in vivo and in vitro results to **TXG-MAPr**.

TXG-MAPr produces module-level **eigengene score (EGS)** tables.

- **Input:** normalised expression or DEG matrices  
- **Output:** eigengene scores for each module (in vivo and in vitro)

Place these EGS tables in the expected data folders (`data/in vivo` or `data/in vitro`) so that KE mapping and model fitting scripts can use them.

---

## 3. Key Event Mapping (KE_mapping/)

In the `KE_mapping/` directory, run:

- `KEmappingIGSEGs.Rmd`

Transcriptomics-based KE mapping is performed using both IGS and EG scores from TXG-MAPr. 
Additionally, PI staining and histopathology data are integrated to map Cell Death (both in vivo and in vitro) and Kidney Failure (in vivo).

- **Input:** Log2FC data, TXG-MAPr EG scores, raw PI staining and histopathology data  
- **Output:** KE-mapped data used as inputs for qAOP model fitting  

---

## 4. Model Fitting (fitting_scripts/)

Use the scripts in `fitting_scripts/` to run all PK and qAOP model fits:

- `fitPK2.R` – in vivo PK model (2 compartments)  
- `fitPK3.R` – in vivo PK model (3 compartment)  
- `Vitro_fit.R` – in vitro qAOP model  
- `Vivo_fit.R` – in vivo qAOP model  

### **Required execution order for in vivo qAOP model**

1. Run `fitPK2.R`  
2. Run `fitPK3.R`  
3. Run `Vivo_fit.R`  

The in vivo qAOP model uses PK parameter posteriors from PK2 and PK3, so these fits **must** be executed first.  
The in vitro qAOP model (`Vitro_fit.R`) can be run independently.

All CmdStan outputs (CSV, RDS, diagnostics) are stored in`Stanresults/`.   

---

## 5. Posterior Analysis, QIVIVE, and BMD (Posterior_analysis/)

After model fitting, navigate to `Posterior_analysis/` to analyse posterior distributions, model predictions, plotting time-course model fits and QIVIVE results:

- `PKplot.Rmd` – posterior and predictive checks for PK models  
- `invitroplot_QIVIVE.Rmd` – in vitro qAOP results and QIVIVE  
- `invivoplot.Rmd` – in vivo qAOP posterior analysis
- `BMD.Rmd` – benchmark dose (BMD) calculations  

These notebooks:

- load Stan fit objects  
- generate diagnostic plots
- generate time-course plots comparing model fit to KE data for all conditions.
- perform qAOP-based QIVIVE 
- Performs BMD analysis  

---

## 6. Outputs (Plots&output/)

All **final** results and figures are saved in:

- `Plots&output/` 

Use this folder for inspection, export, and publication-ready material.

---


