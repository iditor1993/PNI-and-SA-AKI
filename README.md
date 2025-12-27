Prognostic Nutritional Index (PNI) and Sepsis-Associated Acute Kidney Injury: A MIMIC-IV Cohort Study
Description
This repository contains the complete R analysis code and dataset for a retrospective cohort study investigating the association between Prognostic Nutritional Index (PNI) and clinical outcomes in patients with sepsis-associated acute kidney injury (AKI) from the MIMIC-IV database.
Study Objective
To evaluate PNI as a prognostic biomarker for 30-day mortality and to assess its incremental predictive value when combined with established critical illness severity scores in septic AKI patients.
Methods Overview
The analysis pipeline includes:
 
Data Processing: MIMIC-IV extraction, variable selection, missing data handling (multiple imputation), and outlier removal
 
PNI Calculation: Computed as  10 × serum albumin (g/dL) + 0.005 × absolute lymphocyte count (per μL) 
 
Statistical Analysis:
 
Baseline characteristic comparison across PNI quartiles
 
Kaplan-Meier survival analysis with log-rank test
 
Univariable and multivariable Cox proportional hazards regression
 
Restricted cubic spline (RCS) analysis for non-linear dose-response relationships
 
Threshold effect analysis using two-piecewise Cox models
 
Subgroup interaction analysis (forest plots)
 
C-index comparisons for model discrimination
Repository Structure
 
 All code for DYW project.R : Complete analysis script covering data cleaning, statistical modeling, and visualization
 
 mimiciv_export.xlsx : Primary dataset containing clinical variables, laboratory values, severity scores, and outcomes
 
Generated outputs include survival curves, forest plots, RCS plots, and C-index comparison charts
Key Variables
 
Primary Exposure: Prognostic Nutritional Index (PNI)
 
Primary Outcome: 30-day all-cause mortality
 
Covariates: Demographics, comorbidities (Charlson index), severity scores (SOFA, APS III, SAPS II), organ support (MV, CRRT, vasoactive agents), and laboratory parameters
Requirements
 
R packages:  tidyverse ,  survival ,  survminer ,  rms ,  ggrcs ,  mice ,  tableone ,  forestploter ,  pROC 
Citation
If you use this code or data, please cite the original study and acknowledge the MIMIC-IV database.
Usage Notes: The script is fully reproducible with all random seeds set (123) and includes detailed annotations for each analytical step.
