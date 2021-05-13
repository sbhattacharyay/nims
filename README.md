# Assessing the clinical utility of triaxial accelerometry in severe brain injury monitoring and prognosis
A pilot study of the Neurological Injury Motion Sensing (NIMS) project

## Contents

- [Overview](#overview)
- [Abstract](#abstract)
- [Code](#code)
- [License](./LICENSE)
- [Citation](#citation)

## Overview

This repository contains the code underlying the article entitled **Assessing the clinical utility of triaxial accelerometry in severe brain injury monitoring and prognosis** from the Johns Hopkins University [Laboratory of Computational Intensive Care Medicine](https://lcicm.jhmi.edu/). In this file, we present the abstract, to outline the motivation for the work and the findings, and then a brief description of the code with which we generate these finding and achieve this objective.\
\
The code on this repository is commented throughout to provide a description of each step alongside the code which achieves it.

## Abstract

## Code 
All of the code used in this work can be found in the `./scripts` directory as MATLAB (`.m`) files, R (`.R`) files, or Jupyter notebooks (`.ipynb`). Moreover, generalised functions have been saved in the `./scripts/functions` sub-directory and `.py` scripts used to record accelerometry data from the bedside are available in the `./scripts/accel_recording_scripts` sub-directory.

### 1. [Complete Motion Feature Extraction](scripts/01_motion_feature_extraction.m)
In this `.m` script, we iterate through compiled triaxial accelerometry information from each patient, filter each axis with a high-pass (_f_c_ = 0.2 Hz) 4th-order Butterworth filter, and extract 7 different motion features from non-overlapping 5 second windows. Outputs are saved as `.csv` feature tables. We also plot short examples of the accelerometry processing pipeline for Figure 1 in this script.

### 2. [Multiple Imputation of Missing Accelerometery Values](scripts/02_missing_feature_imputation.R)
In this `.R` script, we apply our multiple missing feature imputation algorithm. In the event of totally missing recordings (_n_ = 10/483), we impute upper extremity recordings with linear regression from ipsilateral upper extremity sensors, we impute lower extremity recordings with linear regression from contralateral upper extremity sensors, and bed sensors are imputed with sampling with replacement from the total distribution of bed sensor values. Then, the large majority of missing values were imputed with multiple, normalized time-series imputation with the `Amelia II` package. We create 9 imputations, each stored in a separate `.csv` file.

### 3. [Bed Motion Correction and Collection of Multiple Imputations](scripts/03_bed_movement_correction.R)
In this `.R` script, we correct gross-external movements by properly adjusting for the motion features calculated from the sensor placed at the foot of each patient's bed. Based on a literature-sourced threshold of SMA for human dynamic activity, we define distributions for each feature correspodning to static activity, and correct bed sensoer feature values from extremity sensors accordingly. As a result, we have 9 bed-corrected, imputed feature sets, each stored in a separate `.csv` file.

### 4. [Resampling of GCS data for classification](scripts/04_create_repeated_cv_folds.R)
In this `.R` file, we create repeated cross-validation (5 repeats of 5-fold CV) for each tested observation window based on the available GCS observations for each observation window. We principally use the `createMultiFolds` function from the `caret` package to stratify folds by outcome labels. Folds for each observation window are stored in a `.csv` file in a newly created directory.

### 5. [LOL embedding for dimensionality reduction](scripts/05_dim_reduction.R)
In this `.R` file, we train Linear Optimal Low Rank Projections [(LOL)](https://neurodata.io/lol/) on model training sets and reduce both training and validation sets to low-dimensional spaces prior to model training. Prior to LOL, we normalize feature spaces per the distribution of feature type and sensor combinations. This enables us to use LOL coefficients to compare feature type and sensor significance. 

### 6. [Train and evaluate prediction models and measure feature significance scores](scripts/06_prediction_models.R)
In this `.R` file, we train and validate logistic regression models (GLM) for threshold-level GCSm detection, threshold-level GOSE at discharge prediction, and threshold-level GOSE at 12 months prediction. We train and evaluate models of varying observation windows and target dimensionalities. In this script, we also calculate our feature significance scores. This is equiavalent to the absolute LOL coefficient weighted by the trained linear coefficients of the corresponding logistic regression model.

### 7. [Calculate bootstrapping bias-corrected cross-validation (BBC-CV) area under the receiver operating characteristic curve (AUC) and classification metrics](scripts/07_calculate_metrics.ipynb)
In this `.ipynb` notebook, we calculate AUCs, ROC curves, and classification metrics for each observation window based on the validation set predictions returned by our models. We use bias-corrected bootstrapping for repeated cross-validation [(Repeated BBC-CV)](https://doi.org/10.1007/s10994-018-5714-4) to calculate 95% confidence intervals for the metrics and the ROC curve. This script is programmed to perform bootstrapping in parallel on 10 cores.

### 8. [Calculate model calibration on validation set predictions](scripts/08_model_calibration_calculation.R)
In this `.R` file, we calculate probability calibration curves and associated calibration metrics for each observation window based on the validation set predictions returned by our models. We use bias-corrected bootstrapping for repeated cross-validation [(Repeated BBC-CV)](https://doi.org/10.1007/s10994-018-5714-4) to calculate 95% confidence intervals for the metrics and the calibration curve. This script is programmed to perform bootstrapping in parallel on 10 cores.

### 9. [Case study exploration](scripts/09_case_study_exploration.R)
In this `.R` file, we retrospectively examine predictions of _Pr_(GCSm > 4) in patients (_n_ = 6) who experienced neurological transitions across this threshold to visually determine potential clinical utility of the accelerometry-based system. For each of the 6 patients, we train optimally discriminating detection models (one with a shorter observation window of 27 minutes and one with a longer observation window of 6 hours) on the remaining patient set and validate predictions on the case study patients specifically over a large, continuously overlapping observation window set. We bootstrap across imputations to produce 95% confidence intervals that account for variation due to imputation on the predictions. Then, we prepare the probability trajectories for plotting.

### 10. [Create manuscript tables and perform statistical analyses](scripts/10_tables_and_statistics.R)
In this `.R` file, we construct manuscript tables and perform miscellaneous statistical analyses for different parts (including figures) of the manuscript and supplementary materials. In addition to the classification metrics calculated in script no. 7, we also calculate classification accuracy with repeated BBC-CV in this script.

### 11. [Plot figures for the manuscript](scripts/11_manuscript_figures.R)
In this `.R` file, we produce the figures for the manuscript and the supplementary figures. The large majority of the quantitative figures in the manuscript are produced using the `ggplot` package.

## Citation
