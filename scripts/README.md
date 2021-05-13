## Scripts 
All of the code used in this work can be found in this directory as MATLAB (`.m`) files, R (`.R`) files, or Jupyter notebooks (`.ipynb`). Moreover, generalised functions have been saved in the `./functions` sub-directory and `.py` scripts used to record accelerometry data from the bedside are available in the `./accel_recording_scripts` sub-directory.

### 1. [Complete Motion Feature Extraction](01_motion_feature_extraction.m)
In this `.m` script, we iterate through compiled triaxial accelerometry information from each patient, filter each axis with a high-pass (f_c = 0.2 Hz) 4th-order Butterworth filter, and extract 7 different motion features from non-overlapping 5 second windows. Outputs are saved as `.csv` feature tables. We also plot short examples of the accelerometry processing pipeline for Figure 1 in this script.

### 2. [Multiple Imputation of Missing Accelerometery Values](02_missing_feature_imputation.R)
In this `.R` script, we apply our multiple missing feature imputation algorithm. In the event of totally missing recordings (_n_ = 10/483), we impute upper extremity recordings with linear regression from ipsilateral upper extremity sensors, we impute lower extremity recordings with linear regression from contralateral upper extremity sensors, and bed sensors are imputed with sampling with replacement from the total distribution of bed sensor values. Then, the large majority of missing values were imputed with multiple, normalized time-series imputation with the `Amelia II` package. We create 9 imputations, each stored in a separate `.csv` file.

### 3. [Bed Motion Correction and Collection of Multiple Imputations](03_bed_movement_correction.R)
In this `.R` script, we correct gross-external movements by properly adjusting for the motion features calculated from the sensor placed at the foot of each patient's bed. Based on a literature-sourced threshold of SMA for human dynamic activity, we define distributions for each feature correspodning to static activity, and correct bed sensoer feature values from extremity sensors accordingly. As a result, we have 9 bed-corrected, imputed feature sets, each stored in a separate `.csv` file.

### 4. [Resampling of GCS data for classification](04_create_repeated_cv_folds.R)
In this `.R` file, we create repeated cross-validation (5 repeats of 5-fold CV) for each tested observation window based on the available GCS observations for each observation window. We principally use the `createMultiFolds` function from the `caret` package to stratify folds by outcome labels. Folds for each observation window are stored in a `.csv` file in a newly created directory.

### 5. [LOL embedding for dimensionality reduction](05_dim_reduction.R)
In this `.R` file, we train Linear Optimal Low Rank Projections [(LOL)](https://neurodata.io/lol/) on model training sets and reduce both training and validation sets to low-dimensional spaces prior to model training. Prior to LOL, we normalize feature spaces per the distribution of feature type and sensor combinations. This enables us to use LOL coefficients to compare feature type and sensor significance. 

### 6. [Train and evaluate prediction models and measure feature significance scores](06_prediction_models.R)
In this `.R` file, we train and validate logistic regression models (GLM) for threshold-level GCSm detection, threshold-level GOSE at discharge prediction, and threshold-level GOSE at 12 months prediction. 

### 7. [Calculate bootstrapping bias-corrected cross-validation (BBC-CV) area under the receiver operating characteristic curve (AUC) and classification metrics](07_calculate_metrics.ipynb)

### 8. [Calculate model calibration on validation set predictions](08_model_calibration_calculation.R)

### 9. [Case study exploration](09_case_study_exploration.R)

### 10. [Create manuscript tables and perform statistical analyses](10_tables_and_statistics.R)


### 11. [Plot figures for the manuscript](11_manuscript_figures.R)
