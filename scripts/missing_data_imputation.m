%% Authors: Shubhayu Bhattacharyay, B.S. Candidate, Matthew Wang, B.S. Candidate
% Department of Biomedical Engineering
% Department of Applied Mathematics and Statistics
% Whiting School of Engineering, Johns Hopkins University
% email address: shubhayu@jhu.edu
% March 2019; Last revision: 20-Feb-2020
%% ------------- BEGIN CODE --------------
%Set Directory and Procure all feature data
addpath('functions/')
tic
cd ..
cd motion_feature_data
load('complete_sensor_data.mat');
toc

dim_of_sensors=size(sensors);
sensor_count = dim_of_sensors(1);
feature_count = dim_of_sensors(2);

studyPatientsPY = [2,3,4,5,	6,	7,	8,	9,	10,	11,	12,	13, ...
    14,	15,	16,	17,	18,	19,	20,	21,	22,	23,	24,	26,	27,	28,	29,	30, ....
    31,	32,	33,	34,	35,	36,	37,	38,	39,	40,	41,	42,	43,	49,	51,	46, .....
    47,	48,	50,	52,	53,	54,	55,	56,	57,	59,	60,	61,	62,	63,	64,	65, ......
    67,	68];

[sortedPY,sort_order] = sort(studyPatientsPY);

% Get Time Point Labels (in datenum format)
cd band_power
load('band_power02.mat')
t = bandPower{2, 1};
t_copy = t;
clear bandPower
cd ..
cd ..
cd scripts/
%% Identify rows that are completely missing data (and replace with NaN)
%contains indices of the totally missing time series as follows:
%Row_1 is the sensor number (1 through 7)
%Row_2 is the feature_count (1 through 7) that corresponds to "feature_names"
%Row_3 is the patient number that corresponds to the patient index (1
%through 62)
[totallyMissingIdxs, sensors]=get_totallyMissingIdxs(sensor_count,...
    feature_count,sensors);

%% Characterize percentage of missing data per time-series recording
tic
[missing_percentages,missing_time_series,missingIdxs] = ...
    characterize_missing_data(sensors,t,studyPatientsPY,true);
toc
%% Identifying distributions of the features:
% We first wish to identify whether the feature datapoints have any
% association with time of day. We begin my calculating the mean and
% variance of the feature columns by ignoring NaN values.

TOD_Means = cellfun(@(X) nanmean(X,1),sensors,'UniformOutput',false);
TOD_Std = cellfun(@(X) nanstd(X,1),sensors,'UniformOutput',false);
TOD_Diff = cellfun(@(X) abs(diff(X)),TOD_Means,'UniformOutput',false);

% Mean by time of day figures:
plot_TOD_Means(t,TOD_Means)

% STD by time of day figures:
plot_TOD_Std(t,TOD_Std)

% Diff Plots by time of day
plot_TOD_Diff(t,TOD_Diff)

% Histogram of different sensor/feature combination means by time of day:
plot_TOD_Hist(TOD_Means)

%% Finding no motion thresholds for each feature based on SMA:

%NOTE: section takes about 20 seconds to run

SMA_threshold=0.1;
feature_thresholds=find_noMotion_thresholds(SMA_threshold,sensors,...
    feature_names);

%% Fit Zero-inflated poisson distribution to all feature spaces

% NOTE: takes about 4.5 seconds to run

cd ..
cd clinical_data/
load('clinical_extraction_output.mat')
cd ..
cd scripts/

clearvars dc_dataset yr_dataset dc_dataset_labels ...
    yr_dataset_labels boxcox_lambdas zscore_mus zscore_sigs ....
    imputeDataset

patient_table = fit_zip(sensors,patient_table,feature_thresholds,...
    sort_order,feature_names);
%% Regression Strategy:
% Takes about 3-4 seconds to run
tic
% NOTES:
% - for Wrist and Elbow, we use the other available arm sensor as the best
% predictor

% wrist indices (default 4 and 7):
wrist_Idx = [4,7];
% elbow indices (defult 3 and 6):
elbow_Idx = [3,6];

patient_table=wrist_elbow_Regression(patient_table,wrist_Idx,elbow_Idx,...
    totallyMissingIdxs,feature_names);

% - for ankle, it seems like none of the available predictors pass the
% predictor t-test, so we will simply use wrist, elbow, and interaction
% terms of the available side

ankle_Idx = [2,5];

patient_table = ankle_Regression(patient_table,ankle_Idx,...
    totallyMissingIdxs,feature_names);

% - for bed imputation, I will simply take the average of all the available
% data:
bed_Idx = 1;

patient_table = bed_imputation(patient_table,bed_Idx,totallyMissingIdxs,...
    feature_names);

% - We should use pis and lambdas as predictive features themselves!
toc
%% Totally missing data imputation

sensors=impute_totallyMissingData(sensors,patient_table,...
    totallyMissingIdxs,feature_names,feature_thresholds,....
    studyPatientsPY,sortedPY,t);

%% Quasi-missing data imputation
% for recordings with more than (50%) missing data, we impute by a similar
% manner as the totally missing data:
quasi_threshold=0.5;

sensors = impute_quasiMissingData(sensors,patient_table,quasi_threshold,...
    missing_percentages,missingIdxs,feature_names,feature_thresholds,....
    studyPatientsPY,sortedPY);

%% Remaining missing data imputation

% Warning: takes about ~15 minutes to run

trainWindow=60; % in minutes

sensors = impute_rmngMissingData(sensors,quasi_threshold,trainWindow,...
    missing_percentages,missingIdxs,feature_names,feature_thresholds,t);

cd ..
cd motion_feature_data/

save('imputed_complete_sensor_data.mat','sensors','feature_names');

cd ..
cd scripts/