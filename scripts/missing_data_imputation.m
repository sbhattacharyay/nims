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
path = pwd;
directory = dir;
load('complete_sensor_data.mat');
toc

dim_of_sensors=size(sensors);
sensor_count = dim_of_sensors(1);
feature_count = dim_of_sensors(2);

studyPatientsPY = [2, 3,	4,	5,	6,	7,	8,	9,	10,	11,	12,	13, ...
    14,	15,	16,	17,	18,	19,	20,	21,	22,	23,	24,	26,	27,	28,	29,	30, ....
    31,	32,	33,	34,	35,	36,	37,	38,	39,	40,	41,	42,	43,	49,	51,	46, .....
    47,	48,	50,	52,	53,	54,	55,	56,	57,	59,	60,	61,	62,	63,	64,	65, ......
    67,	68];

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
SMA_threshold=0.1;
feature_thresholds=find_noMotion_thresholds(SMA_threshold,sensors,feature_names);

%% Impute the totally-missing rows
%First, given our newly derived thresholds, we calculate the "no motion"
%for each of our patients per sensor per feature:

%Then, fit the values above the threshold to a decaying exponential

tic
nm_percentages_for_features = {};
lambda_for_features = {};

for j = 1:dim_of_sensors(2)
    temp_feature_data = sensors(:,j);
    k = feature_thresholds(j);
    if j == 2 || j == 5
        nm_percentages_for_features{j,1}=cellfun(@(x) noMotion_th_neg(x,k), ...
            temp_feature_data,'UniformOutput',false);
        lambda_for_features{j,1}=cellfun(@(x) fit_exp(x,k), ...
            temp_feature_data,'UniformOutput',false);
    else
        nm_percentages_for_features{j,1}=cellfun(@(x) noMotion_th(x,k), ...
            temp_feature_data,'UniformOutput',false);
        lambda_for_features{j,1}=cellfun(@(x) fit_exp_neg(x,k), ...
            temp_feature_data,'UniformOutput',false);
    end
end
toc

%% Create a regression to parameters
% We will regress 

%% Cycle through totally missing recordings and impute

rng(1)

dim_TM = size(totallyMissingIdxs);

bp_nm_range = [0 feature_thresholds(1)];
fe_nm_range = [feature_thresholds(2) 1.707];
fp1_nm_range = [0 feature_thresholds(3)];
fp2_nm_range = [0 feature_thresholds(4)];
mf_nm_range = [feature_thresholds(5) 3.2];
sma_nm_range = [0 feature_thresholds(6)];
wav_nm_range = [0 feature_thresholds(7)];

for i = 1:dim_TM(2)
    sensIdx = totallyMissingIdxs(1,i);
    featIdx = totallyMissingIdxs(2,i);
    ptIdx = totallyMissingIdxs(3,i);
    
    draws=rand(1,length(t));
    
    if featIdx == 2 || featIdx == 5
    else 
    end
    
    curr_perc=nm_percentages_for_features{featIdx,1}{sensIdx,1};
    NM_imputes=draws<curr_perc(ptIdx);
    curr_mat=sensors{sensIdx,featIdx};
    curr_mat(ptIdx,NM_imputes)= 0;
end
