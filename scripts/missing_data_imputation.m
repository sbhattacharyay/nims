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
%% Get Time Point Labels (in datenum format)
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

% In order to impute the missing data, we first need to understand how
% "no motion" manifests in each of the given feature types.

% We assume that SMA is our ground truth in calculating No Motion
% thresholds.
SMA_threshold = 0.1;
SMA_nm_percentages = cellfun(@(x) noMotion_th(x,SMA_threshold),...
    sensors(:,6),'UniformOutput',false);

% Warning: Elapsed Time is 19.35 seconds

tic
% Based on the SMA_nm_percentages, we will derive thresholds that allow us
% to most closely match the percentages from SMA:

sma_idx = find(feature_names == "sma");
mf_idx = find(feature_names == "med_freq");

feature_thresholds = NaN(dim_of_sensors(2),1);
feature_thresholds(sma_idx) = SMA_threshold;

% Finding Band Power threshold:
j=1;
temp_feature_data = sensors(:,j);
min_norm = realmax;
opt_k = NaN;

for k = 0:0.001:0.15
    feat_percentages = cellfun(@(x) noMotion_th(x,k),...
        temp_feature_data,'UniformOutput',false);
    curr_norm=norm_calc(feat_percentages,SMA_nm_percentages);
    if curr_norm < min_norm
        opt_k = k;
        min_norm = curr_norm;
    end
end
feature_thresholds(j) = opt_k;

% Finding Freq Entropy threshold:
j=2;
temp_feature_data = sensors(:,j);
min_norm = realmax;
opt_k = NaN;

for k = 1.66:0.001:1.7
    feat_percentages = cellfun(@(x) noMotionFE_th(x,k),...
        temp_feature_data,'UniformOutput',false);
    curr_norm=norm_calc(feat_percentages,SMA_nm_percentages);
    if curr_norm < min_norm
        opt_k = k;
        min_norm = curr_norm;
    end
end
feature_thresholds(j) = opt_k;

% Finding Freq Pairs 1 threshold:
j=3;
temp_feature_data = sensors(:,j);
min_norm = realmax;
opt_k = NaN;

for k = 0.003:0.001:0.02
    feat_percentages = cellfun(@(x) noMotion_th(x,k),...
        temp_feature_data,'UniformOutput',false);
    curr_norm=norm_calc(feat_percentages,SMA_nm_percentages);
    if curr_norm < min_norm
        opt_k = k;
        min_norm = curr_norm;
    end
end
feature_thresholds(j) = opt_k;

% Finding Freq Pairs 2 threshold:
j=4;
temp_feature_data = sensors(:,j);
min_norm = realmax;
opt_k = NaN;

for k = 0.001:0.0001:0.006
    feat_percentages = cellfun(@(x) noMotion_th(x,k),...
        temp_feature_data,'UniformOutput',false);
    curr_norm=norm_calc(feat_percentages,SMA_nm_percentages);
    if curr_norm < min_norm
        opt_k = k;
        min_norm = curr_norm;
    end
end
feature_thresholds(j) = opt_k;

% Finding Med Freq threshold:
j=5;
temp_feature_data = sensors(:,j);
min_norm = realmax;
opt_k = NaN;

for k = 2.5:0.01:3.5
    feat_percentages = cellfun(@(x) noMotionFE_th(x,k),...
        temp_feature_data,'UniformOutput',false);
    curr_norm=norm_calc(feat_percentages,SMA_nm_percentages);
    if curr_norm < min_norm
        opt_k = k;
        min_norm = curr_norm;
    end
end
feature_thresholds(j) = opt_k;

% Finding Wavelet threshold:
j=7;
temp_feature_data = sensors(:,j);
min_norm = realmax;
opt_k = NaN;

for k = 1:0.1:15
    feat_percentages = cellfun(@(x) noMotion_th(x,k),...
        temp_feature_data,'UniformOutput',false);
    curr_norm=norm_calc(feat_percentages,SMA_nm_percentages);
    if curr_norm < min_norm
        opt_k = k;
        min_norm = curr_norm;
    end
end
feature_thresholds(j) = opt_k;
toc
%% Impute the totally-missing rows
%First, given our newly derived thresholds, we calculate the "no motion"
%for each of our patients per sensor per feature:

nm_percentages_for_features = {};

for j = 1:dim_of_sensors(2)
    temp_feature_data = sensors(:,j);
    k = feature_thresholds(j);
    if j == 2 || j == 5
        nm_percentages_for_features{j,1}=cellfun(@(x) noMotionFE_th(x,k), ...
            temp_feature_data,'UniformOutput',false);
    else
        nm_percentages_for_features{j,1}=cellfun(@(x) noMotion_th(x,k), ...
            temp_feature_data,'UniformOutput',false);
    end
end


% Then, we cycle through the completely missing time series
dim_TM = size(totallyMissingIdxs);

bp_nm_range = [0 feature_thresholds(1)];
fe_nm_range = [0 feature_thresholds(6)];
sma_nm_range = [0 feature_thresholds(6)];
sma_nm_range = [0 feature_thresholds(6)];
sma_nm_range=[0 feature_thresholds(6)];
sma_nm_range = [0 feature_thresholds(6)];
wav_nm_range = [0 feature_thresholds(7)];

for i = 1:dim_TM(2)
    sensIdx = totallyMissingIdxs(1,i);
    featIdx = totallyMissingIdxs(2,i);
    ptIdx = totallyMissingIdxs(3,i);
    draws=rand(1,length(t));
    curr_perc=nm_percentages_for_features{featIdx,1}{sensIdx,1};
    NM_imputes=draws<curr_perc(ptIdx);
    curr_mat=sensors{sensIdx,featIdx};
    curr_mat(ptIdx,NM_imputes)= 0;
end
