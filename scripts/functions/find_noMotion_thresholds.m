function [feature_thresholds,noMotion_time_series] = ...
    find_noMotion_thresholds(SMA_threshold,sensors,feature_names)

% In order to impute the missing data, we first need to understand how
% "no motion" manifests in each of the given feature types.

% We assume that SMA is our ground truth in calculating No Motion
% thresholds.

dim_of_sensors=size(sensors);
[n,~] = size(sensors{1,1});

SMA_nm_percentages = cellfun(@(x) noMotion_th(x,SMA_threshold),...
    sensors(:,6),'UniformOutput',false);

stacked_matrix=[];

for i = 1:dim_of_sensors(1)
    curr_matrix = sensors{i,6};
    stacked_matrix = [stacked_matrix; curr_matrix];
end

noMotion_time_series=sum(stacked_matrix < SMA_threshold)/...
    (n*(dim_of_sensors(1)));

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
    feat_percentages = cellfun(@(x) noMotion_th_neg(x,k),...
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
    feat_percentages = cellfun(@(x) noMotion_th_neg(x,k),...
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

end

