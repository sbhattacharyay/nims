%% Authors: Shubhayu Bhattacharyay, B.S. Candidate, Matthew Wang, B.S. Candidate
% Department of Biomedical Engineering
% Department of Applied Mathematics and Statistics
% Whiting School of Engineering, Johns Hopkins University
% email address: shubhayu@jhu.edu
% March 2019; Last revision: 06-Apr-2020
%% ------------- BEGIN CODE --------------

%Procure imputed feature data and corresponding thresholds

addpath('functions/')
tic
load('../motion_feature_data/imputed_complete_sensor_data.mat');
load('../motion_feature_data/feature_thresholds.mat');
toc

dim_of_sensors=size(sensors);
feature_thresholds=num2cell(feature_thresholds');

%% Bed Motion Correction Algorithm:
% 1. In SMA domain, locate periods of time when bed contributes significant
% motion
% 2. Save those indices as markers for bed motion
% 3. Analyze SMA profiles directly preceeding those markers. If any motion
% sensor's SMA profile is significantly active prior to bed motion, then we
% 4. If no motion directly preceeds bed motion in those sensors, subtract
% bed SMA from those points.

bp_nm_range = [0 feature_thresholds(1)];
fe_nm_range = [feature_thresholds(2) 1.707];
fp1_nm_range = [0 feature_thresholds(3)];
fp2_nm_range = [0 feature_thresholds(4)];
mf_nm_range = [feature_thresholds(5) 3.2];
sma_nm_range = [0 feature_thresholds(6)];
wav_nm_range = [0 feature_thresholds(7)];

nm_ranges = {bp_nm_range,fe_nm_range,fp1_nm_range,fp2_nm_range,...
    mf_nm_range,sma_nm_range,wav_nm_range};

bed_corrected_sensors=cell(dim_of_sensors(1)-1,dim_of_sensors(2));

bedActIdx=cellfun(@(x,y)(x>=y),sensors(1,:),feature_thresholds,...
    'UniformOutput',false);
bedActIdx{2} = ~bedActIdx{2};
bedActIdx{5} = ~bedActIdx{5};

for i = 1:length(bed_corrected_sensors(1,:))
    curr_sens = sensors(2:length(sensors(1,:)),i);
    curr_bedIdx = bedActIdx{i};
    curr_bedMat = sensors{1,i};
    curr_range = nm_ranges{i};
    for j = 1:length(curr_sens)
        curr_mat=curr_sens{j};
        if ismember(feature_names(i),["freq_entropy","med_freq"])
            curr_mat(curr_bedIdx)=min(curr_mat(curr_bedIdx)+...
                curr_bedMat(curr_bedIdx),(curr_range{2}-curr_range{1})*....
                rand(size(curr_mat(curr_bedIdx)))+curr_range{1});
            curr_sens{j}=curr_mat;
        else
            curr_mat(curr_bedIdx)=max(curr_mat(curr_bedIdx)-...
                curr_bedMat(curr_bedIdx),feature_thresholds{i}*....
                rand(size(curr_mat(curr_bedIdx))));
            curr_sens{j}=curr_mat;
        end
    end
    bed_corrected_sensors(:,i)=curr_sens;
end

save('../motion_feature_data/bed_corrected_imputed_complete_sensor_data.mat',...
    'bed_corrected_sensors','feature_names');
