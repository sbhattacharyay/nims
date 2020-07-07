% fid = fopen('../all_motion_feature_data/featureset_ptidx_timestamp.txt');
% featureset = textscan(fid, '%f%s%s', 'HeaderLines', 1);
% ptIdx = featureset{1,1};
% 
% timeStamp = cell(length(featureset{1,2}),1);
% for i = 1:length(timeStamp)
%    timeStamp(i) = datenum([featureset{1,2}{i,1} ' ' featureset{1,3}{i,1}], 'yyyy-mm-dd HH:MM:SS');
%    timeStamp{i} = [featureset{1,2}{i,1} ' ' featureset{1,3}{i,1}];
% end
% 
% save('../all_motion_feature_data/featureset_ptidx_timestamp.mat', 'ptIdx', 'timeStamp')

clear; clc
load('../all_motion_feature_data/featureset_ptidx_timestamp.mat')

% Change here
time_before = 5/60/24;
[~,~,gcs_data] = xlsread('../clinical_data/GCS_table.xlsx');
header = gcs_data(1,:);
gcs_pt = cell2mat(gcs_data(2:end,1));
gcs_date = datenum(gcs_data(2:end,2), 'dd/mm/yyyy') + cell2mat(gcs_data(2:end,9));
accel_overlap = logical(cell2mat(gcs_data(2:end,8)));
gcs_data = gcs_data(accel_overlap,:);
gcs_pt = gcs_pt(accel_overlap);
gcs_date = gcs_date(accel_overlap);

[~,~,clinical_data] = xlsread('../clinical_data/patient_clinical_data.csv');
accel_conv = cell2mat(clinical_data(2:70,2));
study_conv = cell2mat(clinical_data(2:70,1));

% for j = 1:length(ptIdx)
%     index = find(accel_conv == ptIdx(j));
%     if isempty(index)
%         continue
%     else
%         ptIdx(j) = study_conv(index);
%     end
% end

for j = 1:length(ptIdx)
   ptIdx(j) = study_conv(ptIdx(j));
end

% for k = 1:length(gcs_pt)
%     if isempty(find(study_conv == gcs_pt(k)))
%         continue
%     else
%         gcs_pt(k) = find(study_conv == gcs_pt(k));
%     end
% end

uniqueid = 1:length(gcs_pt);
ptIdx_uniqueid = zeros(length(ptIdx),1);

for k = 1:length(gcs_pt)
   curr_pt = gcs_pt(k);
   curr_date = gcs_date(k);
   curr_unique_id = uniqueid(k);
   
   index1 = logical(ptIdx == curr_pt);
   time2 = curr_date - time_before;
   index2 = logical(timeStamp >= time2);
   index3 = logical(timeStamp <= curr_date);
   
   tot_index = logical((index1 + index2 + index3) == 3);
   if ~isempty(find(tot_index))
       ptIdx_uniqueid(find(tot_index)) = curr_unique_id;
   end
end

uniqueid = table2cell(array2table(uniqueid'));
gcs_data = horzcat(gcs_data, uniqueid);
header = horzcat(header, {'Unique ID'});
gcs_data = vertcat(header, gcs_data);

save('../all_motion_feature_data/gcs_featureset_id.mat', 'gcs_data', 'ptIdx', 'ptIdx_uniqueid', 'timeStamp')