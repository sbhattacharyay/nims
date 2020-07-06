% fid = fopen('../all_motion_feature_data/featureset_ptidx_timestamp.txt');
% featureset = textscan(fid, '%f%s%s', 'HeaderLines', 1);
% ptIdx = featureset{1,1};
% 
% timeStamp = zeros(length(featureset{1,2}),1);
% for i = 1:length(timeStamp)
%    timeStamp(i) = datenum([featureset{1,2}{i,1} ' ' featureset{1,3}{i,1}], 'yyyy-mm-dd HH:MM:SS');
% end
% 
% save('../all_motion_feature_data/featureset_ptidx_timestamp.mat', 'ptIdx', 'timeStamp')


load('../all_motion_feature_data/featureset_ptidx_timestamp.mat')
time_before = 5/60/24;
[~,~,gcs_data] = xlsread('../clinical_data/GCS_table.xlsx');
gcs_pt = cell2mat(gcs_data(2:end,1));
gcs_date = datenum(gcs_data(2:end,2), 'dd/mm/yyyy') + cell2mat(gcs_data(2:end,9));
accel_overlap = logical(cell2mat(gcs_data(2:end,8)));
gcs_pt = gcs_pt(accel_overlap);
gcs_date = gcs_date(accel_overlap);

[~,~,clinical_data] = xlsread('../clinical_data/patient_clinical_data.xlsx');
study_accel_conv = cell2mat(clinical_data(2:73,1:2));

for j = 1:length(gcs_pt)
   index = find(study_accel_conv(:,1) == gcs_pt(j));
   gcs_pt(j) = study_accel_conv(index,2);
end

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
sum(ptIdx_uniqueid)