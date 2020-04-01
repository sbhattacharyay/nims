function sensors_output = impute_totallyMissingData(sensors,...
    patient_table,totallyMissingIdxs,feature_names,feature_thresholds,....
    studyPatientsPY,sortedPY,t)

tic
rng(1,'twister')

sensors_output=sensors;
dim_TM = size(totallyMissingIdxs);

bp_nm_range = [0 feature_thresholds(1)];
fe_nm_range = [feature_thresholds(2) 1.707];
fp1_nm_range = [0 feature_thresholds(3)];
fp2_nm_range = [0 feature_thresholds(4)];
mf_nm_range = [feature_thresholds(5) 3.2];
sma_nm_range = [0 feature_thresholds(6)];
wav_nm_range = [0 feature_thresholds(7)];

nm_ranges = {bp_nm_range,fe_nm_range,fp1_nm_range,fp2_nm_range,...
    mf_nm_range,sma_nm_range,wav_nm_range};

for i = 1:dim_TM(2)
    curr_sensorIdx = totallyMissingIdxs(1,i);
    curr_featureIdx = totallyMissingIdxs(2,i);
    
    curr_pt_PYIdx = totallyMissingIdxs(3,i);
    curr_pt_PY = sortedPY(totallyMissingIdxs(3,i));
    
    curr_pt_tableIdx = find(studyPatientsPY == curr_pt_PY);
    
    pi_var = ['pi_' char(feature_names(curr_featureIdx)) '_' ...
        num2str(curr_sensorIdx)];
    lambda_var = ['lambda_' char(feature_names(curr_featureIdx)) '_' ...
        num2str(curr_sensorIdx)];
    
    curr_pi=patient_table{curr_pt_tableIdx,pi_var};
    curr_lambda=patient_table{curr_pt_tableIdx,lambda_var};
    
    draws=rand(1,length(t));
    curr_range = nm_ranges{curr_featureIdx};
    
    NM_imputes=draws<curr_pi;
    NM_count = sum(NM_imputes);
    exp_imputes=~NM_imputes;
    exp_count = sum(exp_imputes);
    
    curr_mat=sensors{curr_sensorIdx,curr_featureIdx};
    
    curr_mat(curr_pt_PYIdx,NM_imputes)=(curr_range(2)-curr_range(1))*...
        rand(1,NM_count)+curr_range(1);
    
    if ismember(feature_names(curr_featureIdx),["freq_entropy","med_freq"])
        curr_mat(curr_pt_PYIdx,exp_imputes)=curr_range(1)-random(...
            'Exponential',curr_lambda,1,exp_count);
    else
        curr_mat(curr_pt_PYIdx,exp_imputes)=random('Exponential',...
            curr_lambda,1,exp_count) + curr_range(2);
    end
    sensors{curr_sensorIdx,curr_featureIdx}=curr_mat;
    sensors_output{curr_sensorIdx,curr_featureIdx}=curr_mat;
end
toc
end

