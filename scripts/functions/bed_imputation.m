function patient_table_output = bed_imputation(...
    patient_table,bed_Idx,totallyMissingIdxs,feature_names)

patient_table_output=patient_table;

bed_TMI = totallyMissingIdxs(:,ismember(totallyMissingIdxs(1,:),bed_Idx));
temp_dim = size(bed_TMI);

for i=1:temp_dim(2)
    curr_sensorIdx = bed_TMI(1,i);
    curr_featureIdx = bed_TMI(2,i);
    
    imputeVar1 = ['pi_' char(feature_names(curr_featureIdx)) '_' ...
        num2str(curr_sensorIdx)];
    imputeVar2 = ['lambda_' char(feature_names(curr_featureIdx)) '_' ...
        num2str(curr_sensorIdx)];
    
    curr_miss_Idx1 = isnan(patient_table.(imputeVar1));
    curr_miss_Idx2 = isnan(patient_table.(imputeVar2));
    
    patient_table_output{curr_miss_Idx1,imputeVar1}=nanmean(...
        patient_table.(imputeVar1));
    
    patient_table_output{curr_miss_Idx2,imputeVar2}=nanmean(...
        patient_table.(imputeVar2));
end

