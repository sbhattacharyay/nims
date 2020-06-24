function patient_table_output = ankle_Regression(patient_table,...
    ankle_Idx,totallyMissingIdxs,feature_names)

patient_table_output=patient_table;

% Totally missing indices with a missing ankle parameter
ankle_TMI = totallyMissingIdxs(:,ismember(totallyMissingIdxs(1,:),ankle_Idx));
temp_dim = size(ankle_TMI);

for i=1:temp_dim(2)
    curr_sensorIdx = ankle_TMI(1,i);
    curr_featureIdx = ankle_TMI(2,i);
    
    responseVar1 = ['pi_' char(feature_names(curr_featureIdx)) '_' ...
        num2str(curr_sensorIdx)];
    responseVar2 = ['lambda_' char(feature_names(curr_featureIdx)) '_' ...
        num2str(curr_sensorIdx)];
    
    predictorVar1_1 = ['pi_' char(feature_names(curr_featureIdx)) '_' ...
        num2str(curr_sensorIdx+1)];
    predictorVar1_2 = ['pi_' char(feature_names(curr_featureIdx)) '_' ...
        num2str(curr_sensorIdx+2)];
    
    predictorVar2_1 = ['lambda_' char(feature_names(curr_featureIdx)) '_' ...
        num2str(curr_sensorIdx+1)];
    predictorVar2_2 = ['lambda_' char(feature_names(curr_featureIdx)) '_' ...
        num2str(curr_sensorIdx+2)];
    
    mdl1=fitglm(patient_table,'interactions','Distribution','binomial',...
        'ResponseVar',responseVar1,'PredictorVars',{predictorVar1_1,....
        predictorVar1_2});
    mdl2=fitglm(patient_table,'interactions','Distribution','poisson',...
        'ResponseVar',responseVar2,'PredictorVars',{predictorVar2_1,....
        predictorVar2_2});
    
    curr_miss_Idx1 = isnan(patient_table.(responseVar1));
    curr_miss_Idx2 = isnan(patient_table.(responseVar2));
    
    patient_table_output{curr_miss_Idx1,responseVar1}=predict(mdl1,...
        patient_table(curr_miss_Idx1,{predictorVar1_1,predictorVar1_2}));
    
    patient_table_output{curr_miss_Idx2,responseVar2}=predict(mdl2,...
        patient_table(curr_miss_Idx2,{predictorVar2_1,predictorVar2_2}));
end
end

