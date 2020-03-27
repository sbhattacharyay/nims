function patient_table_output = wrist_elbow_Regression(patient_table,...
    wrist_Idx,elbow_Idx,totallyMissingIdxs,feature_names)

patient_table_output=patient_table;

% Totally missing indices with a missing arm (elbow or wrist) parameter
arm_TMI = totallyMissingIdxs(:,ismember(totallyMissingIdxs(1,:),[wrist_Idx,elbow_Idx]));
temp_dim = size(arm_TMI);

for i=1:temp_dim(2)
    curr_sensorIdx = arm_TMI(1,i);
    curr_featureIdx = arm_TMI(2,i);
    if ismember(curr_sensorIdx,wrist_Idx)
        
        responseVar1 = ['pi_' char(feature_names(curr_featureIdx)) '_' ...
            num2str(curr_sensorIdx)];
        responseVar2 = ['lambda_' char(feature_names(curr_featureIdx)) '_' ...
            num2str(curr_sensorIdx)];
        predictorVar1 = ['pi_' char(feature_names(curr_featureIdx)) '_' ...
            num2str(curr_sensorIdx-1)];
        predictorVar2 = ['lambda_' char(feature_names(curr_featureIdx)) '_' ...
            num2str(curr_sensorIdx-1)];
        mdl1=fitglm(patient_table,'Distribution','binomial','ResponseVar',...
            responseVar1,'PredictorVars',predictorVar1);
        mdl2=fitglm(patient_table,'Distribution','poisson','ResponseVar',...
            responseVar2,'PredictorVars',predictorVar2);
        curr_miss_Idx1 = isnan(patient_table.(responseVar1));
        curr_miss_Idx2 = isnan(patient_table.(responseVar2));
        
        patient_table_output{curr_miss_Idx1,responseVar1}=predict(mdl1,...
            patient_table(curr_miss_Idx1,{predictorVar1}));
        patient_table_output{curr_miss_Idx2,responseVar2}=predict(mdl2,...
            patient_table(curr_miss_Idx2,{predictorVar2}));        
        
    elseif ismember(curr_sensorIdx,elbow_Idx)
        responseVar1 = ['pi_' char(feature_names(curr_featureIdx)) '_' ...
            num2str(curr_sensorIdx)];
        responseVar2 = ['lambda_' char(feature_names(curr_featureIdx)) '_' ...
            num2str(curr_sensorIdx)];
        predictorVar1 = ['pi_' char(feature_names(curr_featureIdx)) '_' ...
            num2str(curr_sensorIdx+1)];
        predictorVar2 = ['lambda_' char(feature_names(curr_featureIdx)) '_' ...
            num2str(curr_sensorIdx+1)];
        mdl1=fitglm(patient_table,'Distribution','binomial','ResponseVar',...
            responseVar1,'PredictorVars',predictorVar1);
        mdl2=fitglm(patient_table,'Distribution','poisson','ResponseVar',...
            responseVar2,'PredictorVars',predictorVar2);
        curr_miss_Idx1 = isnan(patient_table.(responseVar1));
        curr_miss_Idx2 = isnan(patient_table.(responseVar2));
        
        patient_table_output{curr_miss_Idx1,responseVar1}=predict(mdl1,...
            patient_table(curr_miss_Idx1,{predictorVar1}));
        patient_table_output{curr_miss_Idx2,responseVar2}=predict(mdl2,...
            patient_table(curr_miss_Idx2,{predictorVar2}));   
        
    end
end
end

