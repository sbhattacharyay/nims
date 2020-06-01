function output_patient_table = fit_zip(sensors,patient_table,...
    feature_thresholds,feature_names)

dim_of_sensors=size(sensors);

output_patient_table=patient_table;

tic
all_pis = {};
all_lambdas = {};

for j = 1:dim_of_sensors(2)
    temp_feature_data = sensors(:,j);
    k = feature_thresholds(j);
    if j == 2 || j == 5
        all_pis{j,1}=cellfun(@(x) noMotion_th_neg(x,k), ...
            temp_feature_data,'UniformOutput',false);
        all_lambdas{j,1}=cellfun(@(x) fit_exp(x,k), ...
            temp_feature_data,'UniformOutput',false);
    else
        all_pis{j,1}=cellfun(@(x) noMotion_th(x,k), ...
            temp_feature_data,'UniformOutput',false);
        all_lambdas{j,1}=cellfun(@(x) fit_exp_neg(x,k), ...
            temp_feature_data,'UniformOutput',false);
    end
end

for i = 1:dim_of_sensors(2)
    for j = 1:dim_of_sensors(1)
        curr_col_pi = all_pis{i,1}{j,1};
        new_var = sprintf("pi_%s_%d",feature_names(i),j);
        output_patient_table=addvars(output_patient_table,curr_col_pi,...
            'NewVariableNames',new_var);
    end
end


for i = 1:dim_of_sensors(2)
    for j = 1:dim_of_sensors(1)
        curr_col_lam = all_lambdas{i,1}{j,1};
        new_var = sprintf("lambda_%s_%d",feature_names(i),j);
        output_patient_table=addvars(output_patient_table,curr_col_lam,...
            'NewVariableNames',new_var);
    end
end
toc

end

