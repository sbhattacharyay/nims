function sensors_output = impute_quasiMissingData(sensors,...
    patient_table,quasi_threshold,missing_percentages,missingIdxs,....
    feature_names,feature_thresholds,studyPatientsPY,sortedPY)

tic
rng(1,'twister')

quasiMissing=cellfun(@(x) x>quasi_threshold & x~=1, missing_percentages,...
    'UniformOutput', false);
dim_of_sensors=size(sensors);
sensors_output=sensors;

bp_nm_range = [0 feature_thresholds(1)];
fe_nm_range = [feature_thresholds(2) 1.707];
fp1_nm_range = [0 feature_thresholds(3)];
fp2_nm_range = [0 feature_thresholds(4)];
mf_nm_range = [feature_thresholds(5) 3.2];
sma_nm_range = [0 feature_thresholds(6)];
wav_nm_range = [0 feature_thresholds(7)];

nm_ranges = {bp_nm_range,fe_nm_range,fp1_nm_range,fp2_nm_range,...
    mf_nm_range,sma_nm_range,wav_nm_range};


for sensIdx = 1:dim_of_sensors(1)
    for featIdx = 1:dim_of_sensors(2)
        curr_mat = sensors{sensIdx,featIdx};
        curr_range = nm_ranges{featIdx};
        pi_nam=['pi_' char(feature_names(featIdx)) '_' num2str(sensIdx)];
        lambda_nam=['lambda_',char(feature_names(featIdx)),'_',...
            num2str(sensIdx)];
        currQuasiMissing=quasiMissing{sensIdx,featIdx};
        locA=find(currQuasiMissing);
        [~,locB]=ismember(sortedPY(currQuasiMissing),studyPatientsPY);
        pi_vars=patient_table{locB,pi_nam};
        lambda_vars=patient_table{locB,lambda_nam};
        currIdxArray=missingIdxs{sensIdx,featIdx}(currQuasiMissing,:);
        draws=rand(size(currIdxArray));
        for k = 1:length(pi_vars)
            nm_imputes = currIdxArray(k,:) & (draws(k,:)<pi_vars(k));
            exp_imputes = currIdxArray(k,:) & ~nm_imputes;
            curr_mat(locA(k),nm_imputes)=(curr_range(2)-curr_range(1))*...
                rand(1,sum(nm_imputes))+curr_range(1);
            if ismember(feature_names(featIdx),["freq_entropy","med_freq"])
                curr_mat(locA(k),exp_imputes)=curr_range(1)-random(...
                    'Exponential',lambda_vars(k),1,sum(exp_imputes));
            else
                curr_mat(locA(k),exp_imputes)=random('Exponential',...
                    lambda_vars(k),1,sum(exp_imputes)) + curr_range(2);
            end 
        end
        sensors_output{sensIdx,featIdx} = curr_mat;
    end
end
toc
end