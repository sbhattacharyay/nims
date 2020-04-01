function sensors_output= impute_rmngMissingData(sensors,quasi_threshold,...
    trainWindow,missing_percentages,missingIdxs,feature_names,....
    feature_thresholds,t)

rng(1,'twister')

tic
rmngMissing=cellfun(@(x) x<quasi_threshold & x~=0, missing_percentages,...
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
        currRmngMissing=rmngMissing{sensIdx,featIdx};
        locA=find(currRmngMissing);
        currIdxArray=missingIdxs{sensIdx,featIdx}(currRmngMissing,:);
        draws=rand(size(currIdxArray));
        for k = 1:length(locA)
            temp_Idx=find(currIdxArray(k,:)');
            missingTimes = t(temp_Idx');
            nonMissingTimes=t(~currIdxArray(k,:));
            [~,I]=arrayfun(@(x) mink(abs(x-nonMissingTimes),...
                trainWindow*12),missingTimes,'UniformOutput',false);
            [~,trainIdx]=cellfun(@(x) ismember(nonMissingTimes(x),t),I,...
                'UniformOutput',false);
            
            if ismember(feature_names(featIdx),["freq_entropy","med_freq"])
                temp_pis = cellfun(@(x)sum(curr_mat(locA(k),x)>...
                    feature_thresholds(featIdx))/(trainWindow*12),trainIdx);
                nm_imputes=arrayfun(@(x,y)draws(k,y)<x,...
                    temp_pis,temp_Idx);
                exp_imputes=~nm_imputes;
                curr_mat(locA(k),temp_Idx(nm_imputes))=(curr_range(2)-...
                    curr_range(1))*rand(1,sum(nm_imputes))+curr_range(1);
                if sum(exp_imputes) ~= 0
                    temp_exp_fits=cellfun(@(x)fitdist((-curr_mat(locA(k),x(...
                        curr_mat(locA(k),x)<feature_thresholds(featIdx)))....
                        +feature_thresholds(featIdx))','Exponential'),.....
                        trainIdx(exp_imputes));
                    exp_imp_values=arrayfun(@(x)-random(x,1,1)+...
                        feature_thresholds(featIdx),temp_exp_fits)';
                    curr_mat(locA(k),temp_Idx(exp_imputes))=exp_imp_values;
                end
            else
                temp_pis = cellfun(@(x)sum(curr_mat(locA(k),x)<...
                    feature_thresholds(featIdx))/(trainWindow*12),trainIdx);
                nm_imputes=arrayfun(@(x,y)draws(k,y)<x,...
                    temp_pis,temp_Idx);
                exp_imputes=~nm_imputes;
                curr_mat(locA(k),temp_Idx(nm_imputes))=(curr_range(2)-...
                    curr_range(1))*rand(1,sum(nm_imputes))+curr_range(1);
                if sum(exp_imputes) ~= 0
                    temp_exp_fits=cellfun(@(x)fitdist((curr_mat(locA(k),x(...
                        curr_mat(locA(k),x)>feature_thresholds(featIdx)))....
                        -feature_thresholds(featIdx))','Exponential'),.....
                        trainIdx(exp_imputes));
                    exp_imp_values=arrayfun(@(x)random(x,1,1)+...
                        feature_thresholds(featIdx),temp_exp_fits)';                    
                    curr_mat(locA(k),temp_Idx(exp_imputes))=exp_imp_values;
                end
            end
        end
        sensors{sensIdx,featIdx} = curr_mat;        
        sensors_output{sensIdx,featIdx} = curr_mat;
    end
end

toc
end

