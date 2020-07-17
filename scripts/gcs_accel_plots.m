%% Plot 3 hours before GCS recording
% Decoding Quantitative Motor Features for Classification and Prediction
% in Severe Acquired Brain Injury
%
% Shubhayu Bhattacharyay
% Department of Biomedical Engineering
% Department of Applied Mathematics and Statistics
% Whiting School of Engineering, Johns Hopkins University
% email address: shubhayu@jhu.edu

%% Load patient clinical data:

time_preceding = 3; % in hours

dt_window = hours(time_preceding);
no_pts = 3*60*60*10 + 1;

patient_clinical_data = sortrows(readtable('../clinical_data/patient_clinical_data.csv'),'AccelPatientNo_');
ptIdx = (1:height(patient_clinical_data))';
patient_clinical_data = addvars(patient_clinical_data,ptIdx,'After','AccelPatientNo_');

gcs_data = readtable('../clinical_data/GCS_table.xlsx');
gcs_data.timeStamps = datetime(gcs_data.Date) + days(gcs_data.Time_1);

accel_data_files = string({dir('../pure_accel_data/*.csv').name}');

gcs_t_pur_data = array2table(zeros(no_pts,819));
gcs_t_abs_data = array2table(zeros(no_pts,819));
repScore = repelem(strcat("GCSt_",string(3:15))',63,1);
repSensors = repmat(repelem(["Bed","LA","LE","LW","RA","RE","RW"]',9,1),13,1);
repAxes = repmat(repelem(["x","y","z"]',3,1),91,1);
repMeas = repmat(["m","s","n"]',273,1);
gcs_t_pur_data.Properties.VariableNames = join([repScore,repSensors,repAxes,repMeas],"_");
gcs_t_abs_data.Properties.VariableNames = join([repScore,repSensors,repAxes,repMeas],"_");

gcs_m_pur_data = array2table(zeros(no_pts,378));
gcs_m_abs_data = array2table(zeros(no_pts,378));
repScore = repelem(strcat("GCSm_",string(1:6))',63,1);
repSensors = repmat(repelem(["Bed","LA","LE","LW","RA","RE","RW"]',9,1),6,1);
repAxes = repmat(repelem(["x","y","z"]',3,1),42,1);
repMeas = repmat(["m","s","n"]',126,1);
gcs_m_pur_data.Properties.VariableNames = join([repScore,repSensors,repAxes,repMeas],"_");
gcs_m_abs_data.Properties.VariableNames = join([repScore,repSensors,repAxes,repMeas],"_");

gcs_e_pur_data = array2table(zeros(no_pts,252));
gcs_e_abs_data = array2table(zeros(no_pts,252));
repScore = repelem(strcat("GCSe_",string(1:4))',63,1);
repSensors = repmat(repelem(["Bed","LA","LE","LW","RA","RE","RW"]',9,1),4,1);
repAxes = repmat(repelem(["x","y","z"]',3,1),28,1);
repMeas = repmat(["m","s","n"]',84,1);
gcs_e_pur_data.Properties.VariableNames = join([repScore,repSensors,repAxes,repMeas],"_");
gcs_e_abs_data.Properties.VariableNames = join([repScore,repSensors,repAxes,repMeas],"_");

gcs_v_pur_data = array2table(zeros(no_pts,315));
gcs_v_abs_data = array2table(zeros(no_pts,315));
repScore = repelem(strcat("GCSv_",string(1:5))',63,1);
repSensors = repmat(repelem(["Bed","LA","LE","LW","RA","RE","RW"]',9,1),5,1);
repAxes = repmat(repelem(["x","y","z"]',3,1),35,1);
repMeas = repmat(["m","s","n"]',105,1);
gcs_v_pur_data.Properties.VariableNames = join([repScore,repSensors,repAxes,repMeas],"_");
gcs_v_abs_data.Properties.VariableNames = join([repScore,repSensors,repAxes,repMeas],"_");

missingID=[0 NaN Inf];

for ptIdx = 1:length(accel_data_files)
    studyNo = patient_clinical_data.StudyPatientNo_(ptIdx);
    disp(['Study patient no. ',num2str(studyNo),' initated.']);
    currPtData = readtable(fullfile('..','pure_accel_data',accel_data_files(ptIdx)));
    currPtData.timeStamps = table2array(varfun(@(x) datetime(x,...
        'InputFormat','dd-MMM-yyyy HH:mm:ss:SSS','Format','dd-MMM-yyyy HH:mm:ss:SSS'),currPtData,'InputVariables','timeStamps'));
    % Correct totally missing data as NaN
    for w = 1:21
        if all(ismissing(currPtData{:,w},missingID))
            currPtData{:,w} = nan(height(currPtData),1);
        end
    end
    
    filt_GCS = gcs_data(gcs_data.StudyPatientNumber == studyNo,:);
    
    for rowIdx = 1:height(filt_GCS)
        
        currGCStimestamp = filt_GCS.timeStamps(rowIdx);
        currEye = filt_GCS.Eye(rowIdx);
        currVerbal = filt_GCS.Verbal(rowIdx);
        currMotor = filt_GCS.Motor(rowIdx);
        currTotal = filt_GCS.TotalGCS(rowIdx);
        
        timeStamps = ((currGCStimestamp - dt_window):(seconds(0.1)):currGCStimestamp)';
        timeStamps.Format = 'dd-MMM-yyyy HH:mm:ss:SSS';
        segIdx = (1:length(timeStamps))';
        
        if sum(sum(ismember(currPtData.timeStamps,timeStamps))) == 0
            continue
        end
        
        orgData = table(timeStamps,segIdx);
        filt_PtData = outerjoin(currPtData,orgData,'Type','right','Keys','timeStamps');
        abs_filt_ptData = filt_PtData;
        abs_filt_ptData(:,1:21) = varfun(@abs,abs_filt_ptData,'InputVariables',1:21);
        
        for sensIdx = 1:21
            currStream = filt_PtData{:,sensIdx};
            currStreamName = filt_PtData.Properties.VariableNames{sensIdx};
            nmCurrStream = ~ismissing(currStream);
            nmCount = sum(nmCurrStream);
            
            if all(~nmCurrStream)
                continue
            end
            
            gcs_t_name = ['GCSt_',num2str(currTotal),'_',currStreamName];
            gcs_m_name = ['GCSm_',num2str(currMotor),'_',currStreamName];
            gcs_e_name = ['GCSe_',num2str(currEye),'_',currStreamName];
            gcs_v_name = ['GCSv_',num2str(currVerbal),'_',currStreamName];
            
            % Pure Values
            gcs_t_pur_data{nmCurrStream,[gcs_t_name,'_n']} = ...
                nansum([gcs_t_pur_data{nmCurrStream,[gcs_t_name,....
                '_n']},ones(nmCount,1)],2);
            std_condition = (gcs_t_pur_data{:,[gcs_t_name,'_n']}>1);
            oldMeans = gcs_t_pur_data{nmCurrStream&std_condition,[gcs_t_name,'_m']};
            newMeans = gcs_t_pur_data{nmCurrStream,[gcs_t_name,'_m']}...
                +(((gcs_t_pur_data{nmCurrStream,[gcs_t_name,'_n']}).^....
                (-1)).*(currStream(nmCurrStream)-gcs_t_pur_data.....
                {nmCurrStream,[gcs_t_name,'_m']}));
            if sum(std_condition) ~= 0
                gcs_t_pur_data{(nmCurrStream&std_condition),[gcs_t_name,...
                    '_s']} = gcs_t_pur_data{(nmCurrStream&std_condition),....
                    [gcs_t_name,'_s']} + (((gcs_t_pur_data{nmCurrStream& .....
                    std_condition,[gcs_t_name,'_n']}).^(-1)).*((currStream(nmCurrStream&std_condition) ......
                    - oldMeans).^2)) - (((gcs_t_pur_data{nmCurrStream&.......
                    std_condition,[gcs_t_name,'_n']}-1).^(-1)).*........
                    (gcs_t_pur_data{(nmCurrStream&std_condition),.........
                    [gcs_t_name,'_s']}));
            end
            gcs_t_pur_data{nmCurrStream,[gcs_t_name,'_m']} = newMeans;
            
            gcs_m_pur_data{nmCurrStream,[gcs_m_name,'_n']} = ...
                nansum([gcs_m_pur_data{nmCurrStream,[gcs_m_name,....
                '_n']},ones(nmCount,1)],2);
            std_condition = (gcs_m_pur_data{:,[gcs_m_name,'_n']}>1);
            oldMeans = gcs_m_pur_data{nmCurrStream&std_condition,[gcs_m_name,'_m']};
            newMeans = gcs_m_pur_data{nmCurrStream,[gcs_m_name,'_m']}...
                +(((gcs_m_pur_data{nmCurrStream,[gcs_m_name,'_n']}).^....
                (-1)).*(currStream(nmCurrStream)-gcs_m_pur_data.....
                {nmCurrStream,[gcs_m_name,'_m']}));
            if sum(std_condition) ~= 0
                gcs_m_pur_data{(nmCurrStream&std_condition),[gcs_m_name,...
                    '_s']} = gcs_m_pur_data{(nmCurrStream&std_condition),....
                    [gcs_m_name,'_s']} + (((gcs_m_pur_data{nmCurrStream& .....
                    std_condition,[gcs_m_name,'_n']}).^(-1)).*((currStream(nmCurrStream&std_condition) ......
                    -oldMeans).^2)) - (((gcs_m_pur_data{nmCurrStream&.......
                    std_condition,[gcs_m_name,'_n']}-1).^(-1)).*........
                    (gcs_m_pur_data{(nmCurrStream&std_condition),.........
                    [gcs_m_name,'_s']}));
            end
            gcs_m_pur_data{nmCurrStream,[gcs_m_name,'_m']} = newMeans;
            
            gcs_e_pur_data{nmCurrStream,[gcs_e_name,'_n']} = ...
                nansum([gcs_e_pur_data{nmCurrStream,[gcs_e_name,....
                '_n']},ones(nmCount,1)],2);
            std_condition = (gcs_e_pur_data{:,[gcs_e_name,'_n']}>1);
            oldMeans = gcs_e_pur_data{nmCurrStream&std_condition,[gcs_e_name,'_m']};
            newMeans = gcs_e_pur_data{nmCurrStream,[gcs_e_name,'_m']}...
                +(((gcs_e_pur_data{nmCurrStream,[gcs_e_name,'_n']}).^....
                (-1)).*(currStream(nmCurrStream)-gcs_e_pur_data.....
                {nmCurrStream,[gcs_e_name,'_m']}));
            if sum(std_condition) ~= 0
                gcs_e_pur_data{(nmCurrStream&std_condition),[gcs_e_name,...
                    '_s']} = gcs_e_pur_data{(nmCurrStream&std_condition),....
                    [gcs_e_name,'_s']} + (((gcs_e_pur_data{nmCurrStream& .....
                    std_condition,[gcs_e_name,'_n']}).^(-1)).*((currStream(nmCurrStream&std_condition) ......
                    -oldMeans).^2)) - (((gcs_e_pur_data{nmCurrStream&.......
                    std_condition,[gcs_e_name,'_n']}-1).^(-1)).*........
                    (gcs_e_pur_data{(nmCurrStream&std_condition),.........
                    [gcs_e_name,'_s']}));
            end
            gcs_e_pur_data{nmCurrStream,[gcs_e_name,'_m']} = newMeans;
            
            gcs_v_pur_data{nmCurrStream,[gcs_v_name,'_n']} = ...
                nansum([gcs_v_pur_data{nmCurrStream,[gcs_v_name,....
                '_n']},ones(nmCount,1)],2);
            std_condition = (gcs_v_pur_data{:,[gcs_v_name,'_n']}>1);
            oldMeans = gcs_v_pur_data{nmCurrStream&std_condition,[gcs_v_name,'_m']};
            newMeans = gcs_v_pur_data{nmCurrStream,[gcs_v_name,'_m']}...
                +(((gcs_v_pur_data{nmCurrStream,[gcs_v_name,'_n']}).^....
                (-1)).*(currStream(nmCurrStream)-gcs_v_pur_data.....
                {nmCurrStream,[gcs_v_name,'_m']}));
            if sum(std_condition) ~= 0
                gcs_v_pur_data{(nmCurrStream&std_condition),[gcs_v_name,...
                    '_s']} = gcs_v_pur_data{(nmCurrStream&std_condition),....
                    [gcs_v_name,'_s']} + (((gcs_v_pur_data{nmCurrStream& .....
                    std_condition,[gcs_v_name,'_n']}).^(-1)).*((currStream(nmCurrStream&std_condition) ......
                    -oldMeans).^2)) - (((gcs_v_pur_data{nmCurrStream&.......
                    std_condition,[gcs_v_name,'_n']}-1).^(-1)).*........
                    (gcs_v_pur_data{(nmCurrStream&std_condition),.........
                    [gcs_v_name,'_s']}));
            end
            gcs_v_pur_data{nmCurrStream,[gcs_v_name,'_m']} = newMeans;
            
            % Absolute Values
            gcs_t_abs_data{nmCurrStream,[gcs_t_name,'_n']} = ...
                nansum([gcs_t_abs_data{nmCurrStream,[gcs_t_name,....
                '_n']},ones(nmCount,1)],2);
            std_condition = (gcs_t_abs_data{:,[gcs_t_name,'_n']}>1);
            oldMeans = gcs_t_abs_data{nmCurrStream&std_condition,[gcs_t_name,'_m']};
            newMeans = gcs_t_abs_data{nmCurrStream,[gcs_t_name,'_m']}....
                +(((gcs_t_abs_data{nmCurrStream,[gcs_t_name,'_n']}).^(-1))....
                .*(abs(currStream(nmCurrStream))-gcs_t_abs_data{nmCurrStream,[gcs_t_name,'_m']}));
            if sum(std_condition) ~= 0
                gcs_t_abs_data{(nmCurrStream&std_condition),[gcs_t_name,...
                    '_s']} = gcs_t_abs_data{(nmCurrStream&std_condition),....
                    [gcs_t_name,'_s']} + (((gcs_t_abs_data{nmCurrStream& .....
                    std_condition,[gcs_t_name,'_n']}).^(-1)).*((abs(currStream(nmCurrStream&std_condition)) ......
                    -oldMeans).^2)) - (((gcs_t_abs_data{nmCurrStream&.......
                    std_condition,[gcs_t_name,'_n']}-1).^(-1)).*........
                    (gcs_t_abs_data{(nmCurrStream&std_condition),.........
                    [gcs_t_name,'_s']}));
            end
            gcs_t_abs_data{nmCurrStream,[gcs_t_name,'_m']} = newMeans;
            
            gcs_m_abs_data{nmCurrStream,[gcs_m_name,'_n']} = ...
                nansum([gcs_m_abs_data{nmCurrStream,[gcs_m_name,....
                '_n']},ones(nmCount,1)],2);
            std_condition = (gcs_m_abs_data{:,[gcs_m_name,'_n']}>1);
            oldMeans = gcs_m_abs_data{nmCurrStream&std_condition,[gcs_m_name,'_m']};
            newMeans = gcs_m_abs_data{nmCurrStream,[gcs_m_name,'_m']}....
                +(((gcs_m_abs_data{nmCurrStream,[gcs_m_name,'_n']}).^(-1))....
                .*(abs(currStream(nmCurrStream))-gcs_m_abs_data{nmCurrStream,[gcs_m_name,'_m']}));
            if sum(std_condition) ~= 0
                gcs_m_abs_data{(nmCurrStream&std_condition),[gcs_m_name,...
                    '_s']} = gcs_m_abs_data{(nmCurrStream&std_condition),....
                    [gcs_m_name,'_s']} + (((gcs_m_abs_data{nmCurrStream& .....
                    std_condition,[gcs_m_name,'_n']}).^(-1)).*((abs(currStream(nmCurrStream&std_condition)) ......
                    -oldMeans).^2)) - (((gcs_m_abs_data{nmCurrStream&.......
                    std_condition,[gcs_m_name,'_n']}-1).^(-1)).*........
                    (gcs_m_abs_data{(nmCurrStream&std_condition),.........
                    [gcs_m_name,'_s']}));
            end
            gcs_m_abs_data{nmCurrStream,[gcs_m_name,'_m']} = newMeans;
            
            gcs_e_abs_data{nmCurrStream,[gcs_e_name,'_n']} = ...
                nansum([gcs_e_abs_data{nmCurrStream,[gcs_e_name,....
                '_n']},ones(nmCount,1)],2);
            std_condition = (gcs_e_abs_data{:,[gcs_e_name,'_n']}>1);
            oldMeans = gcs_e_abs_data{nmCurrStream&std_condition,[gcs_e_name,'_m']};
            newMeans = gcs_e_abs_data{nmCurrStream,[gcs_e_name,'_m']}....
                +(((gcs_e_abs_data{nmCurrStream,[gcs_e_name,'_n']}).^(-1))....
                .*(abs(currStream(nmCurrStream))-gcs_e_abs_data{nmCurrStream,[gcs_e_name,'_m']}));
            if sum(std_condition) ~= 0
                gcs_e_abs_data{(nmCurrStream&std_condition),[gcs_e_name,...
                    '_s']} = gcs_e_abs_data{(nmCurrStream&std_condition),....
                    [gcs_e_name,'_s']} + (((gcs_e_abs_data{nmCurrStream& .....
                    std_condition,[gcs_e_name,'_n']}).^(-1)).*((abs(currStream(nmCurrStream&std_condition)) ......
                    -oldMeans).^2)) - (((gcs_e_abs_data{nmCurrStream&.......
                    std_condition,[gcs_e_name,'_n']}-1).^(-1)).*........
                    (gcs_e_abs_data{(nmCurrStream&std_condition),.........
                    [gcs_e_name,'_s']}));
            end
            gcs_e_abs_data{nmCurrStream,[gcs_e_name,'_m']} = newMeans;
            
            gcs_v_abs_data{nmCurrStream,[gcs_v_name,'_n']} = ...
                nansum([gcs_v_abs_data{nmCurrStream,[gcs_v_name,....
                '_n']},ones(nmCount,1)],2);
            std_condition = (gcs_v_abs_data{:,[gcs_v_name,'_n']}>1);
            oldMeans = gcs_v_abs_data{nmCurrStream&std_condition,[gcs_v_name,'_m']};
            newMeans = gcs_v_abs_data{nmCurrStream,[gcs_v_name,'_m']}....
                +(((gcs_v_abs_data{nmCurrStream,[gcs_v_name,'_n']}).^(-1))....
                .*(abs(currStream(nmCurrStream))-gcs_v_abs_data{nmCurrStream,[gcs_v_name,'_m']}));
            if sum(std_condition) ~= 0
                gcs_v_abs_data{(nmCurrStream&std_condition),[gcs_v_name,...
                    '_s']} = gcs_v_abs_data{(nmCurrStream&std_condition),....
                    [gcs_v_name,'_s']} + (((gcs_v_abs_data{nmCurrStream& .....
                    std_condition,[gcs_v_name,'_n']}).^(-1)).*((abs(currStream(nmCurrStream&std_condition)) ......
                    -oldMeans).^2)) - (((gcs_v_abs_data{nmCurrStream&.......
                    std_condition,[gcs_v_name,'_n']}-1).^(-1)).*........
                    (gcs_v_abs_data{(nmCurrStream&std_condition),.........
                    [gcs_v_name,'_s']}));
            end
            gcs_v_abs_data{nmCurrStream,[gcs_v_name,'_m']} = newMeans;
            
        end
    end
    
    clear currPtData
    disp(['Study patient no. ',num2str(studyNo),' complete.']);
end

% Correct variance to standard deviation

gcs_t_pur_data{:,2:3:width(gcs_t_pur_data)} = sqrt(gcs_t_pur_data{:,2:3:width(gcs_t_pur_data)});
gcs_t_abs_data{:,2:3:width(gcs_t_abs_data)} = sqrt(gcs_t_abs_data{:,2:3:width(gcs_t_abs_data)});

gcs_m_pur_data{:,2:3:width(gcs_m_pur_data)} = sqrt(gcs_m_pur_data{:,2:3:width(gcs_m_pur_data)});
gcs_m_abs_data{:,2:3:width(gcs_m_abs_data)} = sqrt(gcs_m_abs_data{:,2:3:width(gcs_m_abs_data)});

gcs_e_pur_data{:,2:3:width(gcs_e_pur_data)} = sqrt(gcs_e_pur_data{:,2:3:width(gcs_e_pur_data)});
gcs_e_abs_data{:,2:3:width(gcs_e_abs_data)} = sqrt(gcs_e_abs_data{:,2:3:width(gcs_e_abs_data)});

gcs_v_pur_data{:,2:3:width(gcs_v_pur_data)} = sqrt(gcs_v_pur_data{:,2:3:width(gcs_v_pur_data)});
gcs_v_abs_data{:,2:3:width(gcs_v_abs_data)} = sqrt(gcs_v_abs_data{:,2:3:width(gcs_v_abs_data)});

%% Save data

save('../pure_accel_data/gcs_t_pur_data.mat','gcs_t_pur_data')
save('../pure_accel_data/gcs_t_abs_data.mat','gcs_t_abs_data')

save('../pure_accel_data/gcs_m_pur_data.mat','gcs_m_pur_data')
save('../pure_accel_data/gcs_m_abs_data.mat','gcs_m_abs_data')

save('../pure_accel_data/gcs_e_pur_data.mat','gcs_e_pur_data')
save('../pure_accel_data/gcs_e_abs_data.mat','gcs_e_abs_data')

save('../pure_accel_data/gcs_v_pur_data.mat','gcs_v_pur_data')
save('../pure_accel_data/gcs_v_abs_data.mat','gcs_v_abs_data')

%% Load motion features

time_preceding = 3; % in hours

dt_window = hours(time_preceding);
no_pts = time_preceding*60*12 + 1;

featNames = ["band_power","freq_entropy","freq_pairs1","freq_pairs2","med_freq","sma","wavelets"];

patient_clinical_data = sortrows(readtable('../clinical_data/patient_clinical_data.csv'),'AccelPatientNo_');
ptIdx = (1:height(patient_clinical_data))';
patient_clinical_data = addvars(patient_clinical_data,ptIdx,'After','AccelPatientNo_');

gcs_data = readtable('../clinical_data/GCS_table.xlsx');
gcs_data.timeStamps = datetime(gcs_data.Date) + days(gcs_data.Time_1);

completeSensorData=load('../all_motion_feature_data/complete_sensor_data.mat');
completeSensorData=completeSensorData.sensors;
completeSensorTime=load('../all_motion_feature_data/formatted_times.mat');
completeSensorTime=completeSensorTime.formatted_times;

gcs_t_mf_data = array2table(zeros(no_pts,1911));
repScore = repelem(strcat("GCSt_",string(3:15))',147,1);
repSensors = repmat(repelem(["Bed","LA","LE","LW","RA","RE","RW"]',21,1),13,1);
repFeats = repmat(repelem(featNames',3,1),91,1);
repMeas = repmat(["m","s","n"]',637,1);
gcs_t_mf_data.Properties.VariableNames = join([repScore,repSensors,repFeats,repMeas],"_");

gcs_m_mf_data = array2table(zeros(no_pts,882));
repScore = repelem(strcat("GCSm_",string(1:6))',147,1);
repSensors = repmat(repelem(["Bed","LA","LE","LW","RA","RE","RW"]',21,1),6,1);
repFeats = repmat(repelem(featNames',3,1),42,1);
repMeas = repmat(["m","s","n"]',294,1);
gcs_m_mf_data.Properties.VariableNames = join([repScore,repSensors,repFeats,repMeas],"_");

gcs_e_mf_data = array2table(zeros(no_pts,588));
repScore = repelem(strcat("GCSe_",string(1:4))',147,1);
repSensors = repmat(repelem(["Bed","LA","LE","LW","RA","RE","RW"]',21,1),4,1);
repFeats = repmat(repelem(featNames',3,1),28,1);
repMeas = repmat(["m","s","n"]',196,1);
gcs_e_mf_data.Properties.VariableNames = join([repScore,repSensors,repFeats,repMeas],"_");

gcs_v_mf_data = array2table(zeros(no_pts,735));
repScore = repelem(strcat("GCSv_",string(1:5))',147,1);
repSensors = repmat(repelem(["Bed","LA","LE","LW","RA","RE","RW"]',21,1),5,1);
repFeats = repmat(repelem(featNames',3,1),35,1);
repMeas = repmat(["m","s","n"]',245,1);
gcs_v_mf_data.Properties.VariableNames = join([repScore,repSensors,repFeats,repMeas],"_");

missingID=[0 NaN Inf];

for ptIdx = 1:length(completeSensorData)
    studyNo = patient_clinical_data.StudyPatientNo_(ptIdx);
    filt_GCS = gcs_data(gcs_data.StudyPatientNumber == studyNo,:);
    disp(['Study patient no. ',num2str(studyNo),' initated.']);
    currTmData = datetime(completeSensorTime{ptIdx},'InputFormat',...
        'dd-MMM-yyyy HH:mm:ss','Format','dd-MMM-yyyy HH:mm:ss');
    for rowIdx = 1:height(filt_GCS)
        currGCStimestamp = filt_GCS.timeStamps(rowIdx);
        currEye = filt_GCS.Eye(rowIdx);
        currVerbal = filt_GCS.Verbal(rowIdx);
        currMotor = filt_GCS.Motor(rowIdx);
        currTotal = filt_GCS.TotalGCS(rowIdx);
        timeStamps = ((currGCStimestamp - dt_window):(seconds(5)):currGCStimestamp)';
        segIdx = (1:length(timeStamps))';
        if sum(sum(ismember(currTmData,timeStamps))) == 0
            continue
        end
        %disp(['GCS row no. ',num2str(rowIdx),' intiated.']);
        orgData = table(timeStamps,segIdx);
        for ftIdx = 1:length(featNames)
            currFtName = featNames(ftIdx);
            currPtData = addvars(array2table(completeSensorData{ptIdx,ftIdx}'),currTmData);
            currPtData.Properties.VariableNames(1:7) = ...
                join([["Bed","LA","LE","LW", ....
                "RA","RE","RW"]',repelem(currFtName,7,1)],"_");
            currPtData.Properties.VariableNames(8)="timeStamps";
            % Correct totally missing data as NaN
            for w = 1:(width(currPtData)-1)
                if all(ismissing(currPtData{:,w},missingID))
                    currPtData{:,w} = nan(height(currPtData),1);
                end
            end
            filt_PtData = outerjoin(currPtData,orgData,'Type','right','Keys','timeStamps','MergeKeys',true);
            for sensIdx = 1:7
                currStream = filt_PtData{:,sensIdx};
                currStreamName = filt_PtData.Properties.VariableNames{sensIdx};
                nmCurrStream = ~ismissing(currStream);
                nmCount = sum(nmCurrStream);
                if all(~nmCurrStream)
                    continue
                end
                gcs_t_name = ['GCSt_',num2str(currTotal),'_',currStreamName];
                gcs_m_name = ['GCSm_',num2str(currMotor),'_',currStreamName];
                gcs_e_name = ['GCSe_',num2str(currEye),'_',currStreamName];
                gcs_v_name = ['GCSv_',num2str(currVerbal),'_',currStreamName];
                
                old_n = gcs_t_mf_data{nmCurrStream,[gcs_t_name,'_n']};
                new_n = old_n + 1;
                std_condition = new_n>1;
                oldMeans = gcs_t_mf_data{nmCurrStream,[gcs_t_name,'_m']};
                newMeans = oldMeans+((1./new_n).*(currStream(nmCurrStream)-oldMeans));
                if sum(std_condition) ~= 0
                    nmIdx = find(nmCurrStream);
                    sdIdx = nmIdx(std_condition);
                    oldVar=gcs_t_mf_data{sdIdx,[gcs_t_name,'_s']};
                    newVar=oldVar+((1./new_n(std_condition)).* ...
                        ((currStream(sdIdx)-oldMeans(std_condition)).^2))....
                        -((1./old_n(std_condition)).*oldVar);
                    gcs_t_mf_data{sdIdx,[gcs_t_name,'_s']} = newVar;
                end
                gcs_t_mf_data{nmCurrStream,[gcs_t_name,'_m']} = newMeans;
                gcs_t_mf_data{nmCurrStream,[gcs_t_name,'_n']} = new_n;
                
                old_n = gcs_m_mf_data{nmCurrStream,[gcs_m_name,'_n']};
                new_n = old_n + 1;
                std_condition = new_n>1;
                oldMeans = gcs_m_mf_data{nmCurrStream,[gcs_m_name,'_m']};
                newMeans = oldMeans+((1./new_n).*(currStream(nmCurrStream)-oldMeans));
                if sum(std_condition) ~= 0
                    nmIdx = find(nmCurrStream);
                    sdIdx = nmIdx(std_condition);
                    oldVar=gcs_m_mf_data{sdIdx,[gcs_m_name,'_s']};
                    newVar=oldVar+((1./new_n(std_condition)).* ...
                        ((currStream(sdIdx)-oldMeans(std_condition)).^2))....
                        -((1./old_n(std_condition)).*oldVar);
                    gcs_m_mf_data{sdIdx,[gcs_m_name,'_s']} = newVar;
                end
                gcs_m_mf_data{nmCurrStream,[gcs_m_name,'_m']} = newMeans;
                gcs_m_mf_data{nmCurrStream,[gcs_m_name,'_n']} = new_n;
                
                
                old_n = gcs_e_mf_data{nmCurrStream,[gcs_e_name,'_n']};
                new_n = old_n + 1;
                std_condition = new_n>1;
                oldMeans = gcs_e_mf_data{nmCurrStream,[gcs_e_name,'_m']};
                newMeans = oldMeans+((1./new_n).*(currStream(nmCurrStream)-oldMeans));
                if sum(std_condition) ~= 0
                    nmIdx = find(nmCurrStream);
                    sdIdx = nmIdx(std_condition);
                    oldVar=gcs_e_mf_data{sdIdx,[gcs_e_name,'_s']};
                    newVar=oldVar+((1./new_n(std_condition)).* ...
                        ((currStream(sdIdx)-oldMeans(std_condition)).^2))....
                        -((1./old_n(std_condition)).*oldVar);
                    gcs_e_mf_data{sdIdx,[gcs_e_name,'_s']} = newVar;
                end
                gcs_e_mf_data{nmCurrStream,[gcs_e_name,'_m']} = newMeans;
                gcs_e_mf_data{nmCurrStream,[gcs_e_name,'_n']} = new_n;
                
                old_n = gcs_v_mf_data{nmCurrStream,[gcs_v_name,'_n']};
                new_n = old_n + 1;
                std_condition = new_n>1;
                oldMeans = gcs_v_mf_data{nmCurrStream,[gcs_v_name,'_m']};
                newMeans = oldMeans+((1./new_n).*(currStream(nmCurrStream)-oldMeans));
                if sum(std_condition) ~= 0
                    nmIdx = find(nmCurrStream);
                    sdIdx = nmIdx(std_condition);
                    oldVar=gcs_v_mf_data{sdIdx,[gcs_v_name,'_s']};
                    newVar=oldVar+((1./new_n(std_condition)).* ...
                        ((currStream(sdIdx)-oldMeans(std_condition)).^2))....
                        -((1./old_n(std_condition)).*oldVar);
                    gcs_v_mf_data{sdIdx,[gcs_v_name,'_s']} = newVar;
                end
                gcs_v_mf_data{nmCurrStream,[gcs_v_name,'_m']} = newMeans;
                gcs_v_mf_data{nmCurrStream,[gcs_v_name,'_n']} = new_n;
            end
        end
    end
    clear currPtData
    disp(['Study patient no. ',num2str(studyNo),' complete.']);
end

% Correct variance to standard deviation
gcs_t_mf_data{:,2:3:width(gcs_t_mf_data)} = sqrt(gcs_t_mf_data{:,2:3:width(gcs_t_mf_data)});
gcs_m_mf_data{:,2:3:width(gcs_m_mf_data)} = sqrt(gcs_m_mf_data{:,2:3:width(gcs_m_mf_data)});
gcs_e_mf_data{:,2:3:width(gcs_e_mf_data)} = sqrt(gcs_e_mf_data{:,2:3:width(gcs_e_mf_data)});
gcs_v_mf_data{:,2:3:width(gcs_v_mf_data)} = sqrt(gcs_v_mf_data{:,2:3:width(gcs_v_mf_data)});

%% Save data

save('../all_motion_feature_data/gcs_t_mf_data.mat','gcs_t_mf_data')
save('../all_motion_feature_data/gcs_m_mf_data.mat','gcs_m_mf_data')
save('../all_motion_feature_data/gcs_e_mf_data.mat','gcs_e_mf_data')
save('../all_motion_feature_data/gcs_v_mf_data.mat','gcs_v_mf_data')