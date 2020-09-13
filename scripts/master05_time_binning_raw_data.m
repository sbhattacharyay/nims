%% Master Script 5: Time binning of raw accelerometery data
% Decoding Quantitative Motor Features for Classification and Prediction
% in Severe Acquired Brain Injury
%
% Shubhayu Bhattacharyay
% Department of Biomedical Engineering
% Department of Applied Mathematics and Statistics
% Whiting School of Engineering, Johns Hopkins University
% email address: shubhayu@jhu.edu

%% Set directory and procure file names
tic
addpath('functions/')

patientData = readtable('../clinical_data/patient_clinical_data.xlsx',...
    'TreatAsEmpty',{'.','NA'});
patientData = sortrows(patientData,'AccelPatientNo_');

try
    cd('~/data/accel_sensor_data')
    marcc = true;
    mkdir('~/data/pure_accel_data')
catch
    cd('../accel_sensor_data')
    marcc = false;
    mkdir('../pure_accel_data')
end

d = dir('data*');

studyPatients = cellfun(@(x) (string(x(5:6))),{d.name}');
studyDirs = {d.name}';

toc

%% Organise data into clean 0.1 second chunks

for patIdx = 1:length(studyDirs)
    tic
    
    if marcc == true
        folder_of_interest = ['~/data/accel_sensor_data/' studyDirs{patIdx}];
        disp(['Patient No. ' folder_of_interest(30:31) ' initiated.']);
        curr_AccelPatientNo = str2double(folder_of_interest(30:31));
    else
        folder_of_interest = ['../accel_sensor_data/',studyDirs{patIdx}];
        disp(['Patient No. ' folder_of_interest(26:27) ' initiated.']);
        curr_AccelPatientNo = str2double(folder_of_interest(26:27));
    end
    
    table_row_idx = find(patientData.AccelPatientNo_ == ...
        curr_AccelPatientNo);
    
    try
        load(folder_of_interest)
    catch
        cd(folder_of_interest)
        load('C1.mat')
        load('C2.mat')
        load('C3.mat')
        load('C4.mat')
        load('C5.mat')
        load('C6.mat')
        load('C7.mat')
        data = [C1;C2;C3;C4;C5;C6;C7]';
        cd ..
    end
    toc
    %% Filter Data (bworth 4-th order)
    tic
    
    fc = 0.2;
    fs = 10;
    [b,a] = butter(4,fc/(fs/2),'high');
    data_copy = cell(size(data));
    
    for i = 1:length(data(1,:))
        start_date = patientData.AccelRecordingStartDate(table_row_idx);
        % Date correction for Patient 21, Left Elbow Sensor (#3)
        if patIdx == 20 && i == 3
            start_date = start_date+1;
        end
        curr = data(:,i);
        x = filter(b,a,curr{5});
        y = filter(b,a,curr{6});
        z = filter(b,a,curr{7});
        time = string(curr{12});
        time_cut = extractBefore(time,13);
        dt_info = strcat(string(start_date),{' '},time_cut);
        time_dn = datenum(dt_info,'dd-mmm-yyyy HH:MM:SS:FFF');
        time_diff = (diff(time_dn));
        [max_diffs,date_change_Idx] = findpeaks(-time_diff,'MinPeakHeight',.1,'MinPeakDistance',10);
        time_dn_fixed = time_dn;
        if length(max_diffs) >= 1
            date_change_Idx = date_change_Idx+1;
            disp(['Additional ' num2str(length(date_change_Idx)) ' Day(s) Detected for Current Sensor ' num2str(i)])
            for j = 1:length(date_change_Idx)
                time_dn_fixed(date_change_Idx(j):end) = time_dn_fixed(date_change_Idx(j):end)+1;
            end
        end
        data_copy{5,i} = x;
        data_copy{6,i} = y;
        data_copy{7,i} = z;
        data_copy{12,i}= time_dn_fixed;
    end
    
    %Delete Extraneous Rows
    data_copy(1:4,:) = [];
    data_copy(4:7,:) = [];
    clear data time time_cut time_dn time_diff x y z dt_info time_dn_fixed
    
    toc
    %% Segment entire time series
    tic
    
    windowingStart = datenum(dateshift(datetime(min(cellfun(@(x) x(1), ...
        data_copy(4,:))),'ConvertFrom', 'datenum'),'start','minute',....
        5));
    windowingEnd = datenum(dateshift(datetime(max(cellfun(@(x) x(end),...
        data_copy(4,:))),'ConvertFrom', 'datenum'),'start','minute',....
        -5));
    windowSize = 0.1; %in seconds
    window = windowSize/86400;
    
    timeSplit = windowingStart:window:windowingEnd;
    
    binCount = length(timeSplit)-1;
    toc
    
    compiledSplits = nan(binCount,(length(data_copy)*3));
    
    for j = 1:length(data_copy)
        split_Data = nan(binCount,3);
        
        x_stream = data_copy{1,j};
        y_stream = data_copy{2,j};
        z_stream = data_copy{3,j};
        t_stream = data_copy{4,j};
        
        [NCounts, ~, indices] = histcounts(t_stream,timeSplit);
        totalIdx = indices';
        
        singeValues = find(NCounts==1);
        split_Data(singeValues,1) = x_stream(ismember(indices,singeValues));
        split_Data(singeValues,2) = y_stream(ismember(indices,singeValues));
        split_Data(singeValues,3) = z_stream(ismember(indices,singeValues));
        
        disp(['Splitting multipe value indices for sensor no. ',num2str(j)])
        multValues = find(NCounts>=2);
        split_Data(multValues,1) = arrayfun(@(curr_bin) nanmean(x_stream(totalIdx==curr_bin)),multValues);
        split_Data(multValues,2) = arrayfun(@(curr_bin) nanmean(y_stream(totalIdx==curr_bin)),multValues);
        split_Data(multValues,3) = arrayfun(@(curr_bin) nanmean(z_stream(totalIdx==curr_bin)),multValues);
        disp(['Sensor no. ',num2str(j), ' complete'])
        
        compiledSplits(:,((j-1)*3+1):((j-1)*3+3)) = split_Data;
    end
    
    accelDataTable = array2table(compiledSplits);
    accelDataTable.timeStamps = datestr(timeSplit(2:end)','dd-mmm-yyyy HH:MM:SS:FFF');
    accelDataTable.ptIdx = patIdx*ones(height(accelDataTable),1);
    
    repSensors = repelem(["Bed","LA","LE","LW","RA","RE","RW"]',3,1);
    repAxes = repmat(["x","y","z"]',length(data_copy),1);
    
    accelDataTable.Properties.VariableNames(1:(length(data_copy)*3)) = join([repSensors,repAxes],".");
    
    if marcc == true
        writetable(accelDataTable,['~/data/pure_accel_data/puredata',num2str(curr_AccelPatientNo),'.csv']);
        disp(['Patient No. ' num2str(curr_AccelPatientNo) ' completed.']);
    else
        
        writetable(accelDataTable,['../pure_accel_data/puredata',num2str(curr_AccelPatientNo),'.csv']);
        disp(['Patient No. ' num2str(curr_AccelPatientNo) ' completed.']);
    end
end