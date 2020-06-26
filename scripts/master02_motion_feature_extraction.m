%% Master Script 2: Complete Motion Feature Extraction
% Decoding Quantitative Motor Features for Classification and Prediction
% in Severe Acquired Brain Injury
%
% Shubhayu Bhattacharyay
% Department of Biomedical Engineering
% Department of Applied Mathematics and Statistics
% Whiting School of Engineering, Johns Hopkins University
% email address: shubhayu@jhu.edu

%% Set Directory and Procure File Names
tic
addpath('functions/')

patientData = readtable('../clinical_data/SB_patient_table.xlsx',...
    'TreatAsEmpty',{'.','NA'});
patientData = sortrows(patientData,'AccelPatientNo_');

try
    cd('~/data/accel_sensor_data')
    marcc = true;
catch
    cd('../accel_sensor_data')
    marcc = false;
end

d = dir('data*');

studyPatients = cellfun(@(x) (string(x(5:6))),{d.name}');
studyDirs = {d.name}';

toc
%% Arrange data by time and cut out extraneous times
%WARNING: Elapsed Time: ~111.286133 seconds.
%WARNING: Elapsed Time: ~1104.244113 seconds for mega streams.

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
    
    toc
    %% Feature Extraction
    %Warning! Run Time: >937.769720 seconds.
    tic
    
    SMA= {};
    freqPairs = {};
    medF = {};
    bandPower = {};
    freqEnt = {};
    wavelets = {};
    
    windowingStart = datenum(dateshift(datetime(min(cellfun(@(x) x(1), ...
        data_copy(4,:))),'ConvertFrom', 'datenum'),'start','minute',....
        5));
    windowingEnd = datenum(dateshift(datetime(max(cellfun(@(x) x(end),...
        data_copy(4,:))),'ConvertFrom', 'datenum'),'start','minute',....
        -5));
    windowSize = 5; %in seconds
    window = windowSize/86400;
    
    timeSplit = windowingStart:window:windowingEnd;
    
    binCount = length(timeSplit)-1;
    binIdx = 1:binCount;
    
    for j = 1:length(data_copy)

        split_Data = cell(binCount,4);
        
        x_stream = data_copy{1,j};
        y_stream = data_copy{2,j};
        z_stream = data_copy{3,j};
        t_stream = data_copy{4,j};
        
        [NCounts, ~, indices] = histcounts(t_stream,timeSplit);
        totalIdx = indices';
        
        nonMP = NCounts>40;        
        availBins = find(nonMP);
        
        split_Data(nonMP,1) = arrayfun(@(curr_bin) x_stream(totalIdx==curr_bin),availBins,'UniformOutput', false);
        split_Data(nonMP,2) = arrayfun(@(curr_bin) y_stream(totalIdx==curr_bin),availBins,'UniformOutput', false);
        split_Data(nonMP,3) = arrayfun(@(curr_bin) z_stream(totalIdx==curr_bin),availBins,'UniformOutput', false);
        split_Data(nonMP,4) = arrayfun(@(curr_bin) t_stream(totalIdx==curr_bin),availBins,'UniformOutput', false);
        lens = cellfun(@(x) length(x),split_Data(:,1));
        times = datestr(timeSplit(2:end)');
                
        X = split_Data(:,1);
        Y = split_Data(:,2);
        Z = split_Data(:,3);
        T = split_Data(:,4);
        
        maskx = binIdx(cellfun(@(x) ~isempty(x), X));
        masky = binIdx(cellfun(@(y) ~isempty(y), Y));
        maskz = binIdx(cellfun(@(z) ~isempty(z), Z));
        maskt = binIdx(cellfun(@(t) ~isempty(t), T));
        
        superMask = intersect(maskx,masky);
        superMask = intersect(superMask,maskz);
        superMask = intersect(superMask,maskt);
        
        %SMA
        curr_sma = NaN(binCount,1);
        curr_sma(superMask) = ((window)^-1).*cellfun(@(x,y,z,t) trapz(t,abs(x)+abs(y)+abs(z)), ...
            X(superMask),Y(superMask),Z(superMask),T(superMask));
        SMA = [SMA {curr_sma;times;lens}];
        %Frequency component median pairs
        
        curr_Med = NaN(binCount,2);
        fc = 2.5;
        fs = 10;
        [b,a] = butter(4,fc/(fs/2),'high');
        [d,c] = butter(4,fc/(fs/2),'low');
        
        curr_Med(superMask,1) = cellfun(@(x,y,z,t) rssq([median(filter(d,c,x)),...
            median(filter(d,c,y)),median(filter(d,c,z))]),....
            X(superMask),Y(superMask),Z(superMask),T(superMask));
        
        curr_Med(superMask,2) = cellfun(@(x,y,z,t) rssq([median(filter(b,a,x)),...
            median(filter(b,a,y)),median(filter(b,a,z))]),....
            X(superMask),Y(superMask),Z(superMask),T(superMask));
        
        freqPairs = [freqPairs {curr_Med;times;lens}];
        
        %Median Freq
        curr_medFreq = NaN(binCount,1);
        curr_medFreq(superMask)=cellfun(@(x,y,z,t) rssq([medfreq(x),medfreq(y),medfreq(z)]), ...
            X(superMask),Y(superMask),Z(superMask),T(superMask));
        
        medF = [medF {curr_medFreq;times;lens}];
        
        %Band Power
        curr_bandPower = NaN(binCount,1);
        curr_bandPower(superMask)=cellfun(@(x,y,z,t) rssq([bandpower(x,fs,[0.3 3.5]), ...
            bandpower(y,fs,[0.3 3.5]),bandpower(z,fs,[0.3 3.5])]),....
            X(superMask),Y(superMask),Z(superMask),T(superMask));
        
        bandPower = [bandPower {curr_bandPower;times;lens}];
        
        %Frequency-Domain Entropy
        curr_freqEnt = NaN(binCount,1);
        curr_freqEnt(superMask)=cellfun(@(x,y,z) rssq([pentropy(x,fs,'Instantaneous',false),...
            pentropy(y,fs,'Instantaneous',false),pentropy(z,fs,'Instantaneous',false)]), ....
            X(superMask),Y(superMask),Z(superMask));
        
        freqEnt = [freqEnt {curr_freqEnt;times;lens}];
        %wavelets
        if marcc == true
            cd('~/data/accel_sensor_data')
        else
            cd('../accel_sensor_data')
        end
        
        curr_wvlt = NaN(binCount,1);
        curr_wvlt(superMask)=cellfun(@(x,y,z) get_wavelets(x,y,z), ...
            X(superMask),Y(superMask),Z(superMask));
        
        wavelets = [wavelets {curr_wvlt;times;lens}];
    end
    
    if marcc == true
        save(['~/data/all_motion_feature_data/band_power/band_power' folder_of_interest(30:31) '.mat'],'bandPower','-v7.3');
        save(['~/data/all_motion_feature_data/freq_entropy/freq_entropy' folder_of_interest(30:31) '.mat'],'freqEnt','-v7.3');
        save(['~/data/all_motion_feature_data/freq_pairs/freq_pairs' folder_of_interest(30:31) '.mat'],'freqPairs','-v7.3');
        save(['~/data/all_motion_feature_data/med_freq/med_freq' folder_of_interest(30:31) '.mat'],'medF','-v7.3');
        save(['~/data/all_motion_feature_data/sma/sma' folder_of_interest(30:31) '.mat'],'SMA','-v7.3');
        save(['~/data/all_motion_feature_data/wavelets/wavelets' folder_of_interest(30:31) '.mat'],'wavelets','-v7.3');
        disp(['Patient No. ' folder_of_interest(30:31) ' completed.']);
    else
        save(['../all_motion_feature_data/band_power/band_power' folder_of_interest(26:27) '.mat'],'bandPower','-v7.3');
        save(['../all_motion_feature_data/freq_entropy/freq_entropy' folder_of_interest(26:27) '.mat'],'freqEnt','-v7.3');
        save(['../all_motion_feature_data/freq_pairs/freq_pairs' folder_of_interest(26:27) '.mat'],'freqPairs','-v7.3');
        save(['../all_motion_feature_data/med_freq/med_freq' folder_of_interest(26:27) '.mat'],'medF','-v7.3');
        save(['../all_motion_feature_data/sma/sma' folder_of_interest(26:27) '.mat'],'SMA','-v7.3');
        save(['../all_motion_feature_data/wavelets/wavelets' folder_of_interest(26:27) '.mat'],'wavelets','-v7.3');
        disp(['Patient No. ' folder_of_interest(26:27) ' completed.']);
    end
    toc
end

%% Collection of Motion Features into consolidated files:

cd ../all_motion_feature_data/
directory = dir;
dirNames = string({directory.name});
no_feature_logicals=([directory.isdir] & dirNames ~= "." & ...
    dirNames ~= "..");
feature_names = dirNames(no_feature_logicals);

subDirectory = dir(feature_names(1));
subDirNames = string({subDirectory.name});
no_patient_logicals=endsWith(subDirNames,".mat");

sensors = cell(sum(no_patient_logicals),length(feature_names));
times = cell(sum(no_patient_logicals),1);
lens = cell(sum(no_patient_logicals),1);

placement_labels = ["Bed","LA","LE","LW","RA","RE","RW"]';

for i = 1:length(feature_names)
    folder = [pwd filesep char(feature_names(i))];
    temp_dir = dir(folder);
    tempDirNames = string({temp_dir.name});
    matFileIdxs = find((tempDirNames ~= "." & tempDirNames ~= ".." & ...
        endsWith(tempDirNames,".mat")));
    
    loaded_featCells = arrayfun(@(x) struct2cell(load([folder ...
        filesep temp_dir(x).name])), matFileIdxs,'UniformOutput',false);
    loaded_featMats = cellfun(@(x) x{1},loaded_featCells,...
        'UniformOutput',false);
    [~,n_sensors] = size(loaded_featMats{1});
    for j = 1:length(loaded_featMats)
        curr_featMats = loaded_featMats{j};
        [t_length, corr_scale] = size(curr_featMats{1,1});
        times{j} = curr_featMats{2,1}';
        lens{j} = curr_featMats{3,1}';
        compiled_matrix = [];
        for k = 1:(n_sensors)
            compiled_matrix = [compiled_matrix; curr_featMats{1,k}'];
        end
        sensors{j,i} = compiled_matrix;
    end
end


% Organizing the frequency pairs dataset

%At this point the frequency pairs are all combined in the same (nx2)x12959
%matrix. However, we would like to split up the high and low frequency
%matrices. Thus, the odd rows will be partitioned to one matrix while the
%even rows to another:

[n,~] = size(sensors{1,3});

freqPair1 = cellfun(@(x) x(1:2:(n-1),:),sensors(:,3),'UniformOutput',false);
freqPair2 = cellfun(@(x) x(2:2:n,:),sensors(:,3),'UniformOutput',false);

sensors(:,5:7)=sensors(:,4:6);
sensors(:,3)=freqPair1;
sensors(:,4)=freqPair2;

feature_names = [feature_names(1:2),strcat(feature_names(3),"1"),...
    strcat(feature_names(3),"2"),feature_names(4:6)];

%  Saving the data
save('complete_sensor_data.mat','sensors');
save('times.mat','times');
save('lens.mat','lens');
writematrix(placement_labels,"placement_labels.csv");
writematrix(feature_names,"feature_names.csv");

%% Correction of Time Formatting
times = load('times.mat');
times = times.times;

formatted_times = cellfun(@(x) string(x'),times,'UniformOutput',false);
indexed_times = cell2table(cell(0,2),'VariableNames',{'times','ptIdx'});


for i = 1:length(formatted_times)
   curr_cell = formatted_times{i}; 
   new_cell = [curr_cell,i*ones(length(curr_cell),1)];
   indexed_times = [indexed_times ; array2table(new_cell,'VariableNames',{'times','ptIdx'})];
end

indexed_times.ptIdx = str2double(indexed_times.ptIdx);

writetable(indexed_times,'indexed_times.csv');

cd ../scripts/

license('inuse')

%% Characterization of Missing Data

load('../all_motion_feature_data/complete_sensor_data.mat')
load('../all_motion_feature_data/formatted_times.mat')
addpath('functions/')

[~,sensors] = get_totallyMissingIdxs(sensors);

missing_data = zeros(size(sensors,1), size(sensors,2));
recording_time = cell(size(formatted_times,1),3);

% Calculate the start & end times of each patient recording and the
% duration of recording
for i = 1:size(formatted_times,1)
    recording_time{i,1} = formatted_times{i}(1);
    recording_time{i,2} = formatted_times{i}(end);
    recording_time{i,3} = char(datetime(recording_time{i,2}, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss')... 
    - datetime(recording_time{i,1}, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss'));
end

% Find the fraction of missing data
for row = 1:size(sensors,1)
    curr_data = sensors{row,1};
    for i = 1:size(curr_data,1)
        curr_sensor_data = curr_data(i,:);
        nan_sum = sum(isnan(curr_sensor_data));

        iszero = find(curr_sensor_data == 0);
        count = 0;
        next = 0;
        if length(iszero) > 2
            for j = 1:length(iszero)-2
                if j < next
                    continue
                end
                if (iszero(j+1) - iszero(j) == 1) && (iszero(j+2) - iszero(j+1) == 1)
                    seq = 3;
                    if j == length(iszero)-2
                        count = count + seq;
                        continue
                    end
                    while iszero(j+seq) - iszero(j+seq-1) == 1
                        seq = seq + 1;
                        if iszero(j+seq-1) == iszero(end)
                            break
                        end
                    end
                    count = count + seq;
                    next = j + seq;
                else
                    continue
                end   
            end
        end
        
        tot_sum = nan_sum + count;
        missing_data(row,i) = tot_sum / length(curr_sensor_data);
    end 
end

d = dir('../accel_sensor_data/data*');

studyPatients = str2double(cellfun(@(x) (string(x(5:6))),{d.name}'));
studyPatients = array2table(studyPatients);
studyPatients.Properties.VariableNames = {'Study Patient No.'};

table1 = array2table(missing_data);
table1.Properties.VariableNames = {'Bed','LA','LE','LW','RA','RE','RW'};

table2 = cell2table(recording_time);
table2.Properties.VariableNames = {'Start Timestamp','End Timestamp','Recording Duration'};

finalTable = [studyPatients,table1,table2];
%writetable(finalTable,'../all_motion_feature_data/MissingPercentTables.xlsx');