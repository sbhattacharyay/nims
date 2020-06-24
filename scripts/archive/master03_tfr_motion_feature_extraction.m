%% Master Script 3: Motion Feature Extraction per Time from Recording
% Decoding Quantitative Motor Features for Classification and Prediction
% in Severe Acquired Brain Injury
%
% Shubhayu Bhattacharyay
% Department of Biomedical Engineering
% Department of Applied Mathematics and Statistics
% Whiting School of Engineering, Johns Hopkins University
% email address: shubhayu@jhu.edu
%% ------------- BEGIN CODE --------------
%Set Directory and Procure File Names
tic

addpath('functions/')

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

if marcc == true
    d2 = dir('~/data/tfr_motion_feature_data/band_power/*.mat');
    currentlyDone = cellfun(@(x) string(x(11:12)),{d2.name}');
else
    d2 = dir('../tfr_motion_feature_data/band_power/*.mat');
    currentlyDone = cellfun(@(x) string(x(11:12)),{d2.name}');
end

[featSet,diffIdx] = setdiff(studyPatients,currentlyDone);

toc
%% Arrange data by time and cut out extraneous times
%WARNING: Elapsed Time: ~111.286133 seconds.
%WARNING: Elapsed Time: ~1104.244113 seconds for mega streams.

for patIdx = diffIdx'
    tic
    
    if marcc == true
        folder_of_interest = ['~/data/accel_sensor_data/' studyDirs{patIdx}];
        disp(['Patient No. ' folder_of_interest(30:31) ' initiated.']);
    else
        folder_of_interest = ['../accel_sensor_data/',studyDirs{patIdx}];
        disp(['Patient No. ' folder_of_interest(26:27) ' initiated.']);
    end
    
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
    %%
    tic
    %Filter Data (bworth 4-th order)
    fc = 0.2;
    fs = 10;
    [b,a] = butter(4,fc/(fs/2),'high');    
    data_copy = data;
    score_1 = cellfun(@(x) double(length(x)>=864000), data(12,:));
    score_2 = cellfun(@(x) double(length(x)>=(864000*2)), data(12,:));
    dayLengths = score_1 + score_2;
    time_dn_fixed = {};
    for i = 1:length(data(1,:))
        curr = data(:,i);
        x = filter(b,a,curr{5});
        y = filter(b,a,curr{6});
        z = filter(b,a,curr{7});
        time = string(curr{12});
        time_cut = extractBefore(time,13);
        time_dn = datenum(time_cut,'HH:MM:SS:FFF');
        time_diff = (diff(time_dn));
        [max_diffs,date_change_Idx] = findpeaks(-time_diff,...
            'MinPeakHeight',.1,'MinPeakDistance',10);
        date_change_Idx = date_change_Idx+1;
        time_dn_fixed{i} = time_dn;
        if length(max_diffs) >= 1
            disp('Additional Day Detected for Current Sensor')
            for j=1:length(date_change_Idx)
                time_dn_fixed{i}(date_change_Idx(j):end) = ...
                    time_dn_fixed{i}(date_change_Idx(j):end)+1;
            end
        end
        data_copy{5,i} = x;
        data_copy{6,i} = y;
        data_copy{7,i} = z;
        data_copy{12,i}= time_dn_fixed{i};
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
    
    second = datenum('00:00:01:000','HH:MM:SS:FFF')-datenum('00:00:00:000','HH:MM:SS:FFF');
    windowSize = 5; %in seconds
    window = windowSize.*second;
    
    for j = 1:length(data_copy)
        windowingStart = time_dn_fixed{j}(1);
        windowingEnd = windowingStart+(1/3);
        
        timeSplit =windowingStart:window:windowingEnd;
        
        binCount = length(timeSplit)-1;
        binIdx = 1:binCount;
        
        split_Data = cell(binCount,4);
        D = data_copy{4,j};
        x_stream = data_copy{1,j};
        y_stream = data_copy{2,j};
        z_stream = data_copy{3,j};
        
        [NCounts, ~, indices] = histcounts(D,timeSplit);
        totalIdx = indices';
        
        nonMP = NCounts>40;
        rowSizes=(sum(nonMP,2));
        
        [~,streamList] = sort(rowSizes,'descend');
        
        emptyChoice =  arrayfun(@(idx) sum(nonMP(:,idx)) == 0,binIdx);
        singleChoice = arrayfun(@(idx) sum(nonMP(:,idx)) == 1,binIdx);
        multChoice = arrayfun(@(idx) sum(nonMP(:,idx)) >= 2,binIdx);
        
        singleIdx = binIdx(singleChoice);
        multIdx = binIdx(multChoice);
        
        singleRows = arrayfun(@(sIdx) (find(nonMP(:,sIdx)==1)),singleIdx);
        multRows = arrayfun(@(mIdx)find(nonMP(:,mIdx)==1),multIdx,'UniformOutput', false);
        multRowSelect = cellfun(@(mRows) intersect(streamList,mRows,'stable'),multRows,'UniformOutput',false);
        multRowSelect2= cellfun(@(mRows) mRows(1),multRowSelect);
        
        split_Data(singleIdx,1) = arrayfun(@(sRow,sIdx) x_stream(totalIdx(sRow,:)==sIdx),singleRows,singleIdx,'UniformOutput', false);
        split_Data(multIdx,1) = arrayfun(@(mRow,mIdx) x_stream(totalIdx(mRow,:)==mIdx),multRowSelect2,multIdx,'UniformOutput', false);
        
        split_Data(singleIdx,2) = arrayfun(@(sRow,sIdx) y_stream(totalIdx(sRow,:)==sIdx),singleRows,singleIdx,'UniformOutput', false);
        split_Data(multIdx,2) = arrayfun(@(mRow,mIdx) y_stream(totalIdx(mRow,:)==mIdx),multRowSelect2,multIdx,'UniformOutput', false);
        
        split_Data(singleIdx,3) = arrayfun(@(sRow,sIdx) z_stream(totalIdx(sRow,:)==sIdx),singleRows,singleIdx,'UniformOutput', false);
        split_Data(multIdx,3) = arrayfun(@(mRow,mIdx) z_stream(totalIdx(mRow,:)==mIdx),multRowSelect2,multIdx,'UniformOutput', false);
        
        split_Data(singleIdx,4) = arrayfun(@(sRow,sIdx) D(totalIdx(sRow,:)==sIdx),singleRows,singleIdx,'UniformOutput', false);
        split_Data(multIdx,4) = arrayfun(@(mRow,mIdx) D(totalIdx(mRow,:)==mIdx),multRowSelect2,multIdx,'UniformOutput', false);
        
        lens = cellfun(@(x) length(x),split_Data(:,1));
        times=timeSplit(2:end)';
                
        X = split_Data(:,1);
        Y = split_Data(:,2);
        Z = split_Data(:,3);
        W = split_Data(:,4);
        
        maskx = binIdx(cellfun(@(x) ~isempty(x), X));
        masky = binIdx(cellfun(@(y) ~isempty(y), Y));
        maskz = binIdx(cellfun(@(z) ~isempty(z), Z));
        maskw = binIdx(cellfun(@(w) ~isempty(w), W));
        
        superMask = intersect(maskx,masky);
        superMask = intersect(superMask,maskz);
        superMask = intersect(superMask,maskw);
        
        %SMA
        curr_sma = NaN(binCount,1);
        curr_sma(superMask) = ((window)^-1).*cellfun(@(x,y,z,w) trapz(w,abs(x)+abs(y)+abs(z)), ...
            X(superMask),Y(superMask),Z(superMask),W(superMask));
        SMA = [SMA {curr_sma;times;lens}];
        
        %Frequency component median pairs
        curr_Med = NaN(binCount,2);
        fc = 2.5;
        fs = 10;
        [d,c] = butter(4,fc/(fs/2),'low');
        [b,a] = butter(4,fc/(fs/2),'high');
                
        curr_Med(superMask,1)=cellfun(@(x,y,z,w) rssq([median(filter(d,c,x)),...
            median(filter(d,c,y)),median(filter(d,c,z))]),....
            X(superMask),Y(superMask),Z(superMask),W(superMask));
        
        curr_Med(superMask,2)=cellfun(@(x,y,z,w) rssq([median(filter(b,a,x)),...
            median(filter(b,a,y)),median(filter(b,a,z))]),....
            X(superMask),Y(superMask),Z(superMask),W(superMask));
        
        freqPairs = [freqPairs {curr_Med;times;lens}];
        
        %Median Freq
        curr_medFreq = NaN(binCount,1);
        curr_medFreq(superMask)=cellfun(@(x,y,z,w) rssq([medfreq(x),medfreq(y),medfreq(z)]), ...
            X(superMask),Y(superMask),Z(superMask),W(superMask));
        
        medF = [medF {curr_medFreq;times;lens}];
        
        %Band Power
        curr_bandPower = NaN(binCount,1);
        curr_bandPower(superMask)=cellfun(@(x,y,z,w) rssq([bandpower(x,fs,[0.3 3.5]), ...
            bandpower(y,fs,[0.3 3.5]),bandpower(z,fs,[0.3 3.5])]),....
            X(superMask),Y(superMask),Z(superMask),W(superMask));
        
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
        save(['~/data/tfr_motion_feature_data/band_power/band_power' folder_of_interest(30:31) '.mat'],'bandPower','-v7.3');
        save(['~/data/tfr_motion_feature_data/freq_entropy/freq_entropy' folder_of_interest(30:31) '.mat'],'freqEnt','-v7.3');
        save(['~/data/tfr_motion_feature_data/freq_pairs/freq_pairs' folder_of_interest(30:31) '.mat'],'freqPairs','-v7.3');
        save(['~/data/tfr_motion_feature_data/med_freq/med_freq' folder_of_interest(30:31) '.mat'],'medF','-v7.3');
        save(['~/data/tfr_motion_feature_data/sma/sma' folder_of_interest(30:31) '.mat'],'SMA','-v7.3');
        save(['~/data/tfr_motion_feature_data/wavelets/wavelets' folder_of_interest(30:31) '.mat'],'wavelets','-v7.3');
        disp(['Patient No. ' folder_of_interest(30:31) ' completed.']);
    else
        save(['../tfr_motion_feature_data/band_power/band_power' folder_of_interest(26:27) '.mat'],'bandPower','-v7.3');
        save(['../tfr_motion_feature_data/freq_entropy/freq_entropy' folder_of_interest(26:27) '.mat'],'freqEnt','-v7.3');
        save(['../tfr_motion_feature_data/freq_pairs/freq_pairs' folder_of_interest(26:27) '.mat'],'freqPairs','-v7.3');
        save(['../tfr_motion_feature_data/med_freq/med_freq' folder_of_interest(26:27) '.mat'],'medF','-v7.3');
        save(['../tfr_motion_feature_data/sma/sma' folder_of_interest(26:27) '.mat'],'SMA','-v7.3');
        save(['../tfr_motion_feature_data/wavelets/wavelets' folder_of_interest(26:27) '.mat'],'wavelets','-v7.3');
        disp(['Patient No. ' folder_of_interest(26:27) ' completed.']);
    end
    toc
end

license('inuse')