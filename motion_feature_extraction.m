%% Authors: Shubhayu Bhattacharyay, B.S. Candidate, Matthew Wang, B.S. Candidate
% Department of Biomedical Engineering
% Department of Applied Mathematics and Statistics
% Whiting School of Engineering, Johns Hopkins University
% email address: shubhayu@jhu.edu
% March 2019; Last revision: 23-May-2019
%% ------------- BEGIN CODE --------------
%Set Directory and Procure File Names
tic
cd('/data/rsteven1/sensor_data')
%cd('D:/')
path = pwd;
directory = dir;

studyPatientsPY = [2, 3,	4,	5,	6,	7,	8,	9,	10,	11,	12,	13, ...
    14,	15,	16,	17,	18,	19,	20,	21,	22,	23,	24,	26,	27,	28,	29,	30, ....
    31,	32,	33,	34,	35,	36,	37,	38,	39,	40,	41,	42,	43,	49,	51,	46, .....
    47,	48,	50,	52,	53,	54,	55,	56,	57,	59,	60,	61,	62,	63,	64,	65, ......
    67,	68];

py_folders = {};
for i = 1:length(directory)
    for j = 1:length(studyPatientsPY)
        patientNum = num2str(studyPatientsPY(j),'%02.f');
        if startsWith(directory(i).name,'data') && (string(extractBetween(directory(i).name,5,6))==patientNum)
            py_folders{end+1} = directory(i).name;
        end
    end
end
folderNames = string(py_folders)';
toc
%% Arrange data by time and cut out extraneous times
%WARNING: Elapsed Time: ~111.286133 seconds.
%WARNING: Elapsed Time: ~1104.244113 seconds for mega streams.

for patIdx = 27:length(folderNames)
    tic
    folder_of_interest = ['/data/rsteven1/sensor_data/' folderNames{patIdx}];
    %folder_of_interest = ['D:/' folderNames{patIdx}];
    disp(['Patient No.' folder_of_interest(32:33) ' initiated.']);
    %disp(['Patient No.' folder_of_interest(8:9) ' initiated.']);
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
    
    dayLengths = [];
    data_copy = data;
    score_1 = cellfun(@(x) double(length(x)>=864000), data(12,:));
    score_2 = cellfun(@(x) double(length(x)>=(864000*2)), data(12,:));
    dayLengths = score_1 + score_2;
    
    for i = 1:length(data(1,:))
        curr = data(:,i);
        x = filter(b,a,curr{5});
        y = filter(b,a,curr{6});
        z = filter(b,a,curr{7});
        time = string(curr{12});
        time_cut = extractBefore(time,13);
        time_dn = datenum(time_cut,'HH:MM:SS:FFF');
        time_diff = (diff(time_dn));
        max_diffs = [];
        date_change_Idx = [];
        [max_diffs,date_change_Idx] = findpeaks(-time_diff,'MinPeakHeight',.1,'MinPeakDistance',10);
        date_change_Idx = date_change_Idx+1;
        time_dn_fixed = time_dn;
        if length(max_diffs) >= 1
            disp('holla')
            for j = length(date_change_Idx)
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
    
    windowingStart = datenum('18:00:00:000','HH:MM:SS:FFF');
    windowingEnd = datenum('12:00:00:000','HH:MM:SS:FFF')+1;
    windowSize = 5; %in seconds
    second = datenum('00:00:01:000','HH:MM:SS:FFF')-datenum('00:00:00:000','HH:MM:SS:FFF');
    window = windowSize.*second;
    
    timeSplit2 =windowingStart:window:windowingEnd;
    timeSplit1 = timeSplit2-1;
    timeSplit3 = timeSplit2+1;
    timeSplit4 = timeSplit2+2;
    
    binCount = length(timeSplit2)-1;
    binIdx = 1:binCount;
    
    for j = 1:length(data_copy)
        split_Data = cell(binCount,4);
        D = data_copy{4,j};
        x_stream = data_copy{1,j};
        y_stream = data_copy{2,j};
        z_stream = data_copy{3,j};
        
        [N2, ~, indices2] = histcounts(D,timeSplit2);
        [N1, ~, indices1] = histcounts(D,timeSplit1);
        [N3, ~, indices3] = histcounts(D,timeSplit3);
        [N4, ~, indices4] = histcounts(D,timeSplit4);
        
        totalIdx = [indices1';indices2'; indices3'; indices4'];
        NCounts = [N1;N2;N3;N4];
        
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
        times = [timeSplit1(2:end)'];
        
        totalIdxSet = 1:binCount;
        
        X = split_Data(:,1);
        Y = split_Data(:,2);
        Z = split_Data(:,3);
        W = split_Data(:,4);
        
        maskx = totalIdxSet(cellfun(@(x) ~isempty(x), X));
        masky = totalIdxSet(cellfun(@(y) ~isempty(y), Y));
        maskz = totalIdxSet(cellfun(@(z) ~isempty(z), Z));
        maskw = totalIdxSet(cellfun(@(w) ~isempty(w), W));
        
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
        [b,a] = butter(4,fc/(fs/2),'high');
        [d,c] = butter(4,fc/(fs/2),'low');
        
        curr_Med(superMask,1) = cellfun(@(x,y,z,w) rssq([median(filter(d,c,x)),...
            median(filter(d,c,y)),median(filter(d,c,z))]),....
            X(superMask),Y(superMask),Z(superMask),W(superMask));
        
        curr_Med(superMask,2) = cellfun(@(x,y,z,w) rssq([median(filter(b,a,x)),...
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
        cd('/data/rsteven1/sensor_data')
        %         cd('D:/')
        curr_wvlt = NaN(binCount,1);
        curr_wvlt(superMask)=cellfun(@(x,y,z) get_wavelets(x,y,z), ...
            X(superMask),Y(superMask),Z(superMask));
        
        wavelets = [wavelets {curr_wvlt;times;lens}];
    end
    save(['/data/rsteven1/sensor_features/band_power' folder_of_interest(32:33) '.mat'],'bandPower','-v7.3');
    save(['/data/rsteven1/sensor_features/freq_entropy' folder_of_interest(32:33) '.mat'],'freqEnt','-v7.3');
    save(['/data/rsteven1/sensor_features/freq_pairs' folder_of_interest(32:33) '.mat'],'freqPairs','-v7.3');
    save(['/data/rsteven1/sensor_features/med_freq' folder_of_interest(32:33) '.mat'],'medF','-v7.3');
    save(['/data/rsteven1/sensor_features/sma' folder_of_interest(32:33) '.mat'],'SMA','-v7.3');
    save(['/data/rsteven1/sensor_features/wavelets' folder_of_interest(32:33) '.mat'],'wavelets','-v7.3');
    %     save(['D:/sensor_features/band_power' folder_of_interest(8:9) '.mat'],'bandPower','-v7.3');
    %     save(['D:/sensor_features/freq_entropy' folder_of_interest(8:9) '.mat'],'freqEnt','-v7.3');
    %     save(['D:/sensor_features/freq_pairs' folder_of_interest(8:9) '.mat'],'freqPairs','-v7.3');
    %     save(['D:/sensor_features/med_freq' folder_of_interest(8:9) '.mat'],'medF','-v7.3');
    %     save(['D:/sensor_features/sma' folder_of_interest(8:9) '.mat'],'SMA','-v7.3');
    %     save(['D:/sensor_features/wavelets' folder_of_interest(8:9) '.mat'],'wavelets','-v7.3');
    toc
    disp(['Patient No.' folder_of_interest(32:33) ' completed.']);
    %disp(['Patient No.' folder_of_interest(8:9) ' completed.']);
end

function data = read_files()
% READ_FILES captures the data from the patient's seven sensors.
% read_files() returns an array of cell arrays with the data where each
% cell array corresponds to one sensor.
len = zeros(1,6);

fid = fopen('dataLA.txt');
C2 = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %s');
fclose(fid);
len(1) = length(C2{1});

fid = fopen('dataLE.txt');
C3 = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %s');
fclose(fid);
len(2) = length(C3{1});

fid = fopen('dataLW.txt');
C4 = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %s');
fclose(fid);
len(3) = length(C4{1});

fid = fopen('dataRA.txt');
C5 = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %s');
fclose(fid);
len(4) = length(C5{1});

fid = fopen('dataRE.txt');
C6 = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %s');
fclose(fid);
len(5) = length(C6{1});

fid = fopen('dataRW.txt');
C7 = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %s');
fclose(fid);
len(6) = length(C7{1});


fid = fopen('databed.txt');
C1 = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %s');
fclose(fid);

data = [C1;C2;C3;C4;C5;C6;C7]';
end
function start_time = get_start(data)
% Isolate the start time character arrays.
start_t_char = {};
i = 1;
for C = data
    start_t_char{i} = C{12}{1,1};
    i = i + 1;
end
times = string(start_t_char);
sortedTimes = sort(times);
start_time = sortedTimes(end);
end
function end_time = get_end(data)
end_t_char = {};
i = 1;
for C = data
    end_t_char{i} = C{12}{end,1};
    i = i + 1;
end
times = string(end_t_char);
sortedTimes = sort(times);
end_time = sortedTimes(1);
end
function waveOutput = get_wavelets(x,y,z)
[c1,l1] = wavedec(x,6,'db5');
[c2,l2] = wavedec(y,6,'db5');
[c3,l3] = wavedec(z,6,'db5');
xNorms = 0;
yNorms = 0;
zNorms = 0;

for p = 2:6
    xNorms=xNorms+norm(detcoef(c1,l1,p))^2;
    yNorms=yNorms+norm(detcoef(c2,l2,p))^2;
    zNorms=zNorms+norm(detcoef(c3,l3,p))^2;
end
waveOutput=rssq([xNorms,yNorms,zNorms]);
end