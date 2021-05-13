%%%% Master Script 1: Complete Motion Feature Extraction %%%%
%
% Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
% University of Cambridge
% Johns Hopkins University
% email address: sb2406@cam.ac.uk
%
%%% Contents:
% I. Initialization
% II. Loop through study patients to begin feature extraction
% III. Filter Data (bworth 4-th order)
% IV. Feature Extraction
% V. Example of accelerometry for figure 1
% VI. Plot of frequncy response of Butterworth Filter

%% I. Initialization
% Set Directory and Procure File Names
tic
addpath('functions/')

d = dir('../accel_sensor_data/accel_*');
study_dirs = {d.name}';

% Extract temporal clinical metadata and get full list of study UPIs
patient_temporal_info = readtable('../clinical_data/patient_temporal_info.csv',...
    'TreatAsEmpty',{'.','NA'});
patient_temporal_info = sortrows(patient_temporal_info,'UPI');
study_upis = string(patient_temporal_info.UPI);

toc
%% II. Loop through study patients to begin feature extraction
%WARNING: Elapsed Time: ~111.286133 seconds per patient.
%WARNING: Elapsed Time: ~1104.244113 seconds for mega streams.

% Define sensor placement labels
placement_labels = ["Bed","LA","LE","LW","RA","RE","RW"]';

for curr_upi_idx = 1:size(study_upis,1)
    
    tic
    % Assign current UPI based on current loop index
    curr_upi = char(study_upis(curr_upi_idx));
    
    % Display patient intiation status
    disp(['Patient UPI: ',curr_upi,' initiated.']);
    
    % Find current accelerometry files with current UPI
    curr_upi_file_idx = find(cellfun(@(x) contains(x,curr_upi),study_dirs));
    
    % Load accelerometry of current UPI
    if length(curr_upi_file_idx) == 0 %current UPI does not exist in directory
        warning(['Patient UPI: ' curr_upi ' not found in accelerometry directory'])
        continue
    elseif length(curr_upi_file_idx) == 1 %current UPI exists as a single file
        load(['../accel_sensor_data/' study_dirs{curr_upi_file_idx}])
    elseif length(curr_upi_file_idx) == 7 %current UPI exists as 7 sensor-specific files
        for i = 1:length(curr_upi_file_idx)
            load(['../accel_sensor_data/' study_dirs{curr_upi_file_idx(i)}])
        end
        data = [C1;C2;C3;C4;C5;C6;C7]';
    else
        warning(['Patient UPI: ' curr_upi ' has an erroneous number of files (not 1 or 7)'])
        continue
    end
    
    % Get hours from ICU admission of start and end of recording for current UPI
    table_row_idx = find(string(patient_temporal_info.UPI) == curr_upi);
    curr_hours_from_ICU_adm_start = patient_temporal_info.HoursFromICUAdmissionAccelRecordingStart(table_row_idx);
    curr_hours_from_ICU_adm_end = patient_temporal_info.HoursFromICUAdmissionAccelRecordingEnd(table_row_idx);
    toc
    
    %% III. Filter Data (bworth 4-th order)
    tic
    
    % Set parameters for signal filtering
    fc = 0.2;
    fs = 10;
    [b,a] = butter(4,fc/(fs/2),'high');
    data_copy = cell(size(data));
    
    % Loop through different sensors
    for i = 1:length(data(1,:))
        % Arbitrarily assign start time to Jan 1, 1970 for data alignment
        start_date = datetime('01-Jan-1970');
        
        % Ad-hoc date correction for Patient 21, Left Elbow Sensor (#3)
        if curr_upi == 9kGOg64Y && i == 3
            start_date = start_date+1;
        end
        
        % Extract and filter accelerometry for current sensor
        curr = data(:,i);
        x = filter(b,a,curr{5});
        y = filter(b,a,curr{6});
        z = filter(b,a,curr{7});
        
        % Extract time of day stamps and convert to datetime variable types
        time = string(curr{12});
        time_cut = extractBefore(time,13);
        dt_info = strcat(string(start_date),{' '},time_cut);
        time_dn = datenum(dt_info,'dd-mmm-yyyy HH:MM:SS:FFF');
        
        % Use consecutive datenum differences to detect day changes from
        % the timestamps
        time_diff = diff(time_dn);
        [max_diffs,date_change_Idx] = findpeaks(-time_diff,'MinPeakHeight',.1,'MinPeakDistance',10);
        
        % Repair datenum values based on additional days
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
    %% IV. Feature Extraction
    %Warning! Run Time: >937.769720 seconds.
    tic
    
    % Initialize cells for each feature type
    SMA = {};
    HLF = {};
    MFR = {};
    BPW = {};
    FDE = {};
    WVL = {};
    
    % Set windowing start and end times. We pad both times inwards by 5
    % minutes to avoid possible data contamination at start and end times.
    windowingStart = datenum(dateshift(datetime(min(cellfun(@(x) x(1), ...
        data_copy(4,:))),'ConvertFrom', 'datenum'),'start','minute',....
        5));
    windowingEnd = datenum(dateshift(datetime(max(cellfun(@(x) x(end),...
        data_copy(4,:))),'ConvertFrom', 'datenum'),'start','minute',....
        -5));
    
    % Set windowing duration
    window_size_secs = 5; %in seconds
    window_size_hours = window_size_secs/3600;
    window_size_days = window_size_secs/86400;
    
    % Split windowing duration by window size
    time_split_days = windowingStart:window_size_days:windowingEnd;
    time_split_hours = curr_hours_from_ICU_adm_start:window_size_hours:(curr_hours_from_ICU_adm_start + (window_size_hours*(length(time_split_days)-2)));
    
    % Define bin counts
    bin_count = length(time_split_days)-1;
    RecordingIdx = 1:bin_count;
    
    % Loop through each sensor and split into defined bins and calculate
    % features from the bins
    for j = 1:size(data_copy,2)
        
        % create temporary cell to store split data
        split_Data = cell(bin_count,4);
        
        % extract filtered information of current sensor
        x_stream = data_copy{1,j};
        y_stream = data_copy{2,j};
        z_stream = data_copy{3,j};
        t_stream = data_copy{4,j};
        
        % split time stamps by defined bins
        [NCounts, ~, indices] = histcounts(t_stream,time_split_days);
        totalIdx = indices';
        
        % Define missing bins as those with less than 40/50 points
        nonMP = NCounts>40;
        availBins = find(nonMP);
        
        % Split all the information by the found bin indices
        split_Data(nonMP,1) = arrayfun(@(curr_bin) x_stream(totalIdx==curr_bin),availBins,'UniformOutput', false);
        split_Data(nonMP,2) = arrayfun(@(curr_bin) y_stream(totalIdx==curr_bin),availBins,'UniformOutput', false);
        split_Data(nonMP,3) = arrayfun(@(curr_bin) z_stream(totalIdx==curr_bin),availBins,'UniformOutput', false);
        split_Data(nonMP,4) = arrayfun(@(curr_bin) t_stream(totalIdx==curr_bin),availBins,'UniformOutput', false);
        
        % Find number of datapoints in each bin
        lens = cellfun(@(x) length(x),split_Data(:,1));
        
        % Define feature time points by removing the start time of the
        % first bin
        TimeOfDay = string(datestr(time_split_days(2:end)','HH:MM:SS'));
        
        % Define new cells for X, Y, Z, and T (time stamp)
        X = split_Data(:,1);
        Y = split_Data(:,2);
        Z = split_Data(:,3);
        T = split_Data(:,4);
        
        % Define masks for non empty indices
        maskx = RecordingIdx(cellfun(@(x) ~isempty(x), X));
        masky = RecordingIdx(cellfun(@(y) ~isempty(y), Y));
        maskz = RecordingIdx(cellfun(@(z) ~isempty(z), Z));
        maskt = RecordingIdx(cellfun(@(t) ~isempty(t), T));
        
        % Define super mask for when all variables are not missing
        superMask = intersect(maskx,masky);
        superMask = intersect(superMask,maskz);
        superMask = intersect(superMask,maskt);
        
        % Signal magnitude area (SMA)
        curr_SMA = NaN(bin_count,1);
        curr_SMA(superMask) = ((window_size_days)^-1).*cellfun(@(x,y,z,t) trapz(t,abs(x)+abs(y)+abs(z)), ...
            X(superMask),Y(superMask),Z(superMask),T(superMask));
        SMA = [SMA {curr_SMA;TimeOfDay;lens}];
        
        % Frequency component median pairs
        curr_HLF = NaN(bin_count,2);
        
        fc = 2.5;
        fs = 10;
        
        %%% Low-pass filtered component
        [d,c] = butter(4,fc/(fs/2),'low');
        curr_HLF(superMask,1) = cellfun(@(x,y,z,t) rssq([median(filter(d,c,x)),...
            median(filter(d,c,y)),median(filter(d,c,z))]),....
            X(superMask),Y(superMask),Z(superMask),T(superMask));
        
        %%% High-pass filtered component
        [b,a] = butter(4,fc/(fs/2),'high');
        curr_HLF(superMask,2) = cellfun(@(x,y,z,t) rssq([median(filter(b,a,x)),...
            median(filter(b,a,y)),median(filter(b,a,z))]),....
            X(superMask),Y(superMask),Z(superMask),T(superMask));
        
        HLF = [HLF {curr_HLF;TimeOfDay;lens}];
        
        %Median Freq
        curr_MFR = NaN(bin_count,1);
        curr_MFR(superMask)=cellfun(@(x,y,z,t) rssq([medfreq(x),medfreq(y),medfreq(z)]), ...
            X(superMask),Y(superMask),Z(superMask),T(superMask));
        
        MFR = [MFR {curr_MFR;TimeOfDay;lens}];
        
        % Band power between 0.3 and 3.5 Hz (BPW)
        curr_BPW = NaN(bin_count,1);
        curr_BPW(superMask)=cellfun(@(x,y,z,t) rssq([bandpower(x,fs,[0.3 3.5]), ...
            bandpower(y,fs,[0.3 3.5]),bandpower(z,fs,[0.3 3.5])]),....
            X(superMask),Y(superMask),Z(superMask),T(superMask));
        
        BPW = [BPW {curr_BPW;TimeOfDay;lens}];
        
        % Frequency-domain entropy (FDE)
        curr_FDE = NaN(bin_count,1);
        curr_FDE(superMask)=cellfun(@(x,y,z) rssq([pentropy(x,fs,'Instantaneous',false),...
            pentropy(y,fs,'Instantaneous',false),pentropy(z,fs,'Instantaneous',false)]), ....
            X(superMask),Y(superMask),Z(superMask));
        
        FDE = [FDE {curr_FDE;TimeOfDay;lens}];
        
        % Level 2 – 6 detail coefficients of the 5th-order Daubechies wavelet transform (WVL)
        curr_WVL = NaN(bin_count,1);
        curr_WVL(superMask)=cellfun(@(x,y,z) get_wavelets(x,y,z), ...
            X(superMask),Y(superMask),Z(superMask));
        
        WVL = [WVL {curr_WVL;TimeOfDay;lens}];
    end
    
    % Convert feature data into a tabular format
    %%% BPW
    BPW_table = array2table(cell2mat(BPW(1,:)),'VariableNames',placement_labels);
    BPW_table.UPI = repmat(string(curr_upi),height(BPW_table),1);
    BPW_table.RecordingIdx = RecordingIdx';
    BPW_table.HoursFromICUAdmission = time_split_hours';
    BPW_table.TimeOfDay = BPW{2,1};
    BPW_table.Feature = repmat("BPW",height(BPW_table),1);
    
    %%% FDE
    FDE_table = array2table(cell2mat(FDE(1,:)),'VariableNames',placement_labels);
    FDE_table.UPI = repmat(string(curr_upi),height(FDE_table),1);
    FDE_table.RecordingIdx = RecordingIdx';
    FDE_table.HoursFromICUAdmission = time_split_hours';
    FDE_table.TimeOfDay = FDE{2,1};
    FDE_table.Feature = repmat("FDE",height(FDE_table),1);
    
    %%% HLF_l
    HLF_l_table = array2table(cell2mat(cellfun(@(x) x(:,1),HLF(1,:),'UniformOutput',false)),'VariableNames',placement_labels);
    HLF_l_table.UPI = repmat(string(curr_upi),height(HLF_l_table),1);
    HLF_l_table.RecordingIdx = RecordingIdx';
    HLF_l_table.HoursFromICUAdmission = time_split_hours';
    HLF_l_table.TimeOfDay = HLF{2,1};
    HLF_l_table.Feature = repmat("HLF_l",height(HLF_l_table),1);
    
    %%% HLF_h
    HLF_h_table = array2table(cell2mat(cellfun(@(x) x(:,2),HLF(1,:),'UniformOutput',false)),'VariableNames',placement_labels);
    HLF_h_table.UPI = repmat(string(curr_upi),height(HLF_h_table),1);
    HLF_h_table.RecordingIdx = RecordingIdx';
    HLF_h_table.HoursFromICUAdmission = time_split_hours';
    HLF_h_table.TimeOfDay = HLF{2,1};
    HLF_h_table.Feature = repmat("HLF_h",height(HLF_h_table),1);
    
    %%% MFR
    MFR_table = array2table(cell2mat(MFR(1,:)),'VariableNames',placement_labels);
    MFR_table.UPI = repmat(string(curr_upi),height(MFR_table),1);
    MFR_table.RecordingIdx = RecordingIdx';
    MFR_table.HoursFromICUAdmission = time_split_hours';
    MFR_table.TimeOfDay = MFR{2,1};
    MFR_table.Feature = repmat("MFR",height(MFR_table),1);
    
    %%% SMA
    SMA_table = array2table(cell2mat(SMA(1,:)),'VariableNames',placement_labels);
    SMA_table.UPI = repmat(string(curr_upi),height(SMA_table),1);
    SMA_table.RecordingIdx = RecordingIdx';
    SMA_table.HoursFromICUAdmission = time_split_hours';
    SMA_table.TimeOfDay = SMA{2,1};
    SMA_table.Feature = repmat("SMA",height(SMA_table),1);
    
    %%% WVL
    WVL_table = array2table(cell2mat(WVL(1,:)),'VariableNames',placement_labels);
    WVL_table.UPI = repmat(string(curr_upi),height(WVL_table),1);
    WVL_table.RecordingIdx = RecordingIdx';
    WVL_table.HoursFromICUAdmission = time_split_hours';
    WVL_table.TimeOfDay = WVL{2,1};
    WVL_table.Feature = repmat("WVL",height(WVL_table),1);
    
    % Concatenate all feature tables, correct column order
    all_features_table = [BPW_table; FDE_table; HLF_l_table; HLF_h_table; MFR_table; SMA_table; WVL_table];
    all_features_table= all_features_table(:,[8:12 1:7]);
    
    % If any combination of features/sensors remains all 0 (exactly), replace with
    % all NA
    missingID=[0 NaN Inf];
    unique_features = unique(all_features_table.Feature);
    for curr_feature_idx=1:length(unique_features)
        curr_feature = unique_features(curr_feature_idx);
        curr_feature_rows = find(all_features_table.Feature == curr_feature);
        for curr_sensor_idx = 1:length(placement_labels)
            curr_sensor = placement_labels(curr_sensor_idx);
            curr_feature_sensor_values = all_features_table(curr_feature_rows,curr_sensor);
            if all(ismissing(curr_feature_sensor_values,missingID))
                all_features_table{curr_feature_rows,curr_sensor} = NaN(length(curr_feature_rows),1);
            end
        end
    end
    
    % Sort table appropriately
    all_features_table = sortrows(all_features_table,{'RecordingIdx','Feature'});
    
    % Save table into `features` directory
    writetable(all_features_table,['../features/features_' curr_upi '.csv']);
    
    % Display completion of current UPI
    disp(['Patient UPI: ' curr_upi ' completed.']);
    
    toc
end
%% V. Example of accelerometry for figure 1

% Load example patient and extract salient cells corresponding to
% accelerometry
load('../accel_sensor_data/accel_4oCkC1hc.mat');
data = data([5,6,7,12],:);

% Set parameters for signal filtering
fc = 0.2;
fs = 10;
[b,a] = butter(4,fc/(fs/2),'high');
data_copy = cell(size(data));

% Initialize new
sensorLabels = ["Bed","LA","LE","LW","RA","RE","RW"];

% UNFILTED ACCELEROMETRY:
% Iterate through different sensors and extract snippet for plotting
for i = 1:length(data(1,:))
    % Arbitrarily assign start time to Jan 1, 1970 for data alignment
    start_date = datetime('01-Jan-1970');
    
    % Extract salient time information for plot
    time = string(data{4,i});
    curr_sensor_idx = find((time >= '14:26:15') & (time <= '14:26:40'));
    
    % Extract indexed data of each axis
    x = data{1,i}(curr_sensor_idx);
    y = data{2,i}(curr_sensor_idx);
    z = data{3,i}(curr_sensor_idx);
    
    % Extract time of day stamps and convert to datetime variable types
    time_cut = extractBefore(time,13);
    dt_info = strcat(string(start_date),{' '},time_cut);
    time_dn = datenum(dt_info,'dd-mmm-yyyy HH:MM:SS:FFF');
    time_extract = time_dn(curr_sensor_idx);
    
    % Create current sensor figure
    curr_sensor_figure = figure;
    
    % Create axes for current sensor figure
    axes1 = axes('Parent',curr_sensor_figure);
    hold(axes1,'on');
    
    plot(time_extract,x,'LineWidth',1)
    plot(time_extract,y,'LineWidth',1)
    plot(time_extract,z,'LineWidth',1)
    
    rectangle('Position',[datenum('01-Jan-1970 14:26:34','dd-mmm-yyyy HH:MM:SS') -5.75 1/17280 .5],'FaceColor',[0.15,0.15,0.15])
    
    xlim([datenum('01-Jan-1970 14:26:15','dd-mmm-yyyy HH:MM:SS'),datenum('01-Jan-1970 14:26:40','dd-mmm-yyyy HH:MM:SS')]);
    ylim([-6 6]);
    datetick('x','HH:MM:SS')
    
    ylabel(sensorLabels(i),'FontWeight','bold','FontName','Arial');
    
    hold(axes1,'off');
    
    % Set the remaining axes properties
    set(axes1,'FontName','Arial','FontSize',5,'XAxisLocation','origin','XTick',...
        zeros(1,0),'XTickLabel','');
    set(gcf,'units','inches','position',[7.875,8.53,3.5,1]);
    
    % Filter the three axes of the current sensor
    filt_x = filter(b,a,data{1,i});
    filt_y = filter(b,a,data{2,i});
    filt_z = filter(b,a,data{3,i});
    
    % Save filtered and extracted data into the copy data cell
    data_copy{1,i} = filt_x(curr_sensor_idx);
    data_copy{2,i} = filt_y(curr_sensor_idx);
    data_copy{3,i} = filt_z(curr_sensor_idx);
    data_copy{4,i}= time_extract;
end

% FILTERED ACCELEROMETRY:
% Iterate through different sensors
for i = 1:length(data_copy(1,:))
    % Extract extracted filter accelerometry for current sensor
    x = data_copy{1,i};
    y = data_copy{2,i};
    z = data_copy{3,i};
    t = data_copy{4,i};
    
    % Create current sensor figure
    curr_sensor_figure = figure;
    
    % Create axes for current sensor figure
    axes1 = axes('Parent',curr_sensor_figure);
    hold(axes1,'on');
    
    plot(t,x,'LineWidth',1)
    plot(t,y,'LineWidth',1)
    plot(t,z,'LineWidth',1)
    
    rectangle('Position',[datenum('01-Jan-1970 14:26:34','dd-mmm-yyyy HH:MM:SS') -1.25 1/17280 .125],'FaceColor',[0.15,0.15,0.15])
    
    xlim([datenum('01-Jan-1970 14:26:15','dd-mmm-yyyy HH:MM:SS'),datenum('01-Jan-1970 14:26:40','dd-mmm-yyyy HH:MM:SS')]);
    ylim([-1.5 1.5]);
    datetick('x','HH:MM:SS')
    
    ylabel(sensorLabels(i),'FontWeight','bold','FontName','Arial');
    
    hold(axes1,'off');
    
    % Set the remaining axes properties
    set(axes1,'FontName','Arial','FontSize',5,'XAxisLocation','origin','XTick',...
        zeros(1,0),'XTickLabel','');
    set(gcf,'units','inches','position',[7.875,8.53,3.5,1]);
end


%% VI. Plot of frequncy response of Butterworth Filter
fc = 0.2;
fs = 10;
[b,a] = butter(4,fc/(fs/2),'high');
freqz(b,a)