%% Discover time windows
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

disp([num2str(length(studyDirs)) ' participants found'])

toc

%%

varNames = {'first_start','last_end','longest_hours','bed_start',...
    'bed_end','bed_hours','la_start','la_end','la_hours','le_start',....
    'le_end','le_hours','lw_start','lw_end','lw_hours','ra_start',.....
    'ra_end','ra_hours','re_start','re_end','re_hours','rw_start',......
    'rw_end','rw_hours'};

T = table('Size',[length(studyDirs) 24],'VariableTypes',...
    {'string','string','double','string','string','double','string',....
    'string','double','string','string','double','string','string',.....
    'double','string','string','double','string','string','double',......
    'string','string','double'},'VariableNames',varNames);


for i = 1:length(studyDirs)
    if marcc == true
        folder_of_interest = ['~/data/accel_sensor_data/' studyDirs{i}];
        disp(['Patient No. ' folder_of_interest(30:31) ' initiated.']);
    else
        folder_of_interest = ['../accel_sensor_data/',studyDirs{i}];
        disp(['Patient No. ' folder_of_interest(26:27) ' initiated.']);
    end
    
    tic
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
    
    tic
    for j = 1:length(data(1,:))
        idx_shift = j*3;
        curr = data(:,j);
        time = string(curr{12});
        time_cut = extractBefore(time,13);
        %time_dn = datenum(time_cut,'HH:MM:SS:FFF');
        T{i,idx_shift+1} = time_cut(1);
        T{i,idx_shift+2} = time_cut(end);
        T{i,idx_shift+3} = length(time)/10/60/60;
    end
    toc
    sortedStartTimes=sort(T{i,4:3:24});
    sortedEndTimes=sort(T{i,5:3:24});
    sortedSpans = sort(T{i,6:3:24});
    T{i,1} = sortedStartTimes(1);
    T{i,2} = sortedEndTimes(end);
    T{i,3} = sortedSpans(end);
    
    if marcc == true
        writetable(T,'~/data/accel_sensor_data/accel_time_info.xlsx')
        disp(['Patient No. ' folder_of_interest(30:31) ' complete.']);
    else
        writetable(T,'../accel_sensor_data/accel_time_info.xlsx')
        disp(['Patient No. ' folder_of_interest(26:27) ' complete.']);
    end
    
end

%% Reformat times

T=readtable('../accel_sensor_data/accel_time_info.xlsx');

T{:,1:3:24} = extractBefore(T{:,1:3:24},6);
T{:,2:3:24} = extractBefore(T{:,2:3:24},6);

if marcc == true
    writetable(T,'~/data/accel_sensor_data/accel_time_info.xlsx')
else
    writetable(T,'../accel_sensor_data/accel_time_info.xlsx')
end