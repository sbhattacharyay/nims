%% Master Script 4: Collect Motion Features
% Decoding Quantitative Motor Features for Classification and Prediction
% in Severe Acquired Brain Injury
%
% Shubhayu Bhattacharyay, Matthew Wang
% Department of Biomedical Engineering
% Department of Applied Mathematics and Statistics
% Whiting School of Engineering, Johns Hopkins University
% email address: shubhayu@jhu.edu
%% ------------- BEGIN CODE --------------
% Feature collection and organization script
%% ALL:

addpath('functions/')

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

%% TOD:

cd ../tod_motion_feature_data/
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

%% TFR:

cd ../tfr_motion_feature_data/
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

cd ../scripts/

license('inuse')