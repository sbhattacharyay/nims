clear; clc

studyPatients = [2,3,4,5,6,7,8,9,10,11,12,13,14,15, ...
  16,17,18,19,20,21,22,23,24,26,27,28,29,30,31,32,33, ....
  34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50, .....
  51,52,53,54,55,58,59,60,61,62,63,64,66,67];
studyPatientsPY = [2,3,4,5,6,7,8,9,10,11,12,13, ...
  14,15,16,17,18,19,20,21,22,23,24,26,27,28,29,30, ....
  31,32,33,34,35,36,37,38,39,40,41,42,43,49,51,46, .....
  47,48,50,52,53,54,55,56,57,59,60,61,62,63,64,65, ......
  67,68];

% Get the folder to complete_sensor_data.mat (motion_feature_data)
path1 = uigetdir;
% Get the folder to patient_table.mat (clinical_data)
path2 = uigetdir;

% Create an output directory for the figures
output_dir = 'PCA Rank Death 1yr';
mkdir(output_dir)

% Load in the imputed data and clinical table
load([path1 filesep 'imputed_complete_sensor_data.mat'])
load([path2 filesep 'patient_table.mat'])

% Reorder to resolve the mismatch between accelerometer and clinical data
[~,sortOrder]=sort(studyPatientsPY);
patient_table = patient_table(sortOrder,:);
death = table2array(patient_table(:,12));
gose = table2array(patient_table(:,13));
death_1yr = table2array(patient_table(:,14));
gose_1yr = table2array(patient_table(:,15));

colors = zeros(62,3);
rowsToSetBlue = death_1yr == 0;
rowsToSetRed = death_1yr == 1;
colors(rowsToSetBlue,3) = 1;
colors(rowsToSetRed,1) = 1;

% tot_explained = cell(7,7);
% % Sensor by sensor analysis
% for row = 2:7
%     for col = 1:7
%         close('all'); clc; 
%         feature = feature_names(col);
%         curr_mat = sensors{row,col};
%         
%         % Ranks the matrix
%         curr_mat = sort(curr_mat','descend')';
%         
%         [coeff,score,latent,tsquared,explained] = pca(curr_mat);
% 
%         disp(explained(1:3))
%         tot_explained{row,col} = explained(1:3);
% 
%         scatter3(score(:,1),score(:,2),score(:,3),[],colors)
%         axis equal
%         xlabel('1st Principal Component')
%         ylabel('2nd Principal Component')
%         zlabel('3rd Principal Component')
%         
%         str = [pwd filesep output_dir filesep 'Sensor ', num2str(row), ' ', feature, '.fig'];
%         saveas(gcf,join(str,''))
%     end
% end
% 
% save([pwd filesep output_dir filesep 'explained.mat'],'tot_explained')


% Aggregate all sensors together
band_power = zeros(62,12959*7);
freq_entropy = zeros(62,12959*7);
freq_pairs1 = zeros(62,12959*7);
freq_pairs2 = zeros(62,12959*7);
med_freq = zeros(62,12959*7);
sma = zeros(62,12959*7);
wavelets = zeros(62,12959*7);

counter = 1;
for row = 1:7
    band_power(:,(counter-1)*12959+1:counter*12959) = sensors{row,1};
    freq_entropy(:,(counter-1)*12959+1:counter*12959) = sensors{row,2};
    freq_pairs1(:,(counter-1)*12959+1:counter*12959) = sensors{row,3};
    freq_pairs2(:,(counter-1)*12959+1:counter*12959) = sensors{row,4};
    med_freq(:,(counter-1)*12959+1:counter*12959) = sensors{row,5};
    sma(:,(counter-1)*12959+1:counter*12959) = sensors{row,6};
    wavelets(:,(counter-1)*12959+1:counter*12959) = sensors{row,7};
    counter = counter + 1;
end

% Rank sort
band_power = sort(band_power','descend')';
freq_entropy = sort(freq_entropy','descend')';
freq_pairs1 = sort(freq_pairs1','descend')';
freq_pairs2 = sort(freq_pairs2','descend')';
med_freq = sort(med_freq','descend')';
sma = sort(sma','descend')';
wavelets = sort(wavelets','descend')';

all = {band_power; freq_entropy; freq_pairs1; freq_pairs2; med_freq; sma; wavelets};
tot_explained = zeros(7,3);
% Run the PCA
for i = 1:7
    feature = feature_names(i);
    
    [coeff,score,latent,tsquared,explained] = pca(all{i});
    tot_explained(i,:) = explained(1:3);
    disp(explained(1:3))
    scatter3(score(:,1),score(:,2),score(:,3),[],colors)
    axis equal
    xlabel('1st Principal Component')
    ylabel('2nd Principal Component')
    zlabel('3rd Principal Component')
    
    str = [pwd filesep output_dir filesep feature '.fig'];
    saveas(gcf,join(str,''))
    close('all')
    clc
end
save([pwd filesep output_dir filesep 'explained.mat'],'tot_explained')

% Run RICA model, return 3 features
% Md1 = rica(band_power,3);
% Md2 = rica(freq_entropy,3);
% Md3 = rica(freq_pairs1,3);
% Md4 = rica(freq_pairs2,3);
% Md5 = rica(med_freq,3);
% Md6 = rica(sma,3);
% Md7 = rica(wavelets,3);