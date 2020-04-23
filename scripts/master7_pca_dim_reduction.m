%% Authors: Matthew Wang, B.S. Candidate, Shubhayu Bhattacharyay, B.S. Candidate
% Department of Biomedical Engineering
% Department of Applied Mathematics and Statistics
% Whiting School of Engineering, Johns Hopkins University
% email address: shubhayu@jhu.edu
% March 2019; Last revision: 20-Feb-2020
%% ------------- BEGIN CODE --------------
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

% % Get the folder to complete_sensor_data.mat (motion_feature_data)
path1 = uigetdir;
% % Get the folder to patient_table.mat (clinical_data)
path2 = uigetdir;

% Create an output directory for the figures
output_dir = '../plots/Death_1yr';
mkdir(output_dir)

% Load in the imputed data and clinical table
load('../motion_feature_data/bed_corrected_imputed_complete_sensor_data.mat')
load('../clinical_data/patient_table.mat')

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

%%

% tot_explained = cell(7,7);
% % Sensor by sensor analysis
% for row = 2:7
%     for col = 1:7
%         close('all'); clc; 
%         feature = feature_names(col);
%         curr_mat = sensors{row,col};
%         
%         % Ranks the matrix
%         % curr_mat = sort(curr_mat','descend')';
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

%%

% Aggregate all sensors together
band_power = zeros(62,12959*6);
freq_entropy = zeros(62,12959*6);
freq_pairs1 = zeros(62,12959*6);
freq_pairs2 = zeros(62,12959*6);
med_freq = zeros(62,12959*6);
sma = zeros(62,12959*6);
wavelets = zeros(62,12959*6);

counter = 1;
for row = 1:6
    band_power(:,(counter-1)*12959+1:counter*12959) = bed_corrected_sensors{row,1};
    freq_entropy(:,(counter-1)*12959+1:counter*12959) = bed_corrected_sensors{row,2};
    freq_pairs1(:,(counter-1)*12959+1:counter*12959) = bed_corrected_sensors{row,3};
    freq_pairs2(:,(counter-1)*12959+1:counter*12959) = bed_corrected_sensors{row,4};
    med_freq(:,(counter-1)*12959+1:counter*12959) = bed_corrected_sensors{row,5};
    sma(:,(counter-1)*12959+1:counter*12959) = bed_corrected_sensors{row,6};
    wavelets(:,(counter-1)*12959+1:counter*12959) = bed_corrected_sensors{row,7};
    counter = counter + 1;
end

freq_pairs = [freq_pairs1 freq_pairs2];
% Rank sort
band_power = sort(band_power','descend')';
freq_entropy = sort(freq_entropy','descend')';
freq_pairs = sort(freq_pairs','descend')';
med_freq = sort(med_freq','descend')';
sma = sort(sma','descend')';
wavelets = sort(wavelets','descend')';

all = {band_power; freq_entropy; freq_pairs; med_freq; sma; wavelets};
feature_names = {'band_power'; 'freq_entropy'; 'freq_pairs'; 'med_freq'; 'sma'; 'wavelets'};
tot_explained = zeros(6,3);
% Run the PCA
for i = 1:6
    feature = feature_names{i};
    
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

%%

% % Aggregate just the wrist and ankle
% band_power_WA = zeros(62,12959*4);
% freq_entropy_WA = zeros(62,12959*4);
% freq_pairs1_WA = zeros(62,12959*4);
% freq_pairs2_WA = zeros(62,12959*4);
% med_freq_WA = zeros(62,12959*4);
% sma_WA = zeros(62,12959*4);
% wavelets_WA = zeros(62,12959*4);
% 
% counter = 1;
% for row = [4 7]
%     band_power_WA(:,(counter-1)*12959+1:counter*12959) = sensors{row,1};
%     freq_entropy_WA(:,(counter-1)*12959+1:counter*12959) = sensors{row,2};
%     freq_pairs1_WA(:,(counter-1)*12959+1:counter*12959) = sensors{row,3};
%     freq_pairs2_WA(:,(counter-1)*12959+1:counter*12959) = sensors{row,4};
%     med_freq_WA(:,(counter-1)*12959+1:counter*12959) = sensors{row,5};
%     sma_WA(:,(counter-1)*12959+1:counter*12959) = sensors{row,6};
%     wavelets_WA(:,(counter-1)*12959+1:counter*12959) = sensors{row,7};
%     counter = counter + 1;
% end
% % Rank Sort
% band_power_WA = sort(band_power_WA','descend')';
% freq_entropy_WA = sort(freq_entropy_WA','descend')';
% freq_pairs1_WA = sort(freq_pairs1_WA','descend')';
% freq_pairs2_WA = sort(freq_pairs2_WA','descend')';
% med_freq_WA = sort(med_freq_WA','descend')';
% sma_WA = sort(sma_WA','descend')';
% wavelets_WA = sort(wavelets_WA','descend')';
% 
% all_WA = {band_power_WA, freq_entropy_WA, freq_pairs1_WA, freq_pairs2_WA, med_freq_WA, sma_WA, wavelets_WA};
% 
% % Run RICA model, return 3 features
% Md1 = rica(all_WA{1},3);
