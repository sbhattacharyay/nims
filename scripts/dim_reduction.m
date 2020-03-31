clear; clc
% Get the folder to complete_sensor_data.mat (motion_feature_data)
path1 = uigetdir;
% Get the folder to clinical_extraction_output.mat (clinical_data)
path2 = uigetdir;


output_dir = 'PCA Output Rank';
mkdir(output_dir)

load([path1 filesep 'complete_sensor_data.mat'])
load([path2 filesep 'clinical_extraction_output.mat'])

tot_explained = cell(7,7);

% Sensor by sensor analysis
for row = 2:7
    for col = 1:7
        close('all'); clc; 
        feature = feature_names(col);
        curr_mat = sensors{row,col};
        
        % Ranks the matrix
        curr_mat = sort(curr_mat','descend')';
        % https://courses.engr.illinois.edu/bioe298b/sp2018/Lecture%20Examples/23%20PCA%20slides.pdf
        
        [coeff,score,latent,tsquared,explained] = pca(curr_mat);

        disp(explained(1:3))
        tot_explained{row,col} = explained(1:3);

        scatter3(score(:,1),score(:,2),score(:,3))
        axis equal
        xlabel('1st Principal Component')
        ylabel('2nd Principal Component')
        zlabel('3rd Principal Component')
        
        % Change name if output folder changed
        str = [pwd filesep output_dir filesep 'Sensor ', num2str(row), ' ', feature, '.fig'];
        saveas(gcf,join(str,''))
    end
end

save([pwd filesep output_dir filesep 'explained.mat'],'tot_explained')


% Aggregate all sensors together
band_power = zeros(size(sensors{1},1)*7,size(sensors{1},2));
freq_entropy = zeros(size(sensors{1},1)*7,size(sensors{1},2));
freq_pairs1 = zeros(size(sensors{1},1)*7,size(sensors{1},2));
freq_pairs2 = zeros(size(sensors{1},1)*7,size(sensors{1},2));
med_freq = zeros(size(sensors{1},1)*7,size(sensors{1},2));
sma = zeros(size(sensors{1},1)*7,size(sensors{1},2));
wavelets = zeros(size(sensors{1},1)*7,size(sensors{1},2));

counter = 1;
for row = 1:7
    band_power((counter-1)*62+1:counter*62,:) = sensors{row,1};
    freq_entropy((counter-1)*62+1:counter*62,:) = sensors{row,2};
    freq_pairs1((counter-1)*62+1:counter*62,:) = sensors{row,3};
    freq_pairs2((counter-1)*62+1:counter*62,:) = sensors{row,4};
    med_freq((counter-1)*62+1:counter*62,:) = sensors{row,5};
    sma((counter-1)*62+1:counter*62,:) = sensors{row,6};
    wavelets((counter-1)*62+1:counter*62,:) = sensors{row,7};
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

% Run the PCA
[coeff,score,latent,tsquared,explained] = pca(band_power);
disp(explained(1:3))
scatter3(score(:,1),score(:,2),score(:,3))
axis equal
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('3rd Principal Component')