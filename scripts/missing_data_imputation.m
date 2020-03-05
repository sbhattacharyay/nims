%% Authors: Shubhayu Bhattacharyay, B.S. Candidate, Matthew Wang, B.S. Candidate
% Department of Biomedical Engineering
% Department of Applied Mathematics and Statistics
% Whiting School of Engineering, Johns Hopkins University
% email address: shubhayu@jhu.edu
% March 2019; Last revision: 20-Feb-2020
%% ------------- BEGIN CODE --------------
%Set Directory and Procure all feature data
tic
cd ..
cd motion_feature_data
path = pwd;
directory = dir;
load('complete_sensor_data.mat');
toc

dim_of_sensors=size(sensors);
sensor_count = dim_of_sensors(1);
feature_count = dim_of_sensors(2);

studyPatientsPY = [2, 3,	4,	5,	6,	7,	8,	9,	10,	11,	12,	13, ...
    14,	15,	16,	17,	18,	19,	20,	21,	22,	23,	24,	26,	27,	28,	29,	30, ....
    31,	32,	33,	34,	35,	36,	37,	38,	39,	40,	41,	42,	43,	49,	51,	46, .....
    47,	48,	50,	52,	53,	54,	55,	56,	57,	59,	60,	61,	62,	63,	64,	65, ......
    67,	68];
%% Get Time Point Labels (in datenum format)
cd band_power
load('band_power02.mat')
t = bandPower{2, 1};
t_copy = t;
clear bandPower
cd ..
%% Identify rows that are completely missing data (and replace with NaN)
totallyMissingIdxs = [];
%contains indices of the totally missing time series as follows:
%Row_1 is the sensor number (1 through 7)
%Row_2 is the feature_count (1 through 7) that corresponds to "feature_names"
%Row_3 is the patient number that corresponds to the patient index (1
%through 62)

missingID=[0 NaN Inf];

for i=1:sensor_count
    for j=1:feature_count
        currMatrix = sensors{i,j};
        [ptRows,~] = size(currMatrix);
        for k = 1:ptRows
            if all(ismissing((currMatrix(k,:)),missingID))
                totallyMissingIdxs = [totallyMissingIdxs [i;j;k]];
                currMatrix(k,:) = NaN(1,length(currMatrix(k,:)));
            end
            sensors{i,j} = currMatrix;
        end
    end
end
%% Identifying distributions of the features:
% We first wish to identify whether the feature datapoints have any
% association with time of day. We begin my calculating the mean and
% variance of the feature columns by ignoring NaN values.

TOD_Means = cellfun(@(X) nanmean(X,1),sensors,'UniformOutput',false);
TOD_Std = cellfun(@(X) nanstd(X,1),sensors,'UniformOutput',false);
TOD_Diff = cellfun(@(X) abs(diff(X)),TOD_Means,'UniformOutput',false);
%% Mean by time of day figures:
%BandPower:
figure
subplot(7,1,1);
plot(t,TOD_Means{1,1});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Means{2,1});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Means{3,1});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Means{4,1});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Means{5,1});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Means{6,1});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Means{7,1});
datetick('x',15)
axis tight
suptitle('BandPower Means')

%Freq_Entropy

figure
subplot(7,1,1);
plot(t,TOD_Means{1,2});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Means{2,2});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Means{3,2});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Means{4,2});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Means{5,2});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Means{6,2});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Means{7,2});
datetick('x',15)
axis tight
suptitle('FreqEntropy Means')

%Freq_Pairs1

figure
subplot(7,1,1);
plot(t,TOD_Means{1,3});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Means{2,3});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Means{3,3});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Means{4,3});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Means{5,3});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Means{6,3});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Means{7,3});
datetick('x',15)
axis tight
suptitle('FreqPairs1 Means')

%Freq_Pairs2

figure
subplot(7,1,1);
plot(t,TOD_Means{1,4});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Means{2,4});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Means{3,4});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Means{4,4});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Means{5,4});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Means{6,4});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Means{7,4});
datetick('x',15)
axis tight
suptitle('FreqPairs2 Means')

%Med_Freq

figure
subplot(7,1,1);
plot(t,TOD_Means{1,5});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Means{2,5});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Means{3,5});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Means{4,5});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Means{5,5});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Means{6,5});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Means{7,5});
datetick('x',15)
axis tight
suptitle('MedFreq Means')

%SMA

figure
subplot(7,1,1);
plot(t,TOD_Means{1,6});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Means{2,6});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Means{3,6});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Means{4,6});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Means{5,6});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Means{6,6});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Means{7,6});
datetick('x',15)
axis tight
suptitle('SMA Means')

%Wavelets

figure
subplot(7,1,1);
plot(t,TOD_Means{1,7});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Means{2,7});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Means{3,7});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Means{4,7});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Means{5,7});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Means{6,7});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Means{7,7});
datetick('x',15)
axis tight
suptitle('Wavelet Means')
%% STD Plots
% STD by time of day figures:
% Mean by time of day figures:
%BandPower:
figure
subplot(7,1,1);
plot(t,TOD_Std{1,1});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Std{2,1});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Std{3,1});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Std{4,1});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Std{5,1});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Std{6,1});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Std{7,1});
datetick('x',15)
axis tight
suptitle('BandPower STD')

%Freq_Entropy

figure
subplot(7,1,1);
plot(t,TOD_Std{1,2});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Std{2,2});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Std{3,2});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Std{4,2});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Std{5,2});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Std{6,2});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Std{7,2});
datetick('x',15)
axis tight
suptitle('FreqEntropy STD')

%Freq_Pairs1

figure
subplot(7,1,1);
plot(t,TOD_Std{1,3});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Std{2,3});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Std{3,3});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Std{4,3});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Std{5,3});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Std{6,3});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Std{7,3});
datetick('x',15)
axis tight
suptitle('FreqPairs1 STD')

%Freq_Pairs2

figure
subplot(7,1,1);
plot(t,TOD_Std{1,4});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Std{2,4});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Std{3,4});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Std{4,4});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Std{5,4});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Std{6,4});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Std{7,4});
datetick('x',15)
axis tight
suptitle('FreqPairs2 STD')

%Med_Freq

figure
subplot(7,1,1);
plot(t,TOD_Std{1,5});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Std{2,5});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Std{3,5});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Std{4,5});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Std{5,5});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Std{6,5});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Std{7,5});
datetick('x',15)
axis tight
suptitle('MedFreq STD')

%SMA

figure
subplot(7,1,1);
plot(t,TOD_Std{1,6});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Std{2,6});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Std{3,6});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Std{4,6});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Std{5,6});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Std{6,6});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Std{7,6});
datetick('x',15)
axis tight
suptitle('SMA STD')

%Wavelets

figure
subplot(7,1,1);
plot(t,TOD_Std{1,7});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Std{2,7});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Std{3,7});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Std{4,7});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Std{5,7});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Std{6,7});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Std{7,7});
datetick('x',15)
axis tight
suptitle('Wavelet STD')
%% Diff Plots by time of day
t = t(1:end-1);
%BandPower:
figure
subplot(7,1,1);
plot(t,TOD_Diff{1,1});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Diff{2,1});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Diff{3,1});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Diff{4,1});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Diff{5,1});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Diff{6,1});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Diff{7,1});
datetick('x',15)
axis tight
suptitle('BandPower Diff')

%Freq_Entropy

figure
subplot(7,1,1);
plot(t,TOD_Diff{1,2});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Diff{2,2});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Diff{3,2});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Diff{4,2});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Diff{5,2});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Diff{6,2});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Diff{7,2});
datetick('x',15)
axis tight
suptitle('FreqEntropy Diff')

%Freq_Pairs1

figure
subplot(7,1,1);
plot(t,TOD_Diff{1,3});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Diff{2,3});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Diff{3,3});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Diff{4,3});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Diff{5,3});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Diff{6,3});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Diff{7,3});
datetick('x',15)
axis tight
suptitle('FreqPairs1 Diff')

%Freq_Pairs2

figure
subplot(7,1,1);
plot(t,TOD_Diff{1,4});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Diff{2,4});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Diff{3,4});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Diff{4,4});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Diff{5,4});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Diff{6,4});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Diff{7,4});
datetick('x',15)
axis tight
suptitle('FreqPairs2 Diff')

%Med_Freq

figure
subplot(7,1,1);
plot(t,TOD_Diff{1,5});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Diff{2,5});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Diff{3,5});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Diff{4,5});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Diff{5,5});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Diff{6,5});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Diff{7,5});
datetick('x',15)
axis tight
suptitle('MedFreq Diff')

%SMA

figure
subplot(7,1,1);
plot(t,TOD_Diff{1,6});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Diff{2,6});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Diff{3,6});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Diff{4,6});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Diff{5,6});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Diff{6,6});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Diff{7,6});
datetick('x',15)
axis tight
suptitle('SMA Diff')

%Wavelets

figure
subplot(7,1,1);
plot(t,TOD_Diff{1,7});
datetick('x',15)
axis tight

subplot(7,1,2);
plot(t,TOD_Diff{2,7});
datetick('x',15)
axis tight

subplot(7,1,3);
plot(t,TOD_Diff{3,7});
datetick('x',15)
axis tight

subplot(7,1,4);
plot(t,TOD_Diff{4,7});
datetick('x',15)
axis tight

subplot(7,1,5);
plot(t,TOD_Diff{5,7});
datetick('x',15)
axis tight

subplot(7,1,6);
plot(t,TOD_Diff{6,7});
datetick('x',15)
axis tight

subplot(7,1,7);
plot(t,TOD_Diff{7,7});
datetick('x',15)
axis tight
suptitle('Wavelet Diff')

t = t_copy;
%% Histogram of different sensor/feature combination means by time of day:
%BandPower:
figure
subplot(7,1,1);
histogram(TOD_Means{1,1});

axis tight

subplot(7,1,2);
histogram(TOD_Means{2,1});

axis tight

subplot(7,1,3);
histogram(TOD_Means{3,1});

axis tight

subplot(7,1,4);
histogram(TOD_Means{4,1});

axis tight

subplot(7,1,5);
histogram(TOD_Means{5,1});

axis tight

subplot(7,1,6);
histogram(TOD_Means{6,1});

axis tight

subplot(7,1,7);
histogram(TOD_Means{7,1});

axis tight
suptitle('BandPower Means Hist')

%Freq_Entropy

figure
subplot(7,1,1);
histogram(TOD_Means{1,2});

axis tight

subplot(7,1,2);
histogram(TOD_Means{2,2});

axis tight

subplot(7,1,3);
histogram(TOD_Means{3,2});

axis tight

subplot(7,1,4);
histogram(TOD_Means{4,2});

axis tight

subplot(7,1,5);
histogram(TOD_Means{5,2});

axis tight

subplot(7,1,6);
histogram(TOD_Means{6,2});

axis tight

subplot(7,1,7);
histogram(TOD_Means{7,2});

axis tight
suptitle('FreqEntropy Means Hist')

%Freq_Pairs1

figure
subplot(7,1,1);
histogram(TOD_Means{1,3});

axis tight

subplot(7,1,2);
histogram(TOD_Means{2,3});

axis tight

subplot(7,1,3);
histogram(TOD_Means{3,3});

axis tight

subplot(7,1,4);
histogram(TOD_Means{4,3});

axis tight

subplot(7,1,5);
histogram(TOD_Means{5,3});

axis tight

subplot(7,1,6);
histogram(TOD_Means{6,3});

axis tight

subplot(7,1,7);
histogram(TOD_Means{7,3});

axis tight
suptitle('FreqPairs1 Means Hist')

%Freq_Pairs2

figure
subplot(7,1,1);
histogram(TOD_Means{1,4});

axis tight

subplot(7,1,2);
histogram(TOD_Means{2,4});

axis tight

subplot(7,1,3);
histogram(TOD_Means{3,4});

axis tight

subplot(7,1,4);
histogram(TOD_Means{4,4});

axis tight

subplot(7,1,5);
histogram(TOD_Means{5,4});

axis tight

subplot(7,1,6);
histogram(TOD_Means{6,4});

axis tight

subplot(7,1,7);
histogram(TOD_Means{7,4});

axis tight
suptitle('FreqPairs2 Means Hist')

%Med_Freq

figure
subplot(7,1,1);
histogram(TOD_Means{1,5});

axis tight

subplot(7,1,2);
histogram(TOD_Means{2,5});

axis tight

subplot(7,1,3);
histogram(TOD_Means{3,5});

axis tight

subplot(7,1,4);
histogram(TOD_Means{4,5});

axis tight

subplot(7,1,5);
histogram(TOD_Means{5,5});

axis tight

subplot(7,1,6);
histogram(TOD_Means{6,5});

axis tight

subplot(7,1,7);
histogram(TOD_Means{7,5});

axis tight
suptitle('MedFreq Means Hist')

%SMA

figure
subplot(7,1,1);
histogram(TOD_Means{1,6});

axis tight

subplot(7,1,2);
histogram(TOD_Means{2,6});

axis tight

subplot(7,1,3);
histogram(TOD_Means{3,6});

axis tight

subplot(7,1,4);
histogram(TOD_Means{4,6});

axis tight

subplot(7,1,5);
histogram(TOD_Means{5,6});

axis tight

subplot(7,1,6);
histogram(TOD_Means{6,6});

axis tight

subplot(7,1,7);
histogram(TOD_Means{7,6});

axis tight
suptitle('SMA Means Hist')

%Wavelets

figure
subplot(7,1,1);
histogram(TOD_Means{1,7});

axis tight

subplot(7,1,2);
histogram(TOD_Means{2,7});

axis tight

subplot(7,1,3);
histogram(TOD_Means{3,7});

axis tight

subplot(7,1,4);
histogram(TOD_Means{4,7});

axis tight

subplot(7,1,5);
histogram(TOD_Means{5,7});

axis tight

subplot(7,1,6);
histogram(TOD_Means{6,7});

axis tight

subplot(7,1,7);
histogram(TOD_Means{7,7});

axis tight
suptitle('Wavelet Means Hist')
%% Finding no motion thresholds for each feature based on SMA:

% In order to impute the missing data, we first need to understand how
% "no motion" manifests in each of the given feature types.

% We assume that SMA is our ground truth in calculating No Motion
% thresholds.
SMA_threshold = 0.1;
SMA_nm_percentages = cellfun(@(x) noMotion_th(x,SMA_threshold),...
    sensors(:,6),'UniformOutput',false);

% Warning: Elapsed Time is 19.35 seconds

tic
% Based on the SMA_nm_percentages, we will derive thresholds that allow us
% to most closely match the percentages from SMA:

sma_idx = find(feature_names == "sma");
mf_idx = find(feature_names == "med_freq");

feature_thresholds = NaN(dim_of_sensors(2),1);
feature_thresholds(sma_idx) = SMA_threshold;

% Finding Band Power threshold:
j=1;
temp_feature_data = sensors(:,j);
min_norm = realmax;
opt_k = NaN;

for k = 0:0.001:0.15
    feat_percentages = cellfun(@(x) noMotion_th(x,k),...
        temp_feature_data,'UniformOutput',false);
    curr_norm=norm_calc(feat_percentages,SMA_nm_percentages);
    if curr_norm < min_norm
        opt_k = k;
        min_norm = curr_norm;
    end
end
feature_thresholds(j) = opt_k;

% Finding Freq Entropy threshold:
j=2;
temp_feature_data = sensors(:,j);
min_norm = realmax;
opt_k = NaN;

for k = 1.66:0.001:1.7
    feat_percentages = cellfun(@(x) noMotionFE_th(x,k),...
        temp_feature_data,'UniformOutput',false);
    curr_norm=norm_calc(feat_percentages,SMA_nm_percentages);
    if curr_norm < min_norm
        opt_k = k;
        min_norm = curr_norm;
    end
end
feature_thresholds(j) = opt_k;

% Finding Freq Pairs 1 threshold:
j=3;
temp_feature_data = sensors(:,j);
min_norm = realmax;
opt_k = NaN;

for k = 0.003:0.001:0.02
    feat_percentages = cellfun(@(x) noMotion_th(x,k),...
        temp_feature_data,'UniformOutput',false);
    curr_norm=norm_calc(feat_percentages,SMA_nm_percentages);
    if curr_norm < min_norm
        opt_k = k;
        min_norm = curr_norm;
    end
end
feature_thresholds(j) = opt_k;

% Finding Freq Pairs 2 threshold:
j=4;
temp_feature_data = sensors(:,j);
min_norm = realmax;
opt_k = NaN;

for k = 0.001:0.0001:0.006
    feat_percentages = cellfun(@(x) noMotion_th(x,k),...
        temp_feature_data,'UniformOutput',false);
    curr_norm=norm_calc(feat_percentages,SMA_nm_percentages);
    if curr_norm < min_norm
        opt_k = k;
        min_norm = curr_norm;
    end
end
feature_thresholds(j) = opt_k;

% Finding Med Freq threshold:
j=5;
temp_feature_data = sensors(:,j);
min_norm = realmax;
opt_k = NaN;

for k = 2.5:0.01:3.5
    feat_percentages = cellfun(@(x) noMotionFE_th(x,k),...
        temp_feature_data,'UniformOutput',false);
    curr_norm=norm_calc(feat_percentages,SMA_nm_percentages);
    if curr_norm < min_norm
        opt_k = k;
        min_norm = curr_norm;
    end
end
feature_thresholds(j) = opt_k;

% Finding Wavelet threshold:
j=7;
temp_feature_data = sensors(:,j);
min_norm = realmax;
opt_k = NaN;

for k = 1:0.1:15
    feat_percentages = cellfun(@(x) noMotion_th(x,k),...
        temp_feature_data,'UniformOutput',false);
    curr_norm=norm_calc(feat_percentages,SMA_nm_percentages);
    if curr_norm < min_norm
        opt_k = k;
        min_norm = curr_norm;
    end
end
feature_thresholds(j) = opt_k;
toc
%% Impute the totally-missing rows
%First, given our newly derived thresholds, we calculate the "no motion"
%for each of our patients per sensor per feature:

nm_percentages_for_features = {};

for j = 1:dim_of_sensors(2)
    temp_feature_data = sensors(:,j);
    k = feature_thresholds(j);
    if j == 2 || j == 5
        nm_percentages_for_features{j,1}=cellfun(@(x) noMotionFE_th(x,k), ...
            temp_feature_data,'UniformOutput',false);
    else
        nm_percentages_for_features{j,1}=cellfun(@(x) noMotion_th(x,k), ...
            temp_feature_data,'UniformOutput',false);
    end
end


% Then, we cycle through the completely missing time series
dim_TM = size(totallyMissingIdxs);

bp_nm_range = [0 feature_thresholds(1)];
fe_nm_range = [0 feature_thresholds(6)];
sma_nm_range = [0 feature_thresholds(6)];
sma_nm_range = [0 feature_thresholds(6)];
_nm_range = [0 feature_thresholds(6)];
sma_nm_range = [0 feature_thresholds(6)];
wav_nm_range = [0 feature_thresholds(7)];

for i = 1:dim_TM(2)
    sensIdx = totallyMissingIdxs(1,i);
    featIdx = totallyMissingIdxs(2,i);
    ptIdx = totallyMissingIdxs(3,i);
    draws=rand(1,length(t));
    curr_perc=nm_percentages_for_features{featIdx,1}{sensIdx,1};
    NM_imputes=draws<curr_perc(ptIdx);
    curr_mat=sensors{sensIdx,featIdx};
    curr_mat(ptIdx,NM_imputes)= ?????;
end
