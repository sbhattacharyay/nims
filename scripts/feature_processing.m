%% Author: Shubhayu Bhattacharyay
% Department of Biomedical Engineering
% Department of Applied Mathematics and Statistics
% Whiting School of Engineering, Johns Hopkins University
% email address: shubhayu@jhu.edu
% March 2019; Last revision: 09-June-2019
%% ------------- BEGIN CODE --------------
%Set Directory and Procure File Names
tic
cd 'C:\Users\Shubhayu\Documents\Shubhayu_Projects\Johns Hopkins University\ICU Motion Study\MATLAB\motion_feature_data'
path = pwd;
directory = dir;
[a,b,c,d,e,f]=directory(3:end).name;
featureNames=[string(a);string(b);string(c);string(d);string(e);string(f)];

studyPatientsPY = [2	3	4	5	6	7	8	9	10	11	12	13 ...
    14	15	16	17	18	19	20	21	22	23	24	26	27	28	29	30 ....
    31	32	33	34	35	36	37	38	39	40	41	42	43	49	51	46 .....
    47	48	50	52	53	54	55	56	57	59	60	61	62	63	64	65 ......
    67	68]';

fileNames = arrayfun(@(x) [num2str(x,'%02.f') '.mat'],studyPatientsPY,...
    'UniformOutput',false);
cd(path)
%% Load All Features
allFeatures = {length(featureNames),2};
allReshapedFeatures = {length(featureNames),2};
allTimes = {length(featureNames),1};

for i = 1:length(featureNames)
    cd(strcat(path,'\',featureNames(i)))
    currFeatures = cell(length(studyPatientsPY),7);
    currLens = cell(length(studyPatientsPY),7);
    for j = 1:length(fileNames)
        true_name = strcat(featureNames(i),fileNames{j});
        T = load(true_name);
        feat=struct2cell(T);
        feat=feat{1};
        currFeatures(j,:) = feat(1,:);
        currLens(j,:) = feat(3,:);
    end
    
    currFeatures = currFeatures';
    currLens = currLens';
    
    allFeatures{i,1} = currFeatures;
    allFeatures{i,2} = currLens;
    
    reshapedCurr = reshape(currFeatures(2:end,:),(6.*length(studyPatientsPY)),1);
    reshapedLens = reshape(currLens(2:end,:),(6.*length(studyPatientsPY)),1);
    
    allReshapedFeatures{i,1} = reshapedCurr;
    allReshapedFeatures{i,2} = reshapedLens;
    
    currTime = feat{2,1};
    allTimes{i,1} = currTime;
end
%% Identify threshold for motion/no-motion for each feature space.
% Also idnetify distribution for each feature space

bpFeatures = allFeatures{1,1};
feFeatures = allFeatures{2,1};
fpFeatures = allFeatures{3,1};
mfFeatures = allFeatures{4,1};
smaFeatures= allFeatures{5,1};
wavFeatures= allFeatures{6,1};

combinedBP = [];
combinedFE = [];
combinedFP = [];
combinedMF = [];
combinedSMA = [];
combinedWAV = [];

dims = size(bpFeatures(2:end,:));
totalCellNo = dims(1)*dims(2);

for i = 1:totalCellNo
    combinedBP = [combinedBP; bpFeatures{i+62}];
    combinedFE = [combinedFE; feFeatures{i+62}];
    combinedFP = [combinedFP; fpFeatures{i+62}];
    combinedMF = [combinedMF; mfFeatures{i+62}];
    combinedSMA= [combinedSMA; smaFeatures{i+62}];
    combinedWAV =[combinedWAV; wavFeatures{i+62}];
end
%% Easy imputation
MissingIdxBP = cellfun(@(x) isnan(x),bpFeatures,'UniformOutput',false);
MissingIdxFE = cellfun(@(x) isnan(x),feFeatures,'UniformOutput',false);
MissingIdxFP = cellfun(@(x) isnan(x),fpFeatures,'UniformOutput',false);
MissingIdxMF = cellfun(@(x) isnan(x),mfFeatures,'UniformOutput',false);
MissingIdxSMA= cellfun(@(x) isnan(x),smaFeatures,'UniformOutput',false);
MissingIdxWAV= cellfun(@(x) isnan(x),wavFeatures,'UniformOutput',false);

cellfun(@(x,y) x(y),bpFeatures,MissingIdxBP,'UniformOutput',false);
%% Identifying Covariates for Motion and No-Motion Imputation
% First, I will validate in SMA
% Getting vector of $\pi$ values
shape_of_Features = size(corr_SMAFeatures);
flatDim = shape_of_Features(1)*shape_of_Features(2);

MissingIdx = cellfun(@(x) isnan(x),corr_SMAFeatures,'UniformOutput',false);

NoMotionPis = cellfun(@(x) sum(x<threshold)/length(x(~isnan(x))),corr_SMAFeatures);
NoMotionPis(NoMotionPis == 1) = 0.999;
Y_pi = log(NoMotionPis./(1-NoMotionPis));
Y_pi = reshape(Y_pi,[],1);
X_pi = imputeDataset;
X_pi = repmat(X_pi,6,1);
mdl_logitPi = stepwiselm(X_pi,Y_pi,'interactions','CategoricalVars',[2],'VarNames',[clinNames 'logit(pi)']);
gamma_est = mdl_logitPi.Coefficients.Estimate;
Y_pi_hat = predict(mdl_logitPi,X_pi);
pi_hat = sigmf(Y_pi_hat,[1 0]);

% Now we look at $\lambda$
MotionLambdas = cellfun(@(x) nanmean(x(x>=threshold)),corr_SMAFeatures);
Y_L = log(MotionLambdas);
Y_L = reshape(Y_L,[],1);
X_L = imputeDataset;
X_L = repmat(X_L,6,1);
mdl_logLam = stepwiselm(X_L,Y_L,'interactions','CategoricalVars',[2],'VarNames',[clinNames 'log(lambda)']);
beta_est = mdl_logLam.Coefficients.Estimate;
Y_L_hat = predict(mdl_logLam,X_L);
L_hat = exp(Y_L_hat);

imputationRun = corr_SMAFeatures;

for i = 1:flatDim
    currSet = imputationRun{i};
    placeholder = nanmean(currSet);
    currMissIdx = find(MissingIdx{i});
    currPi_Hat = pi_hat(i);
    currL_Hat = L_hat(i);
    no_motion_selection = rand(length(currMissIdx),1) < currPi_Hat;
    poiss_selection = currMissIdx(~no_motion_selection);
    currSet(currMissIdx(no_motion_selection)) = (threshold.*rand(sum(no_motion_selection),1));
    for j = 1:length(poiss_selection)
        currPoissIdx = poiss_selection(j);
        currPoissIdx-20:currPoissIdx+20;
    end
end
%% Identifying Distribution for No-Motion and Motion

% Validate distribution in each frequency space for each patient?

% Distinguish between no motion and motion?
smaFeatures = allFeatures{5,1};
currSet = imputationRun{2};
currSetLuv = currSet(currSet>.1) - .1;
currSetNoLuv = currSet(currSet<=.1);

%% Bed Data Correction Algorithm:
% 1. In SMA domain, locate periods of time when bed contributes significant
% motion
% 2. Save those indices as markers for bed motion
% 3. Analyze SMA profiles directly preceeding those markers. If any motion
% sensor's SMA profile is significantly active prior to bed motion, then we
% 4. If no motion directly preceeds bed motion in those sensors, subtract
% bed SMA from those points.

time = allTimes{5,1};
smaFeatures = allFeatures{5,1};
bedSMA = smaFeatures(1,:);
bedActivationIndices = cellfun(@(x) find(diff(x) >= threshold),bedSMA,'UniformOutput',false);
bedActivationTimes = cellfun(@(x) time(x+1),bedActivationIndices,'UniformOutput',false);

corr_SMAFeatures = cell(length(studyPatientsPY),6);
for i = 1:length(studyPatientsPY)
    currSMA = smaFeatures(2:end,i);
    bedSMA = smaFeatures{1,i};
    currIdx=bedActivationIndices{i};
    for j =1:length(currSMA)
        tempSMA = currSMA{j};
        tempSMA(currIdx) = max(tempSMA(currIdx)-bedSMA(currIdx),threshold*rand);
        corr_SMAFeatures{i,j} = tempSMA;
    end
end
%% No Motion and Missing Data Characterization
threshold = 0.1;  %in g

reshapedSMA = allReshapedFeatures{5,1};
reshapedSMALens = allReshapedFeatures{5,2};
time = allTimes{5,1};

matSMA = [];
matLens = [];

for k = 1:(6.*length(studyPatientsPY))
    matSMA = [matSMA; reshapedSMA{k}'];
    matLens = [matLens; reshapedSMALens{k}'];
end

missingCurve = [];
noMotionCurve = [];

for j = 1:length(time)
    noMotionCurve = [noMotionCurve sum(matSMA(:,j)<threshold)./(6.*length(studyPatientsPY))];
    missingCurve = [missingCurve sum(matLens(:,j)<40)./(6.*length(studyPatientsPY))];
end

noMotionCurve = noMotionCurve+missingCurve;
overallMissPerc = trapz(time,missingCurve);
overallNMPerc =trapz(time,noMotionCurve)-trapz(time,missingCurve);
%% Overall Missingness and No Motion Profile

figure
area(time,missingCurve.*100,'EdgeColor','none','FaceColor',[0.9 0.9 0.9],'FaceAlpha',.8);
hold on
plot(time,missingCurve.*100,'r','LineWidth',2);
plot(time,noMotionCurve.*100,'b');
ylim([0 100]);
datetick('x','HH PM');
ax = gca;
ax.FontSize = 14;
xlabel('Time of Day','FontWeight','bold','FontSize',20);
ylabel('Percentage of Sensor Information ({\itn} = 372)','FontWeight','bold','FontSize',20);
box on;
ytickformat('percentage');
set(gca, 'FontName', 'Myriad Pro');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 10], 'PaperUnits', 'Inches', 'PaperSize', [10, 10])
%% Missing Profiles by Time of Day
figure
patLens = [];

for i = 1:length(studyPatientsPY)
    startIdx = 6.*(i-1)+1;
    subLens = matLens(startIdx:startIdx+5,:)<40;
    patLens = [patLens; sum(subLens)];
end
x = [time(1) time(end)];
y = [1 62];
colormap(flipud(pink))

image(x,y,patLens,'CDataMapping','scaled')
ax = gca;
ax.FontSize = 14;
datetick('x','HH PM');
set(gca, 'FontName', 'Myriad Pro');
h=colorbar;
ylabel(h, 'No. of Sensors Missing','FontSize',14)
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 10], 'PaperUnits', 'Inches', 'PaperSize', [10, 10])
xlabel('Time of Day','FontWeight','bold','FontSize',20);
ylabel('Patient No.','FontWeight','bold','FontSize',20);

%% Feature Visualization Script (SMA)

currFeatures = allFeatures{5,1}';
currLens = allFeatures{5,1}';

favSMAFeatures = currFeatures(favorable,:);
unfavSMAFeatures = currFeatures(~favorable,:);

survSMAFeatures= currFeatures(~death,:);
deathSMAFeatures= currFeatures(death,:);

time = allTimes{5,1};

sensorLabels = ["Bed","LA","LE","LW","RA","RE","RW"];

% Specify Time Window to Plot (in 'HH:MM') !!!Remember to add one for day
% times!!!
startPlotTime = '22:00';
endPlotTime = '22:30';

startPlotTime = datenum(startPlotTime,'HH:MM');
endPlotTime = datenum(endPlotTime,'HH:MM');

% Plotting Favorables
figure
for k =1:7
    flat = cell2mat(favSMAFeatures(:,k));
    refData = reshape(flat,12958,sum(favorable));
    meanSignal = nanmean(refData,2);
    stdSignal = nanstd(refData,0,2);
    upper = (meanSignal + stdSignal);
    lower = max(meanSignal - stdSignal,0);
    
    subplot(7,1,k);
    doi = [time', fliplr(time')];
    loi = [lower', fliplr(upper')];
    fill(doi,loi,'b','LineStyle','none');
    alpha(.1);
    hold on
    plot(time,meanSignal,'b','linewidth',.2);
    datetick('x','keeplimits');
    ylim([0 0.2]);
    xlim([startPlotTime endPlotTime]);
    ylabel(sensorLabels(k));
    box on
    set(gca, 'FontName', 'Myriad Pro');
end

% Plotting Unfavorables
figure

for k =1:7
    flat = cell2mat(unfavSMAFeatures(:,k));
    refData = reshape(flat,12958,length(studyPatientsPY)-sum(favorable));
    meanSignal = nanmean(refData,2);
    stdSignal = nanstd(refData,0,2);
    upper = (meanSignal + stdSignal);
    lower = max(meanSignal - stdSignal,0);
    
    subplot(7,1,k);
    doi = [time', fliplr(time')];
    loi = [lower', fliplr(upper')];
    fill(doi,loi,'r','LineStyle','none');
    alpha(.1);
    hold on
    plot(time,meanSignal,'r','linewidth',.2);
    datetick('x','keeplimits');
    ylim([0 0.2]);
    xlim([startPlotTime endPlotTime]);
    ylabel(sensorLabels(k));
    box on
    set(gca, 'FontName', 'Myriad Pro');
end
% %Plotting Survivals
figure
for k =1:7
    flat = cell2mat(survSMAFeatures(:,k));
    refData = reshape(flat,12958,length(studyPatientsPY)-sum(death));
    meanSignal = nanmean(refData,2);
    stdSignal = nanstd(refData,0,2);
    upper = (meanSignal + stdSignal);
    lower = max(meanSignal - stdSignal,0);
    
    subplot(7,1,k);
    doi = [time', fliplr(time')];
    loi = [lower', fliplr(upper')];
    fill(doi,loi,'b','LineStyle','none');
    alpha(.1);
    hold on
    plot(time,meanSignal,'b','linewidth',.2);
    datetick('x','keeplimits');
    ylim([0 0.2]);
    xlim([startPlotTime endPlotTime]);
    ylabel(sensorLabels(k));
    box on
    set(gca, 'FontName', 'Myriad Pro');
end
% %Plotting Deaths
figure
for k =1:7
    flat = cell2mat(deathSMAFeatures(:,k));
    refData = reshape(flat,12958,sum(death));
    meanSignal = nanmean(refData,2);
    stdSignal = nanstd(refData,0,2);
    upper = (meanSignal + stdSignal);
    lower = max(meanSignal - stdSignal,0);
    
    subplot(7,1,k);
    doi = [time', fliplr(time')];
    loi = [lower', fliplr(upper')];
    fill(doi,loi,'r','LineStyle','none');
    alpha(.1);
    hold on
    plot(time,meanSignal,'r','linewidth',.2);
    datetick('x','keeplimits');
    ylim([0 0.2]);
    xlim([startPlotTime endPlotTime]);
    ylabel(sensorLabels(k));
    box on
    set(gca, 'FontName', 'Myriad Pro');
end
%% Functions
function imputedSignal = accelmissing(data,m,cvNM,cvMotion,threshold)

missing = cellfun(@(y) isnan(y), data, 'UniformOutput', false);
passiveMotion = cellfun(@(x) x<threshold, data,'UniformOutput', false);

activeMotion = cell(size(data));
lambda_hats = zeros(length(studyPatientsPY),6);
for i = 1:length(studyPatientsPY)
    for j = 1:6
        currentDist = data{i,j};
        lambda_hats(i,j)=poissfit(sort((currentDist(currentDist>threshold))));
    end
end
lambda_hats(isnan(lambda_hats)) = threshold.*rand(length(lambda_hats(isnan(lambda_hats))),1);

Betas=pinv([ones(length(cvMotion),1) cvMotion])*log(lambda_hats);

pi_hats = cellfun(@(x) nansum(x)/numel(x(~isnan(x))),passiveMotion);

Gammas=pinv([ones(length(cvNM),1) cvNM])*(log(pi_hats) - log(1-pi_hats));

B_hats = [];

% tic
% Gammas = [];
% for i = 1:length(data(1,:))
%     sensorSubset1 = cell2mat(cellfun(@transpose, passiveMotion(:,i), 'uniform', 0));
%     for j = 1:length(sensorSubset1(1,:))
%         Y_gamma=sensorSubset1(:,j);
%         gamma_i = glmfit(cvNM,Y_gamma,'binomial','link','logit');
%         Gammas = [Gammas gamma_i];
%     end
% end
% toc
%
% tempGamma = cell(length(data(1,:)),1);
% for i = 1:length(data(1,:))
%     range=(1+(i-1).*length(Gammas)/6):(i.*length(Gammas)/6);
%     tempGamma{i}=Gammas(:,range);
% end
% Gammas = tempGamma;
%
% for i = 1:length(data(1,:))
%     beta = Betas(:,i);
%     cut = Gammas{i,1};
%     gamma = cut(:,1);
%     B_hat = [gamma, beta]';
%     V=cov(B_hat);
%     z=randn(2,8);
%     B_dot=B_hat+z*V^(-1/2);
%     gamma_dot =B_dot(1,:);
%     beta_dot =B_dot(2,:);
%     for j =length(data(:,1))
%         currMiss = missing(j,i);
%     end
% end

end