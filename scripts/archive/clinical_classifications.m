%% Author: Shubhayu Bhattacharyay
% Department of Biomedical Engineering
% Department of Applied Mathematics and Statistics
% Whiting School of Engineering, Johns Hopkins University
% email address: shubhayu@jhu.edu
% March 2019; Last revision: 13-March-2019
%% ------------- BEGIN CODE --------------
%% Determine fold-count
% k = length(folds(1,:));
%% APACHE II
clins = [ gcs_scores_en ];
gose_aucs = [];
figure
set(gca, 'FontName', 'Myriad Pro');
for i = 1:length(clins(1,:))
[X,Y,T,AUC]=perfcurve(favorable',clins(:,i)',1);
plot(X,Y,'linewidth',2);
hold on
gose_aucs = [gose_aucs AUC];
end
xq = linspace(0, 1, 101);
plot(xq,linspace(0,1,101),'k','LineStyle','--','LineWidth',1.5);
box on
xlabel('1 - Specificity','Fontsize',14);
ylabel('Sensitivity','Fontsize',14);
% legend({'AMPAC','APACHE','GCS','HLM','JFK'},'Location','northwest','NumColumns',2);
% legend boxoff    

title('Clinical Predictors of GOSE Outcomes','Fontsize',18);
pbaspect([1 1 1])

figure

clins = [apache_scores];
death_aucs = [];
set(gca, 'FontName', 'Myriad Pro');

for i = 1:length(clins(1,:))
[X,Y,T,AUC2]=perfcurve(death',clins(:,i)',1);
plot(X,Y,'linewidth',2);
hold on
death_aucs = [death_aucs AUC2];
end
xq = linspace(0, 1, 101);
plot(xq,linspace(0,1,101),'k','LineStyle','--','LineWidth',1.5);
box on
% legend({'AMPAC','APACHE','GCS','HLM','JFK'},'Location','northwest','NumColumns',2);
xlabel('1 - Specificity','Fontsize',14);
ylabel('Sensitivity','Fontsize',14);
% legend boxoff    
title('Clinical Predictors of Death','Fontsize',18);
pbaspect([1 1 1])
% yolo = table(["AMPAC","APACHE","GCS","HLM","JFK"]');
% yolo = addvars(yolo,death_aucs');
% yolo = addvars(yolo,gose_aucs');
%% GLM
glmAUCs = [];
glmXs = [];
glmYs = [];
xq = linspace(0, 1, 101);

for i = 1:k
    tset = logical(folds(:,i));
    dataTraining = dataset(tset,:);
    labelTraining = dataTraining(:,end);
    dataTraining = dataTraining(:,1:end-1);
    
    dataTesting = dataset(~tset,:);
    labelTesting = dataTesting(:,end);
    dataTesting = dataTesting(:,1:end-1);
    
    mdlGLM = fitglm(dataTraining,labelTraining,'Distribution','binomial');
    ypred = predict(mdlGLM,dataTesting);
    [X,Y,T,AUC]=perfcurve(labelTesting',ypred',1);
    
    glmAUCs = [glmAUCs AUC];
    
    uniq  = [true; diff(X) ~= 0];
    Value = interp1(X(uniq), Y(uniq), xq);
    Value(1) = 0;
    glmYs = [glmYs Value'];
end


meanGLM_AUC = mean(glmAUCs);
stdGLM_AUC = std(glmAUCs);

meanGLMY = mean(glmYs,2);

stdGLMY= std(glmYs,0,2);
upper = min(meanGLMY + stdGLMY,1);
lower = max(meanGLMY - stdGLMY,0);

doi = [xq, fliplr(xq)];
loi = [lower', fliplr(upper')];
fill(doi,loi,'b','LineStyle','none');
alpha(.1);
hold on
plot(xq,linspace(0,1,101),'k','LineStyle','--','LineWidth',1.5);
plot(xq,meanGLMY,'b','linewidth',2);
box on
xlabel('1 - Specificity');
ylabel('Sensitivity');
str = {strcat('AUC = ',num2str(meanGLM_AUC),'(',num2str(stdGLM_AUC),')')};
text(.1,.9,str)
% figure
% plot(base_X,meanGLMY,'b');
%% SVM
svmAUCs = [];
svmXs = [];
svmYs = [];

for i = 1:k
    tset = logical(folds(:,i));
    dataTraining = dataset(tset,:);
    labelTraining = dataTraining(:,end);
    dataTraining = dataTraining(:,1:end-1);
    
    dataTesting = dataset(~tset,:);
    labelTesting = dataTesting(:,end);
    dataTesting = dataTesting(:,1:end-1);
    mdlSVM = fitcsvm(dataTraining,labelTraining,'OptimizeHyperparameters','auto',...
        'HyperparameterOptimizationOptions',...
        struct('AcquisitionFunctionName','expected-improvement-plus'));
    [label, scores] = predict(mdlSVM,dataTesting);
    [X,Y,T,AUC]=perfcurve(labelTesting',scores(:,2)',1);
    svmAUCs = [svmAUCs AUC];
    
    uniq  = [true; diff(X) ~= 0];
    Value = interp1(X(uniq), Y(uniq), xq);
    Value(1) = 0;
    svmYs = [svmYs Value'];
end

meanSVM_AUC = mean(svmAUCs);
stdSVM_AUC = std(svmAUCs);

meanSVMY = mean(svmYs,2);

stdSVMY= std(svmYs,0,2);
upper = min(meanSVMY + stdSVMY,1);
lower = max(meanSVMY - stdSVMY,0);

doi = [xq, fliplr(xq)];
loi = [lower', fliplr(upper')];
fill(doi,loi,'r','LineStyle','none');
alpha(.1);
hold on
plot(xq,linspace(0,1,101),'k','LineStyle','--','LineWidth',1.5);
plot(xq,meanSVMY,'r','linewidth',2);
box on
xlabel('1 - Specificity');
ylabel('Sensitivity');
str = {strcat('AUC = ',num2str(meanSVM_AUC),'(',num2str(stdSVM_AUC),')')};
text(.1,.9,str)
%% k-NN
knnAUCs = [];
knnXs = [];
knnYs = [];

for i = 1:k
    tset = logical(folds(:,i));
    dataTraining = dataset(tset,:);
    labelTraining = dataTraining(:,end);
    dataTraining = dataTraining(:,1:end-1);
    
    dataTesting = dataset(~tset,:);
    labelTesting = dataTesting(:,end);
    dataTesting = dataTesting(:,1:end-1);
    mdlKNN = fitcknn(dataTraining,labelTraining,'OptimizeHyperparameters','auto',...
        'HyperparameterOptimizationOptions',...
        struct('AcquisitionFunctionName','expected-improvement-plus'));
    [label, scores] = predict(mdlKNN,dataTesting);
    [X,Y,T,AUC]=perfcurve(labelTesting',scores(:,2)',1);
    knnAUCs = [knnAUCs AUC];
    
    uniq  = [true; diff(X) ~= 0];
    Value = interp1(X(uniq), Y(uniq), xq);
    Value(1) = 0;
    knnYs = [knnYs Value'];
end
meanKNN_AUC = mean(knnAUCs);
stdKNN_AUC = std(knnAUCs);

meanKNNY = mean(knnYs,2);

stdKNNY= std(knnYs,0,2);
upper = min(meanKNNY + stdKNNY,1);
lower = max(meanKNNY - stdKNNY,0);

doi = [xq, fliplr(xq)];
loi = [lower', fliplr(upper')];
fill(doi,loi,'g','LineStyle','none');
alpha(.1);
hold on
plot(xq,linspace(0,1,101),'k','LineStyle','--','LineWidth',1.5);
plot(xq,meanKNNY,'g','linewidth',2);
box on
xlabel('1 - Specificity');
ylabel('Sensitivity');
str = {strcat('AUC = ',num2str(meanKNN_AUC),'(',num2str(stdKNN_AUC),')')};
text(.1,.9,str)
%% Random Forest
% for i = 1:k
%     tset = logical(folds(:,i));
%     dataTraining = dataset(tset,:);
%     labelTraining = dataTraining(:,end);
%     dataTraining = dataTraining(:,1:end-1);
%
%     dataTesting = dataset(~tset,:);
%     labelTesting = dataTesting(:,end);
%     dataTesting = dataTesting(:,1:end-1);
%     Mdl = fitcknn(dataTraining,labelTraining,'OptimizeHyperparameters','auto',...
%     'HyperparameterOptimizationOptions',...
%     struct('AcquisitionFunctionName','expected-improvement-plus'));
%     [label, scores] = predict(Mdl,dataTesting);
%     [X,Y,T,AUC]=perfcurve(labelTesting',scores(:,2)',1);
%     knnAUCs = [svmAUCs AUC];
% end
% meanKNN_AUC = mean(knnAUCs);
% stdKNN_AUC = std(knnAUCs);
%% ------------- END OF CODE --------------