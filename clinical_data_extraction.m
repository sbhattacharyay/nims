%% 1). Clinical Data Extraction
% Decoding Quantitative Motor Features for Classification and Prediction
% in Severe Acquired Brain Injury
%
% Shubhayu Bhattacharyay
% Department of Biomedical Engineering
% Department of Applied Mathematics and Statistics
% Whiting School of Engineering, Johns Hopkins University
% email address: shubhayu@jhu.edu
% March 2019; Last revision: 29-January-2020
%% ------------- BEGIN CODE --------------
% Load data from most recent patient clinical data spreadsheet
cd 'C:\Users\Shubhayu\OneDrive - Johns Hopkins\BIMS Project\clinical_data'
patientData = readtable('Final_QMACP.xlsx','TreatAsEmpty',{'.','NA'});

% Assign variable names to patient datasheet
patientData.Properties.VariableNames([64 66 17 18 20 21 22 23 1 2 7 25 ...
    35 62 53 34 3 43 44 76]) = {'gose' 'death' 'hlm_en' 'hlm_dis' ....
    'ampac_en_r' 'ampac_en_t' 'ampac_dis_r' 'ampac_dis_t' .....
    'pNum' 'gender' 'los' 'apache' 'gcs_en' 'jfk_dis' 'jfk_en' ......
    'gcs_m_en' 'age' 'gcs_m_dis' 'gcs_dis' 'gose_12mo'};

% Remove extraneous cells in spreadsheet
excluded = isnan(patientData.pNum);
patientData(excluded,:) = [];

% Indices of patients to use in study. To see exclusion criteria, please
% refer to Materials and Methods section
studyPatients = [2	3	4	5	6	7	8	9	10	11	12	13	14	15 ...
    16	17	18	19	20	21	22	23	24	26	27	28	29	30	31	32	33 ....
    34	35	36	37	38	39	40	41	42	43	44	45	46	47	48	49	50 .....
    51	52	53	54	55	58	59	60	61	62	63	64	66	67];

%% Extract Outcome Information from Patient Dataset

fav_threshold = 4; % GOSE >= fav_threshold is favorable

% Extract outcome scores (GOSE) at discharge
gose_scores = patientData.gose;
gose_scores(3) = {'3'}; % Fix single ambiguous GOSE score in dataset
gose_scores = str2double(gose_scores(studyPatients));

% Extract outcome scores (GOSE) at 12 months (for patients for which data
% is available).
yr_outcome_code = str2num(cell2mat(patientData.Outcome12Months));
yr_outcome_code = yr_outcome_code(studyPatients);

gose_12months = patientData.gose_12mo;
gose_12months = gose_12months(studyPatients);

noYearOutcomeAvailable = isnan(gose_12months);

% Define functional outcomes (favorable outcomes as GOSE >= 4)
dc_outcomes = double(gose_scores>= fav_threshold); %discharge functional outcomes (F vs. UF)
yr_outcomes = double(gose_12months>= fav_threshold);
yr_outcomes(noYearOutcomeAvailable) = NaN;

% Mortality outcomes at discharge and 1 year
dc_death = str2double(patientData.death(studyPatients));
yr_death = double(gose_12months==1);
yr_death(noYearOutcomeAvailable) = NaN;

%% Clinical Data Preprocessing and Imputation

% Apply LLF-optimized box-cox transform to quantitative clinical variables
% and normalize dataset to enable hypothesis testing:

% Age:
[age,lam_age]=boxcox((patientData.age(studyPatients)));
[age,mu_age,sig_age]=zscore(age);

% GCS at enrollment:
[gcs_scores_en,lam_gcs_en]=boxcox(str2double(patientData.gcs_en(studyPatients)));
[gcs_scores_en,mu_gcs_en,sig_gcs_en]=zscore(gcs_scores_en);

% APACHE at enrollment:
[apache_scores,lam_apache]=boxcox(str2double(patientData.apache(studyPatients)));
[apache_scores,mu_apache,sig_apache]=zscore(apache_scores);

% GCS at discharge:
[gcs_scores_dis,lam_gcs_dis]=boxcox(str2double(patientData.gcs_dis(studyPatients)));
[gcs_scores_dis,mu_gcs_dis,sig_gcs_dis]=zscore(gcs_scores_dis);

% Gender (Males are 1, Females are -1):
females = categorical(patientData.gender(studyPatients))~= 'M';
gender = (-1).^females;

% Impute missing clinical variables with Weighted k-NN (code below)
dataset = [age gender apache_scores gcs_scores_en gcs_scores_dis];
imputeDataset = wKNN_impute(dataset);

% % AMPAC at enrollment:
% [ampac_en_r,lam_ampac]=boxcox(cell2mat(patientData.ampac_en_r(studyPatients)));
% [ampac_en_r,mu_ampac,sig_ampac]=zscore(ampac_en_r);
%
% % HLM at enrollment:
% [hlm_en,lam_hlm_en]=boxcox(patientData.hlm_en(studyPatients));
% [hlm_en,mu_hlm_en,sig_hlm_en]=zscore(hlm_en);
%
% % JFK at enrollment:
% [jfk_en,lam_jfk_en]=boxcox(patientData.jfk_en(studyPatients));
% [jfk_en,mu_jfk_en,sig_jfk_en]=zscore(jfk_en);
%% Prepare predictor matrices:

% Penultimate column is death and final column is GOSE
dc_dataset = [imputeDataset(:,1:(end-1)) dc_death dc_outcomes];

% Remove patients who lack 1 year outcomes to form 1-year outcomes dataset
yr_dataset = [imputeDataset(~noYearOutcomeAvailable,:) ...
    yr_death(~noYearOutcomeAvailable) ....
    yr_outcomes(~noYearOutcomeAvailable)];

dc_dataset_labels=["Age" "Gender" "APACHE" "GCS_{en}" "Death" "GOSE"];
yr_dataset_labels=["Age" "Gender" "APACHE" "GCS_{en}" "GCS_{dis}" ...
    "Death" "GOSE"];

boxcox_lambdas = [lam_age,lam_apache,lam_gcs_en,lam_gcs_dis];
zscore_mus = [mu_age,mu_apache,mu_gcs_en,mu_gcs_dis];
zscore_sigs = [sig_age,sig_apache,sig_gcs_en,sig_gcs_dis];

clearvars -except dc_dataset yr_dataset dc_dataset_labels ...
    yr_dataset_labels boxcox_lambdas zscore_mus zscore_sigs ....
    imputeDataset studyPatients

save('clinical_extraction_output.mat')
%% Visualize distribution of transformed (Box-Cox) clinical variables
figure
for i =1:5
    subplot(2,3,i);
    histogram(imputeDataset(:,i),30);
    title(yr_dataset_labels(i));
end

%% Visualize outcome distributions
figure
subplot(2,2,1);
pie(categorical(dc_dataset(:,(end-1))))
title('Died during hospital stay')

subplot(2,2,2);
pie(categorical(dc_dataset(:,(end))))
title('GOSE >= 4 at discharge')

subplot(2,2,3);
pie(categorical(yr_dataset(:,(end-1))))
title('Died within 1 year of discharge')

subplot(2,2,4);
pie(categorical(yr_dataset(:,(end))))
title('GOSE >= 4 within 1 year of discharge')
%% Weighted KNN Function
function filledMatrix = wKNN_impute(data)
corrs =corrcoef(data,'Rows','pairwise');
distances = zeros(length(data));
for i = 1:length(data)
    v_1 = data(i,:);
    for j = 1:length(data)
        v_2= data(j,:);
        if i == j
            continue
        else
            d_ij = abs(v_1 - v_2);
            missVar = isnan(d_ij);
            temp_corrs = corrs(~missVar,~missVar);
            distances(i,j) = sqrt(d_ij(~missVar)*temp_corrs*d_ij(~missVar)');
            distances(j,i) = sqrt(d_ij(~missVar)*temp_corrs*d_ij(~missVar)');
        end
    end
end
filledMatrix = data;
[missRow, ~] = find(isnan(data));
for i = 1:length(missRow)
    [missRow, missCol] = find(isnan(data));
    varMissing = missCol(i);
    currPat = missRow(i);
    missCol(i) = 0;
    missRow(i) = 0;
    tempDistances = distances;
    tempDistances(missRow(find(missCol == varMissing)),:)=NaN;
    tempDistances(:,missRow(find(missCol == varMissing)))=NaN;
    [B,I] = mink(tempDistances(currPat,:),8);
    B(find(I == currPat))=[];
    I(find(I == currPat))=[];
    weights = 1./(B.^2);
    weights = weights./sum(weights);
    imputeValue = dot(weights,data(I,varMissing));
    filledMatrix(currPat,varMissing)=imputeValue;
end
end
%------------- END OF CODE --------------
