%% Master Script 1: Clinical Data Extraction
% Decoding Quantitative Motor Features for Classification and Prediction
% in Severe Acquired Brain Injury
%
% Shubhayu Bhattacharyay
% Department of Biomedical Engineering
% Department of Applied Mathematics and Statistics
% Whiting School of Engineering, Johns Hopkins University
% email address: shubhayu@jhu.edu
%% ------------- BEGIN CODE --------------
% Load data from most recent patient clinical data spreadsheet
addpath('functions/')

patientData = readtable('../clinical_data/SB_patient_table.xlsx',...
    'TreatAsEmpty',{'.','NA'});
excluded = isnan(patientData.AccelPatientNo_) | patientData.StudyPatientNo_==1;
patientData(excluded,:) = [];
patientData=sortrows(patientData,'AccelPatientNo_');

% Indices of patients to use in study. To see exclusion criteria, please
% refer to Materials and Methods section

GOSE_fav_12mo=NaN(height(patientData),1);
GOSE_fav_12mo(patientData.GOS_E12Months>=5)=1;
GOSE_fav_12mo(patientData.GOS_E12Months<5)=0;

mRS_fav_12mo=NaN(height(patientData),1);
mRS_fav_12mo(patientData.mRS12Months<=2)=1;
mRS_fav_12mo(patientData.mRS12Months>2)=0;

patientData=addvars(patientData,GOSE_fav_12mo);
patientData=addvars(patientData,mRS_fav_12mo);

writetable(patientData,'../clinical_data/patient_clinical_data.csv');

%% Clinical Data Preprocessing and Imputation

% Apply LLF-optimized box-cox transform to quantitative clinical variables
% and normalize dataset to enable hypothesis testing:

% Age:
[Age,lam_age]=boxcox(patientData.NCCUAdmissionAge_y_);
[Age,mu_age,sig_age]=zscore(Age);

% GCS at enrollment:
[GCST,lam_gcs_en]=boxcox(patientData.GCST);
[GCST,mu_gcs_en,sig_gcs_en]=zscore(GCST);

% APACHE at enrollment:
[APACHE,lam_apache]=boxcox(patientData.APACHEII);
[APACHE,mu_apache,sig_apache]=zscore(APACHE);

% APACHE mortality risk:
[APACHE_risk,lam_apache_risk]=boxcox(patientData.APACHEIIMortalityRisk);
[APACHE_risk,mu_apache_risk,sig_apache_risk]=zscore(APACHE_risk);

% % GCS at discharge:
% [gcs_scores_dis,lam_gcs_dis]=boxcox((patientData.gcs_dis));
% [gcs_scores_dis,mu_gcs_dis,sig_gcs_dis]=zscore(gcs_scores_dis);

% Gender (Males are 1, Females are -1):
females = categorical(patientData.Sex)~= 'M';
Sex = (-1).^females;

% Diagnosis codes
CVA = -((-1).^(patientData.CVA));
ICH = -((-1).^(patientData.ICH));
SAH = -((-1).^(patientData.SAH));
BT = -((-1).^(patientData.BrainTumorOrLesion));
SDH = -((-1).^(patientData.SDHOrEDH));
TBI = -((-1).^(patientData.TBI));

PY_pNum=patientData.AccelPatientNo_;
pNum=patientData.StudyPatientNo_;
Death=patientData.DiedDuringThisHospitalStay_;
Death_1yr=patientData.Death12Months;
fav_GOSE_1yr=patientData.GOSE_fav_12mo;
fav_mRS_1yr = patientData.mRS_fav_12mo;

tf_patient_covariates = table(PY_pNum,pNum,Age,Sex,CVA,ICH,SAH,BT,...
    SDH,TBI,APACHE,GCST,APACHE_risk,Death,Death_1yr,fav_GOSE_1yr,....
    fav_mRS_1yr);

writetable(tf_patient_covariates,'../clinical_data/tf_patient_covariates.csv');

%% Visualize distribution of transformed (Box-Cox) clinical variables
figure

subplot(2,2,1);
histogram(Age,30);
title('Age')

subplot(2,2,2);
histogram(APACHE,30);
title('APACHE')

subplot(2,2,3);
histogram(GCST,30);
title('GCS T')

subplot(2,2,4);
histogram(APACHE_risk,30);
title('APACHE Risk')

license('inuse')
%------------- END OF CODE --------------

