% Convert Patient Clinical Data spreadsheet (XLSX) to CSV

patientData = readtable('patient_clinical_data.xlsx',...
    'TreatAsEmpty',{'.','NA'});
excluded = isnan(patientData.AccelPatientNo_) | patientData.StudyPatientNo_==1;
patientData(excluded,:) = [];
patientData=sortrows(patientData,'AccelPatientNo_');

patientData.AccelRecordingStartTime = datetime( ...
    patientData.AccelRecordingStartTime,'convertfrom','excel');
patientData.AccelRecordingEndTime = datetime( ...
    patientData.AccelRecordingEndTime,'convertfrom','excel');
patientData.AccelDataMountTime = datetime( ...
    patientData.AccelDataMountTime,'convertfrom','excel');

patientData.AccelRecordingStartTime.Format = 'HH:mm';
patientData.AccelRecordingEndTime.Format = 'HH:mm';
patientData.AccelDataMountTime.Format = 'HH:mm';

writetable(patientData,'patient_clinical_data.csv');