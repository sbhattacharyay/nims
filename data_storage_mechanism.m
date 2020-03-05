tic

studyPatientsPY = [2	3	4	5	6	7	8	9	10	11	12	13 ...
    14	15	16	17	18	19	20	21	22	23	24	26	27	28	29	30 ....
    31	32	33	34	35	36	37	38	39	40	41	42	43	49	51	46 .....
    47	48	50	52	53	54	55	56	57	59	60	61	62	63	64	65 ......
    67	68]';

cd D:\
directory = dir;
alreadyLoaded = [];
for i = 1:length(directory)
    if startsWith(directory(i).name,'data')
        alreadyLoaded = [alreadyLoaded; extractBetween(directory(i).name,5,6)];
    end
end
alreadyLoaded = double(string(alreadyLoaded));

patientsLeft = setdiff(studyPatientsPY,alreadyLoaded);

cd Z:\
path = pwd;
directory = dir;

py_folders = zeros(length(directory),1);
for i = 1:length(directory)
    for j = 1:length(studyPatientsPY)
        patientNum = num2str(studyPatientsPY(j),'%02.f');
        if startsWith(directory(i).name,'PY') && endsWith(directory(i).name,patientNum)
            py_folders(i) = i;
        end
    end
end

py_folders = find(py_folders);
folderNames = cell(length(py_folders),1);
for j = 1:length(py_folders)
    folderNames{j} = directory(py_folders(j)).name;
end
patientNums = sort(studyPatientsPY)';
toc

for patIdx = 1:length(patientsLeft)
    tic
    lub = patientsLeft(patIdx);
    idxOfFocus = find(studyPatientsPY == lub);
    folder_of_interest = ['Z:\' folderNames{idxOfFocus}];    
    destinationFolderPathway =['D:\data' folder_of_interest(end-1:end) '.mat'];
    disp(['Patient No.' folder_of_interest(end-1:end) ' initiated.']);
    cd(folder_of_interest);
    try
        copyfile('data.mat',destinationFolderPathway)
    catch
        disp(['Patient No.' folder_of_interest(end-1:end) ' is a problem case.']);
    end
    toc
end