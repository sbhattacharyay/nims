function data = read_files()
% READ_FILES captures the data from the patient's seven sensors.
% read_files() returns an array of cell arrays with the data where each
% cell array corresponds to one sensor.
len = zeros(1,6);

fid = fopen('dataLA.txt');
C2 = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %s');
fclose(fid);
len(1) = length(C2{1});

fid = fopen('dataLE.txt');
C3 = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %s');
fclose(fid);
len(2) = length(C3{1});

fid = fopen('dataLW.txt');
C4 = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %s');
fclose(fid);
len(3) = length(C4{1});

fid = fopen('dataRA.txt');
C5 = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %s');
fclose(fid);
len(4) = length(C5{1});

fid = fopen('dataRE.txt');
C6 = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %s');
fclose(fid);
len(5) = length(C6{1});

fid = fopen('dataRW.txt');
C7 = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %s');
fclose(fid);
len(6) = length(C7{1});


fid = fopen('databed.txt');
C1 = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %s');
fclose(fid);

data = [C1;C2;C3;C4;C5;C6;C7]';
end