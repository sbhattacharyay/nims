%% Authors: Shubhayu Bhattacharyay, B.S. Candidate, Matthew Wang, B.S. Candidate 
% Department of Biomedical Engineering 
% Department of Applied Mathematics and Statistics
% Whiting School of Engineering, Johns Hopkins University
% email address: shubhayu@jhu.edu
% March 2019; Last revision: 23-May-2019
%% ------------- BEGIN CODE --------------
%Plot Accel Time Series
Fs = 10;                                          %# sampling frequency (Hz)

%Set Time Interval to Plot
initTime1 = '00:00:00';
finTime1 = '12:00:00';
initTime2 = '18:00:00';
finTime2 = '23:59:59:999';

sensorLabels = ["LA","LE","LW","RA","RE","RW","Bed"];

for k = 1:length(data_copy)-1
startTime = datenum('00:00','HH:MM');%recording start time
subplot(7,1,k)
hold on;
t = (0:length(data_copy{1,k+1})-1)./ Fs;%time in seconds
t = t/3600/24 + startTime; %time in days (serial date)
plot(t,data_copy{1,k+1},'r')
plot(t,data_copy{2,k+1},'b')
plot(t,data_copy{3,k+1},'k')
ylabel(sensorLabels(k))
ylim([-2 2]);
xlim([datenum(initTime,'HH:MM:SS'),datenum(finTime,'HH:MM:SS')]);
xticks([datenum(initTime,'HH:MM:SS'),datenum(finTime,'HH:MM:SS')]);
xticklabels(string(datestr([datenum(initTime,'HH:MM:SS'),datenum(finTime,'HH:MM:SS')],'HH:MM:SS'))');
set(gca, 'FontName', 'Myriad Pro');
end

subplot(7,1,7)
hold on;
t = (0:length(data_copy{1,1})-1)./ Fs;            %# time in seconds
t = t/3600/24 + startTime;                        %# time in days (serial date)
plot(t,data_copy{1,1},'r')
plot(t,data_copy{2,1},'b')
plot(t,data_copy{3,1},'k')
xlim([datenum(initTime,'HH:MM:SS'),datenum(finTime,'HH:MM:SS')]);
ylabel('Bed')
ylim([-2 2]);
xticks([datenum(initTime,'HH:MM:SS'),datenum(finTime,'HH:MM:SS')]);
xticklabels(string(datestr([datenum(initTime,'HH:MM:SS'),datenum(finTime,'HH:MM:SS')],'HH:MM:SS'))');
set(gca, 'FontName', 'Myriad Pro');
%% Plot Feature Space
Fs = 10;                                          %# sampling frequency (Hz)

%Set Time Interval to Plot
initTime = '07:00:00';
finTime = '17:00:00';
sensorLabels = ["LA","LE","LW","RA","RE","RW","Bed"];

for k = 1:length(SMA)-1
startTime = datenum('00:00','HH:MM');%recording start time
subplot(7,1,k)
hold on;
plot(SMA{2,k+1},SMA{1,k+1},'r')
ylabel(sensorLabels(k))
xlim([datenum(initTime,'HH:MM:SS'),datenum(finTime,'HH:MM:SS')]);
datetick('x','HH:MM:SS');
% xticks([datenum(initTime,'HH:MM:SS'),datenum(finTime,'HH:MM:SS')]);
% xticklabels(string(datestr([datenum(initTime,'HH:MM:SS'),datenum(finTime,'HH:MM:SS')],'HH:MM:SS'))');
axis tight
end               

subplot(7,1,7)
hold on;
t = (0:length(SMA{1,1})-1)./ Fs;%time in seconds
t = t/3600/24 + startTime; %time in days (serial date)
plot(SMA{2,k+1},SMA{1,1},'r');
xlim([datenum(initTime,'HH:MM:SS'),datenum(finTime,'HH:MM:SS')]);
ylabel('Bed')
datetick('x','HH:MM:SS');
axis tight
% xticks([datenum(initTime,'HH:MM:SS'),datenum(finTime,'HH:MM:SS')]);
% xticklabels(string(datestr([datenum(initTime,'HH:MM:SS'),datenum(finTime,'HH:MM:SS')],'HH:MM:SS'))');