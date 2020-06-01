%% Author: Shubhayu Bhattacharyay, B.S. Candidate 
% Department of Biomedical Engineering 
% Department of Applied Mathematics and Statistics
% Whiting School of Engineering, Johns Hopkins University
% email address: shubhayu@jhu.edu
% March 2019; Last revision: 11-March-2019
%% ------------- BEGIN CODE --------------
%Correlation Tests
figure
%Plot Accel Time Series
Fs = 10;                                          %# sampling frequency (Hz)

%Set Time Interval to Plot
initTime = '10:30:10';
finTime = '10:30:25';
% initTime2 = '18:00:00';
% finTime2 = '23:59:59:999';

sensorLabels = ["LA","LE","LW","RA","RE","RW","Bed"];

for k = 1:6
startTime = datenum('00:00','HH:MM');%recording start time
subplot(7,1,k)
hold on;
t = (0:length(data_copy{5,k+1})-1)./ Fs;%time in seconds
t = t/3600/24 + startTime; %time in days (serial date)
plot(t,data_copy{5,k+1},'LineWidth',2)
plot(t,data_copy{6,k+1},'LineWidth',2)
plot(t,data_copy{7,k+1},'LineWidth',2)
plot(t,zeros(length(t),1),'k')
ylabel(sensorLabels(k),'FontWeight','bold','FontSize',20)
ylim([-5 5])
xlim([datenum(initTime,'HH:MM:SS'),datenum(finTime,'HH:MM:SS')]);
xticks([datenum(initTime,'HH:MM:SS'),datenum(finTime,'HH:MM:SS')]);
xticklabels(string(datestr([datenum(initTime,'HH:MM:SS'),datenum(finTime,'HH:MM:SS')],'HH:MM:SS'))');
set(gca, 'FontName', 'Myriad Pro');
box off
set(gca,'XColor','none')
end

subplot(7,1,7)
hold on;
t = (0:length(data_copy{1,1})-1)./ Fs;            %# time in seconds
t = t/3600/24 + startTime;                        %# time in days (serial date)
plot(t,data_copy{5,1},'LineWidth',2)
plot(t,data_copy{6,1},'LineWidth',2)
plot(t,data_copy{7,1},'LineWidth',2)
plot(t,zeros(length(t),1),'k')
xlim([datenum(initTime,'HH:MM:SS'),datenum(finTime,'HH:MM:SS')]);
ylabel('Bed','FontWeight','bold','FontSize',20)
ylim([-5 5])
xticks([datenum(initTime,'HH:MM:SS'),datenum(finTime,'HH:MM:SS')]);
xticklabels(string(datestr([datenum(initTime,'HH:MM:SS'),datenum(finTime,'HH:MM:SS')],'HH:MM:SS'))');
set(gca, 'FontName', 'Myriad Pro');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 4, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [4, 9.125])
box off
set(gca,'XColor','none')

%%
figure
k = 1;
    curve = diff(smaFeatures{patOfInterest,k});
    plot(time(2:end),curve,'b','linewidth',1);
    hold on
    plot(time, ones(length(time),1).*threshold,'k--','linewidth',2)
      plot(time, ones(length(time),1).*0,'k','linewidth',2)

    datetick('x','keeplimits');
    ylim([0 .5]);
xlim([datenum('07:10','HH:MM')+1 datenum('07:15','HH:MM')+1]);
ylabel([sensorLabels(k) ' SMA'],'FontWeight','bold','FontSize',20);
xticklabels(string(datestr([datenum('07:10','HH:MM'),datenum('07:15','HH:MM')],'HH:MM:SS'))');

    box off
set(gca,'XColor','none')

    set(gca, 'FontName', 'Myriad Pro');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 4, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [4, 9.125])
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 3], 'PaperUnits', 'Inches', 'PaperSize', [10, 3])
%% Correlation Matrix
runningMatrix = zeros(6,6);

for i = 1:length(studyPatientsPY)
    if i ~= 6
        runningMatrix =runningMatrix+ corrcoef(cell2mat(corr_SMAFeatures(i,:)),'Rows','pairwise');
    end
end

meanCorrs = runningMatrix./61;

%%
figure
sensorLabels = ["LA","LE","LW","RA","RE","RW"];
image(meanCorrs,'CDataMapping','scaled');
h=colorbar;
ylabel(h, 'Correlation')
set(gca,'ytick',[1:6],'yticklabel',sensorLabels,'FontWeight','bold','FontSize',20);
set(gca,'xtick',[1:6],'xticklabel',sensorLabels,'FontWeight','bold','FontSize',20);
set(gca,'XAxisLocation','top');
set(gca, 'FontName', 'Myriad Pro');

