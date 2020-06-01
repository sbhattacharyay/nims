function plotMissingCurve(missing_time_series,noMotion_time_series,...
    t,featureType,n)

totRec = n*7;

rootDir = generateRootDir('Missing_Data');

figure

plot(t,missing_time_series*100,'r','LineWidth',3);
hold on
area(t,missing_time_series*100,'FaceAlpha',0.5,'EdgeAlpha',0,'FaceColor',[.827,.827,.827])
plot(t,100*(noMotion_time_series+missing_time_series),'b','LineWidth',3);

axis tight
ylim([0,100])
ytickformat('percentage')

ax = gca;
ax.FontSize = 14;
if featureType==1
    datetick('x','HH PM');
    xlabel('Time of Day','FontWeight','bold','FontSize',20);
else
    x_tick_pts=datenum(hours(0):hours(1):hours(8));
    xticks(x_tick_pts)
    datetick('x','HH:MM','keepticks');
    xlabel('Time from Start of Recording','FontWeight','bold',...
        'FontSize',20);
end
set(gca, 'FontName', 'Myriad Pro');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 10], 'PaperUnits', 'Inches', 'PaperSize', [10, 10])
ylabel(['Percentage of Total Recordings (n = ',num2str(totRec),')'],....
    'FontWeight','bold','FontSize',20);

if featureType == 1
    saveas(gcf,[rootDir filesep 'tod_missing_curve.png'])
else
    saveas(gcf,[rootDir filesep 'tfr_missing_curve.png'])
end
end

