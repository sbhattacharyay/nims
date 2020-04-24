function [missing_percentages,missing_time_series,missingIdxs] = ...
    characterize_missing_data(sensors,t,studyPatientsPY,display_plots,featureType)

dim_of_sensors=size(sensors);

missing_percentages=cellfun(@(x) sum(isnan(x),2)/length(t),sensors,...
    'UniformOutput',false);

stacked_matrix=[];
sensors_missing=zeros(length(studyPatientsPY),length(t));

for i = 1:(dim_of_sensors(1))
    curr_matrix = sensors{i,6};
    stacked_matrix = [stacked_matrix; curr_matrix];
    sensors_missing = sensors_missing+isnan(curr_matrix);
end

missingID=[0 NaN Inf];

missing_time_series=sum(isnan(stacked_matrix))/(length(studyPatientsPY)*...
    (dim_of_sensors(1)-1));
missingIdxs=cellfun(@(x) ismissing(x,missingID),sensors,'UniformOutput',false);

if display_plots == true
    figure
    x = [t(1) t(end)];
    y = [1 length(studyPatientsPY)];
    colormap(flipud(pink))
    
    image(x,y,sensors_missing,'CDataMapping','scaled')
    ax = gca;
    ax.FontSize = 14;
    if featureType==1
        datetick('x','HH PM');
        xlabel('Time of Day','FontWeight','bold','FontSize',20);
    else
        x_tick_pts=datenum(hours(0):hours(1):hours(8));
        xticks(x_tick_pts)
        datetick('x','HH:MM','keepticks');
        xlabel('Time from Start of Recording','FontWeight','bold','FontSize',20);
    end
    set(gca, 'FontName', 'Myriad Pro');
    h=colorbar;
    ylabel(h, 'No. of Sensors Missing','FontSize',14)
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 10], 'PaperUnits', 'Inches', 'PaperSize', [10, 10])
    ylabel('Patient No.','FontWeight','bold','FontSize',20);
end
end

