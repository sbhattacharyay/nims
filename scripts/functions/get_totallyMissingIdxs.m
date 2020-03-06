function [totallyMissingIdxs,output_sensors]=get_totallyMissingIdxs ...
    (sensor_count,feature_count,sensors)

%contains indices of the totally missing time series as follows:
%Row_1 is the sensor number (1 through 7)
%Row_2 is the feature_count (1 through 7) that corresponds to "feature_names"
%Row_3 is the patient number that corresponds to the patient index (1
%through 62)

totallyMissingIdxs = [];
missingID=[0 NaN Inf];
output_sensors = sensors;

for i=1:sensor_count
    for j=1:feature_count
        currMatrix = sensors{i,j};
        [ptRows,~] = size(currMatrix);
        for k = 1:ptRows
            if all(ismissing((currMatrix(k,:)),missingID))
                totallyMissingIdxs = [totallyMissingIdxs [i;j;k]];
                currMatrix(k,:) = NaN(1,length(currMatrix(k,:)));
            end
            output_sensors{i,j} = currMatrix;
        end
    end
end
end
