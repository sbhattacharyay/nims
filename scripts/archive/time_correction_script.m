times = load('times.mat');
times = times.times;

formatted_times = cellfun(@(x) string(x'),times,'UniformOutput',false);
indexed_times = cell2table(cell(0,2),'VariableNames',{'times','ptIdx'});


for i = 1:length(formatted_times)
   curr_cell = formatted_times{i}; 
   new_cell = [curr_cell,i*ones(length(curr_cell),1)];
   indexed_times = [indexed_times ; array2table(new_cell,'VariableNames',{'times','ptIdx'})];
end

indexed_times.ptIdx = str2double(indexed_times.ptIdx);

writetable(indexed_times,'indexed_times.csv');