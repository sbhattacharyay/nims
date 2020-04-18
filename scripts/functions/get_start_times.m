function sortedTimes = get_start_times(data)
% Isolate the start time character arrays.
start_t_char = {};
i = 1;
for C = data
    start_t_char{i} = C{12}{1,1};
    i = i + 1;
end
times = string(start_t_char);
sortedTimes = sort(times);
end