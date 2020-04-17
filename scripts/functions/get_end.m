function end_time = get_end(data)
end_t_char = {};
i = 1;
for C = data
    end_t_char{i} = C{12}{end,1};
    i = i + 1;
end
times = string(end_t_char);
sortedTimes = sort(times);
end_time = sortedTimes(1);
end