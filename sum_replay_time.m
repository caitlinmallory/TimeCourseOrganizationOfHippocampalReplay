function duration = sum_replay_time(replay_timePoints)

duration = 0;
for n = 1:size(replay_timePoints,1)
    duration = duration + (replay_timePoints(n,2)-replay_timePoints(n,1))./30000;
end

if duration == 0
    duration = nan;
end