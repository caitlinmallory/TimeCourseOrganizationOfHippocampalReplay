function timepoint_in_epoch = check_whether_timepoint_is_within_bounds(timepoint,time_bounds)

% time bounds is an nX2 matrix containing the start and stop time points
% you want to check. Each row is one epoch.

timepoint_in_epoch = zeros(size(time_bounds,1),1);
for i = 1:size(time_bounds,1)
    if timepoint >= time_bounds(i,1) & timepoint <= time_bounds(i,2)
       timepoint_in_epoch(i) = 1;
    end
end

if any(timepoint_in_epoch == 1)
    timepoint_in_epoch = 1;
else
    timepoint_in_epoch = 0;
end