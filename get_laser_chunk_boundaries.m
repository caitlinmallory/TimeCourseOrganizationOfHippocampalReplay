function [laser_on_windows, laser_off_windows, long_laser_off_windows] = get_laser_chunk_boundaries(laser_timestamps,spikeSampRate, duty_cycle)

  tolerance = 0.1;


if isempty(duty_cycle)
    duty_cycle = inf;
end

if isempty(laser_timestamps)
    
    long_laser_off_windows = [];
    laser_off_windows = [];
    laser_on_windows = [];
    
else
    if laser_timestamps(1,2) == 0
        laser_timestamps(1,:) = [];
    end
    
    laser_on_inds = find(laser_timestamps(:,2) == 1);
    laser_off_inds = find(laser_timestamps(:,2) == 0);
    
    laser_on_windows = nan(length(laser_on_inds),2);
    laser_off_windows = nan(length(laser_off_inds),2);
    
    for i = 1:length(laser_on_inds)
        
        laser_on_windows(i,1) = laser_timestamps(laser_on_inds(i),1);
        laser_on_windows(i,2) = laser_timestamps(laser_off_inds(find(laser_off_inds > laser_on_inds(i),1,'first'),1)) - 1;
        
    end
    
    for i = 1:length(laser_off_inds)
        laser_off_windows(i,1) = laser_timestamps(laser_off_inds(i),1);
        if ~isempty(find(laser_on_inds > laser_off_inds(i),1,'first'))
            laser_off_windows(i,2) = laser_timestamps(laser_on_inds(find(laser_on_inds > laser_off_inds(i),1,'first'),1)) - 1;
        else
            laser_off_windows(i,:) = [];
        end
    end
    
     laser_on_inds = find((laser_on_windows(:,2)-laser_on_windows(:,1) ) > spikeSampRate*duty_cycle - spikeSampRate*tolerance & ...
        (laser_on_windows(:,2)-laser_on_windows(:,1) < spikeSampRate*duty_cycle + spikeSampRate*tolerance));

    laser_off_inds = find((laser_off_windows(:,2)-laser_off_windows(:,1) ) > spikeSampRate*duty_cycle - spikeSampRate*tolerance & ...
        (laser_off_windows(:,2)-laser_off_windows(:,1) < spikeSampRate*duty_cycle + spikeSampRate*tolerance));

    laser_on_windows = laser_on_windows(laser_on_inds,:);
    laser_off_windows = laser_off_windows(laser_off_inds,:);
    
    long_laser_off_windows = [];
%     long_laser_off_inds = find(laser_off_windows(:,2)-laser_off_windows(:,1) > spikeSampRate*duty_cycle );
%     
%     long_laser_off_windows = laser_off_windows(long_laser_off_inds,:);
%     laser_off_windows(long_laser_off_inds,:) = [];
    
    
end