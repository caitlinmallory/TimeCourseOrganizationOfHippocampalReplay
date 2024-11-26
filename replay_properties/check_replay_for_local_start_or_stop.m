function [local_replay_start, start_distance_from_rat, local_replay_stop, end_distance_from_rat] = check_replay_for_local_start_or_stop(replay_position,ratPos,spatialDim,num_decoding_bins_to_search,distanceThr)

% replay_position = center of mass of the replay event over time. For linear track
%, use the map (left or right) with greater posterior.
% ratPos = rat's actual x,y position, in bins (only x matters for linear
% track)
% time_window_to_search_for_local_start = time window (s) over which to
% examine the start of the replay (1 bin seems too restrictive).
% distanceThr = the distance (in bins) from the rat the rat can start and
% still be considered a local event.
local_replay_start = 0;
local_replay_stop = 0;

% Check if the replay starts near the animal's current location

%replay_location_in_start_window = mean(replay_position(1:num_decoding_bins_to_average,:),1);
distance_from_rat_replay_start = nan(min(size(replay_position,1),num_decoding_bins_to_search),1);

for i = 1:min(size(replay_position,1),num_decoding_bins_to_search)
    
    if spatialDim == 1
        distance_from_rat_replay_start(i) = abs(ratPos(1)-replay_position(i));
    else
        distance_from_rat_replay_start(i) = sqrt((ratPos(1)-replay_position(i,1))^2 + (ratPos(2)-replay_position(i,2))^2);
    end
    
end

if any(distance_from_rat_replay_start < distanceThr)
    local_replay_start = 1;
end

% Take the average distance over the first few bins
start_distance_from_rat = nanmean(distance_from_rat_replay_start);

% Check if the replay ends near the animal's current location

%replay_location_in_start_window = mean(replay_position(1:num_decoding_bins_to_average,:),1);
distance_from_rat_replay_stop = nan(min(size(replay_position,1),num_decoding_bins_to_search),1);

for i = 1:min(size(replay_position,1),num_decoding_bins_to_search)
    
    if spatialDim == 1
        distance_from_rat_replay_stop(i) = abs(ratPos(1)-replay_position(end-i+1));
    else
        distance_from_rat_replay_stop(i) = sqrt((ratPos(1)-replay_position(end-i+1,1))^2 + (ratPos(2)-replay_position(end-i+1,2))^2);
    end
    
end

if any(distance_from_rat_replay_stop < distanceThr)
    local_replay_stop = 1;
end

% Take the average distance over the first few bins
end_distance_from_rat = nanmean(distance_from_rat_replay_stop);