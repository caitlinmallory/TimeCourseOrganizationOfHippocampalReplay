function event_boundaries = load_event_boundaries(event_start_time,event_stop_time,decoding_time_bin_centers)

% Find the indicies in decoding_time_bin_centers corresponding to the start
% and stop times of the candidate event
event_boundaries(1) = find(hist(event_start_time,decoding_time_bin_centers) == 1);
event_boundaries(2) = find(hist(event_stop_time,decoding_time_bin_centers) == 1);
