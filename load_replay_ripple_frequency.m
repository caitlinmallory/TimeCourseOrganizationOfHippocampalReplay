function ripple_frequency = load_replay_ripple_frequency(full_event_timePoints_exact,cropped_event_timePoints,zscored_LFP)


% Pull out the LFP within this replay event
LFP_sub = compute_dataTemporalConcatenation(zscored_LFP,full_event_timePoints_exact);
ripple_filtered_lfp = LFP_sub(:,3);

if isempty(LFP_sub)
    ripple_frequency = nan;
else

% Find peaks in the ripple-filtered LFP:
[a] = findpeaks(ripple_filtered_lfp);
num_peaks = length(a);

% Ripple frequency = number peaks (cycles)/length of time
ripple_frequency = num_peaks/((full_event_timePoints_exact(2)-full_event_timePoints_exact(1))/30000);
end