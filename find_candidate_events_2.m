function [events, event_amplitude_zscore] = find_candidate_events_2(amplitude_zscore,timestamps)
% In this version, if there are two peaks that either occur without a dip
% below the lo_std_cutoff, or if there are two peaks that come within
% event_min_peak_separation, the time of the ripple will be considered the
% second peak (instead of average).
spikeSampRate = 30000;

load Analysis_Information

events_lo_std_cutoff = 0;
events_hi_std_cutoff = 3;
events_min_peak_separation = 0.07;
events_min_length = 0.0;
events_max_length = inf;

[event_amplitude_zscore,peak_idx] = findpeaks(amplitude_zscore,'minpeakheight',events_hi_std_cutoff);
% Find  the start and stop timepoints for each event
num_events = length(peak_idx);
events=zeros(num_events,3);
events(:,3)= timestamps(peak_idx);
events(:,4) = event_amplitude_zscore;

for n = 1:num_events
    start_idx = peak_idx(n);
    while amplitude_zscore(start_idx) >  events_lo_std_cutoff && start_idx > 1
        start_idx = start_idx - 1;
    end
    events(n,1) = timestamps(start_idx);
    
    stop_idx = peak_idx(n);
    while amplitude_zscore(stop_idx) >  events_lo_std_cutoff && stop_idx < length(amplitude_zscore)
        stop_idx = stop_idx + 1;
    end
    events(n,2) = timestamps(stop_idx);
end


% If 2 peaks occur without a dip below zero, this method will count them as
% 2 events instead of one. Delete duplicates. Use whichever peak was
% highest as the power, and peak location.

for n = 2:num_events
    if events(n,1) == events(n-1,1)
        events(n-1,1) = 0;
        events(n,3) = events(n,3);
        [max_val,max_ind] = max([events(n,4),events(n-1,4)]);
        peak_times = [events(n,3),events(n-1,3)];
        events(n,3) = peak_times(max_ind);
        events(n,4) = max_val;
    end
end

inds_to_remove = find(events(:,1) == 0);
events(inds_to_remove,:) = [];

num_events = size(events,1);
% Merge events that occur too close to one another
for n = 2:num_events
    if events(n,3) - events(n-1,3) <= events_min_peak_separation*spikeSampRate
 
        % Find whichever peak was larger, and save that as the power
        [max_val,max_ind] = max([events(n,4),events(n-1,4)]);
        peak_times = [events(n,3),events(n-1,3)];
        
        events(n,3) = peak_times(max_ind);
        events(n,4) = max_val;
        events(n,1) = events(n-1,1);
        events(n-1,1) = 0;
    end
end

% In merging, the start time of the merged event was set to zero: delete
% these events:
idx_to_remove = find(events(:,1) == 0);
events(idx_to_remove,:) = [];
% event_amplitude_zscore(idx_to_remove) = [];

% Delete events that were too short or too long.
events_to_keep = find(events(:,2)-events(:,1)...
    >= events_min_length*spikeSampRate &...
    (events(:,2)-events(:,1))...
    <= events_max_length*spikeSampRate);

% event_amplitude_zscore = event_amplitude_zscore(events_to_keep);
events = events(events_to_keep,:);

event_amplitude_zscore = events(:,4);
% outputs:

%return the event_amplitude_zscore_mean (if two events were merged, this
%compute the mean amplitude of the 2 events);

%Caitlin added the following on 11/26/20:
%recalculate the power and location of peaks after merging.
% for n = 1:size(events,1)
%     event_start_idx = find(timestamps == events(n,1));
%     event_stop_idx = find(timestamps == events(n,2));
%     event_timestamps = timestamps(event_start_idx:event_stop_idx);
%     [event_amplitude_averaged_across_peaks,peak_idx] = findpeaks(amplitude_zscore(event_start_idx:event_stop_idx),'minpeakheight',events_hi_std_cutoff);
%     event_amplitude_zscore_mean(n) = mean(event_amplitude_averaged_across_peaks);
%     events(n,3) = event_timestamps(round(mean(peak_idx)));
% end

% figure()
% plot(timestamps./spikeSampRate,amplitude_zscore)
% hold on
% plot(events(:,1)./spikeSampRate,zeros(length(events),1),'.g')
% hold on
% plot(events(:,2)./spikeSampRate,zeros(length(events),1),'.r')
% hold on
% plot(events(:,3)./spikeSampRate,zeros(length(events),1),'.m')
