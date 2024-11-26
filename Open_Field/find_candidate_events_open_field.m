if isfile("zscored_sd_pyr.mat")
load zscored_sd_pyr.mat
[spike_events, ~] = find_candidate_events_2(zscored_sd(:,3),zscored_sd(:,1));
else
    spike_events = [nan nan nan nan];
end

if isfile('zscored_ripple_power.mat')
load zscored_ripple_power
[ripple_events, ~] = find_candidate_events_2(zscored_ripple_power(:,7),zscored_ripple_power(:,1));
else
    ripple_events = [nan nan nan nan];
end

candidateEvents = struct();
candidateEvents.ripple_events = ripple_events;
candidateEvents.spike_events = spike_events;

save('candidateEvents','candidateEvents');


