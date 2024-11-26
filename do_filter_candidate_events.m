
clearvars -except dayFiles day directory rat windows hand_clustered_only
load Experiment_Information
load Analysis_Information
load clusters
load Position_Data
load Behavior_Data
load laser_state
load spikeDensity
load candidateEvents

Run_Times = Experiment_Information.Run_Times;
Sleep_Times = Experiment_Information.Sleep_Times;
spikeSampRate = Experiment_Information.spikeSampRate;
spatialDim = Experiment_Information.spatialDim;

if spatialDim == 1
    directionalDecoding = 1;
    numSpatialBins = [1 numSpatialBins(2)];
    load binDecoding_02
else
    load binDecoding_08
end

times_list = cat(2, (cat(2, Run_Times{:})), Sleep_Times{:});
Times_day = [min(times_list) max(times_list)];

% pull out the start times for each segment of the day- this will be used
% to compute the time into the session that each replay occured.
segment_start_times = [];
for i = 1:length(Experiment_Information.Segments)
    segment_start_times = [segment_start_times; Experiment_Information.Segments(i).Times(1)];
end

for sessionNum_decoder = 1:size(Run_Times,1)

    timeBins = decoder_binDecoding(sessionNum_decoder).timeBins;
    decoding_time_bin_centers = mean(decoder_binDecoding(sessionNum_decoder).timeBins,2);

    %extract events:
    %For spike density events, ripple events, and laser events, this is
    %converting from start and and stop times to start and stop indices
    %(with resect to a vector of decoding bin times).
    event_boundaries = {};
    event_boundaries{1}(:,1) = find(hist(candidateEvents.ripple_events(:,1),decoding_time_bin_centers) == 1);
    event_boundaries{1}(:,2) = find(hist(candidateEvents.ripple_events(:,2),decoding_time_bin_centers) == 1);
    event_boundaries{2}(:,1) = find(hist(candidateEvents.spike_events(:,1),decoding_time_bin_centers) == 1);
    event_boundaries{2}(:,2) = find(hist(candidateEvents.spike_events(:,2),decoding_time_bin_centers) == 1);

    for i = [1 2]

        event_boundaries_sub = event_boundaries{i};
        if i == 1
            candidateEvents.filtered_ripple_events = [];
        end
        if i == 2
            candidateEvents.filtered_spike_events = [];
        end


        for j = 1:size(event_boundaries_sub,1)


            full_event_boundaries = [event_boundaries_sub(j,1) event_boundaries_sub(j,2)];
            indSub = event_boundaries_sub(j,1):event_boundaries_sub(j,2);

            ce = filter_candidate_events(Experiment_Information,decoder_binDecoding,sessionNum_decoder,indSub);

            if i == 2

                candidateEvents.filtered_spike_events(j).replay = ce.x(:,2);
                candidateEvents.filtered_spike_events(j).replay_NaN = ce.x_NaN(:,2);
                candidateEvents.filtered_spike_events(j).best_map = ce.best_map;
                candidateEvents.filtered_spike_events(j).full_event_boundaries = full_event_boundaries;
                candidateEvents.filtered_spike_events(j).boundaries = ce.boundaries;
                candidateEvents.filtered_spike_events(j).bestMap = ce.best_map;
                candidateEvents.filtered_spike_events(j).posterior_diff = ce.posterior_diff;
            
            elseif i == 1
 
                candidateEvents.filtered_ripple_events(j).replay = ce.x(:,2);
                candidateEvents.filtered_ripple_events(j).replay_NaN = ce.x_NaN(:,2);
                candidateEvents.filtered_ripple_events(j).best_map = ce.best_map;
                candidateEvents.filtered_ripple_events(j).full_event_boundaries = full_event_boundaries;
                candidateEvents.filtered_ripple_events(j).boundaries = ce.boundaries;
                candidateEvents.filtered_ripple_events(j).bestMap = ce.best_map;
                candidateEvents.filtered_ripple_events(j).posterior_diff = ce.posterior_diff;
            end
        end



    end

end


if exist('candidateEvents.mat') == 2
    save('candidateEvents.mat','candidateEvents','-append')
else
    save('candidateEvents.mat','candidateEvents')
end
