clearvars -except dayFiles day directory rat windows hand_clustered_only

append_candidateEvents = 1;
compute_ripple_power_sd_metrics = 0;
compute_ripple_frequency = 0;
compute_cell_participation_metrics = 0;
compute_spike_timing_metrics = 0;
compute_LFP_surrounding_ripples = 0;
compute_delta_power = 0;
compute_behavioral_metrics = 0;
compute_posterior_metrics = 0;
compute_rank_order_correlation = 1;

durationThr = 0.05; % Only applies to John's filtering method. If the event is less than 50 ms, don't bother analyzing.
params.spikeSampRate = 30000;
% Cell Participation Parameters
params.restrict_to_well_isolated_cells = 0;
params.min_peak_fr_thr = 0;
params.min_num_spikes_to_participate_in_replay = 1;
params.num_fields_thr = 0;
% Behavioral Parameters:
params.speedThr = 5; % cm/s
params.end_zone_size = 30; % cm
params.num_bins_to_mask_track_ends = 15; % number of bins to exclude on either side of the track when calculating total posterior in each map.
params.num_bins_to_mask_local_posterior = 15; % number of bins to exlude around the rat's current location when calculated total posterior in each map.
params.num_time_bins_to_check_for_local_start = 4;
params.distance_thr_for_local_start = 15;  % max number of bins between rat's current location and the start of the replay sequence to be considered local.
params.posBinWidth = 2;
params.davidson_distance_from_line_to_search = 6; % num of bins above OR below the line of best fit to sum posterior over.


if append_candidateEvents == 1
    load decoder_candidateEvents.mat
else
    decoder_events = struct();
end

load Experiment_Information
load Analysis_Information
load binDecoding_02
load clusters
load laser_state
load spikeDensity
load candidateEvents
load Position_Data
load Behavior_Data

if compute_ripple_power_sd_metrics==1
    % Load zscored ripple power and sd amplitude
    load zscored_delta_power
    load zscored_ripple_power
    zscored_LFP_ripple = zscored_ripple_power;
    load zscored_ripple_power_15;
    zscored_LFP_15_ripple = zscored_ripple_power;
    load zscored_sd
    zscored_spikeDensity = zscored_sd;
    load zscored_sd_pyr
    zscored_spikeDensity_pyr = zscored_sd;
    load zscored_sd_pyr_smoothing_15;
    zscored_spikeDensity_pyr_smoothing_15 = zscored_sd;
    load zscored_sd_inhibitory_smoothing_15;
    zscored_spikeDensity_inhibitory_smoothing_15 = zscored_sd;
end
%load LFPs
ripple_tetrode = Experiment_Information.ripple_refs;
lfp_file = ['LFP_Data' num2str(ripple_tetrode) '.mat'];
load(lfp_file);

sessionNum_decoder = 1;
decoding_timeBins = decoder_binDecoding(sessionNum_decoder).timeBins;
decoding_timeBin_centers = mean(decoder_binDecoding(sessionNum_decoder).timeBins,2);
params.decodingWindowShift = decoder_binDecoding(1).shiftSizeDecoding; % time in sec
params.decodingWindowSize = decoder_binDecoding(1).windowSizeDecoding; % time in sec


start_stop_times_day= [min(cat(2, (cat(2, Experiment_Information.Run_Times{:})), Experiment_Information.Sleep_Times{:})) max(cat(2, (cat(2, Experiment_Information.Run_Times{:})), Experiment_Information.Sleep_Times{:}))];
times_day = load_timeBins_cm(start_stop_times_day,params.decodingWindowShift*params.spikeSampRate,params.decodingWindowSize*params.spikeSampRate);

%laser state
laser_state = compute_dataTemporalConcatenation(laser_state,start_stop_times_day);
laser_state = compute_dataInterpolation(laser_state,decoding_timeBin_centers,[]);
laser_state(:,2) = round(laser_state(:,2));


%positions
Position_Data = load_positions_full(Experiment_Information.Run_Times,Experiment_Information.Sleep_Times,decoding_timeBin_centers,Position_Data,4);
Position_Data(isnan(Position_Data(:,5)),5) = 0;
Position_Data_scaled = compute_locsToBins(Position_Data,numSpatialBins,x_edges,y_edges); % Position data in bins instead of cm

%Lap behavior
Pass_Info = load_positions_full(Experiment_Information.Run_Times,Experiment_Information.Sleep_Times,decoding_timeBin_centers,Pass_Info,[]);

% pull out the start times for each segment of the day- this will be used
% to compute the time into the session that each replay occured.
segment_start_times = nan(length(Experiment_Information.Segments),1);
for i = 1:length(Experiment_Information.Segments)
    segment_start_times(i) = [Experiment_Information.Segments(i).Times(1)];
end

for i = [6]

    if i == 1
        candidate_events = candidateEvents.ripple_events;
        candidate_event_parent = 'ripple_events';
    elseif  i == 2
        candidate_events = candidateEvents.spike_events;
        candidate_event_parent = 'spike_events';
    elseif i == 4
        candidate_events = [[candidateEvents.filtering.boundaries_direction_1; candidateEvents.filtering.boundaries_direction_2]  [ones(length(candidateEvents.filtering.boundaries_direction_1),1); 2*ones(length(candidateEvents.filtering.boundaries_direction_2),1)]];
    elseif i == 5
        candidate_events = candidateEvents.filtered_ripple_events;
        candidate_event_parent = 'ripple_events';
    elseif i == 6
        candidate_events = candidateEvents.filtered_spike_events;
        candidate_event_parent = 'spike_events';
    end

    if append_candidateEvents==0
        decoder_events(i).replayEvents = struct();
    end
    for j = 1:length(candidate_events)

        % For most replay analyses, just having the start and stop time of the
        % event should be enough to compute metrics.

        % full_event_boundaries are the indices of decoding_timeBin_centers
        % marking the start and stop of the full candidate event
        % cropped_event_boundaries are the indices of decoding_timeBin_centers
        % marking the start and stop of the best replay segment within the candidate event
        % full_event_timeBins are all the decoding_timeBins for the
        % candidate event
        % cropped_event_timeBins are the decoding_timeBins for the best
        % replay segment within the candidate event (NaN's removed).
        % full_event_timePoints are the start and stop times of the
        % candidate event
        % cropped_even_time_points are the start and stop times of the best
        % replay segment within the candidate event.


        if i == 1 || i == 2
            full_event_timePoints_exact = [candidateEvents.(candidate_event_parent)(j,1) candidateEvents.(candidate_event_parent)(j,2)];
            full_event_boundaries = load_event_boundaries(candidate_events(j,1),candidate_events(j,2),decoding_timeBin_centers);
            cropped_event_boundaries = full_event_boundaries;
            event_NaN = [];
            full_event_inds = (full_event_boundaries(1):full_event_boundaries(2))';
            cropped_event_inds = full_event_inds;
            cropped_event_inds(isnan(event_NaN)) = [];
            full_event_timeBins = decoding_timeBins(full_event_inds,:);
            cropped_event_timeBins = decoding_timeBins(full_event_inds,:);
            cropped_event_timeBins(isnan(event_NaN),:) = [];
        end
        if i == 4 % John's method of replay detection
            %create a list of filtered-event boundaries to consider (column 1 =
            %start times indicies, column 2 = stop indicies, column 3 = in which map this
            %time segment has smoothly moving posterior).

            % for the filtering method, I've already saved:
            %         1) the the start and stop indices of each detected event
            %         2) the decoding bin times
            %         3) the com at each decoding bin time
            %         4) the posterior spread at each decoding bin time
            %         5) the peak posterior at each decoding bin time
            %         6) the indicies of the bins that should be removed due to high posterior spread or high animal velocity


            full_event_timePoints_exact = [candidateEvents.filtering.x_direction_1(candidate_events(j,1),1) candidateEvents.filtering.x_direction_1(candidate_events(j,2),1)];
            full_event_boundaries = candidate_events(j,1:2);
            filtering_map = candidate_events(j,3);
            cropped_event_boundaries = candidate_events(j,:); % in this case full and cropped are the same thing
            if filtering_map ==1 % event_detected in map 1
                event = candidateEvents.filtering.x_direction_1(full_event_boundaries(1):full_event_boundaries(2),2); % x position
                event_NaN = candidateEvents.filtering.x_NaN_direction_1(full_event_boundaries(1):full_event_boundaries(2),2); % x position
            else % event detected in map 2
                event = candidateEvents.filtering.x_direction_2(full_event_boundaries(1):full_event_boundaries(2),2); % x position
                event_NaN = candidateEvents.filtering.x_NaN_direction_2(full_event_boundaries(1):full_event_boundaries(2),2); % x position
            end
            full_event_inds = (full_event_boundaries(1):full_event_boundaries(2))';
            cropped_event_inds = full_event_inds;
            cropped_event_inds(isnan(event_NaN)) = [];
            full_event_timeBins = candidateEvents.filtering.timeBins(full_event_inds,:);
            cropped_event_timeBins = full_event_timeBins;
            cropped_event_timeBins(isnan(event_NaN),:) = [];
        end

        if i == 5 || i == 6
            full_event_timePoints_exact = [candidateEvents.(candidate_event_parent)(j,1) candidateEvents.(candidate_event_parent)(j,2)];
            full_event_boundaries = candidate_events(j).full_event_boundaries;
            cropped_event_boundaries = candidate_events(j).boundaries;
            event_NaN = candidate_events(j).replay_NaN;
            filtering_map = candidate_events(j).bestMap;
            full_event_inds = (full_event_boundaries(1):full_event_boundaries(2))';
            cropped_event_inds = full_event_inds;
            cropped_event_inds(isnan(event_NaN)) = [];
            full_event_timeBins = decoding_timeBins(full_event_inds,:);
            cropped_event_timeBins = decoding_timeBins(full_event_inds,:);
            cropped_event_timeBins(isnan(event_NaN),:) = [];
        end

        full_event_timePoints = [full_event_timeBins(1,1) full_event_timeBins(end,2)];

        in_a_segment_of_interest = zeros(length(Experiment_Information.Segments),1);
        for session = 1:length(Experiment_Information.Segments)
            in_a_segment_of_interest(session) = full_event_timePoints(1) >= Experiment_Information.Segments(session).Times(1) & full_event_timePoints(2) <= Experiment_Information.Segments(session).Times(2);
        end

        if sum(in_a_segment_of_interest)==0
            good_candidate_event = 0;
            cropped_event_timeBins = [];
            cropped_event_timePoints = [];
            decoder_events(i).replayEvents(j).duration = nan;
            decoder_events(i).replayEvents(j).timeBins =  [];
            decoder_events(i).replayEvents(j).timePoints =  [nan nan];
            continue
        end

        %% Add the duration of the underlying candidate event, the best sequence within the candidate event
        decoder_events(i).replayEvents(j).duration_og = (full_event_timePoints(2)-full_event_timePoints(1))/params.spikeSampRate;
        decoder_events(i).replayEvents(j).timeBins_og =  full_event_timeBins;
        decoder_events(i).replayEvents(j).timePoints_og =  full_event_timePoints;

        % If cropped_event_boundaries is empty, there was no good segment
        % within this candidate event. You may still want to consider the
        % animal's location, but it doesn't make sense to measure any other
        % features about the event.

        if isempty(cropped_event_boundaries)
            good_candidate_event = 0;
            cropped_event_timeBins = [];
            cropped_event_timePoints = [];
            decoder_events(i).replayEvents(j).duration = nan;
            decoder_events(i).replayEvents(j).timeBins =  [];
            decoder_events(i).replayEvents(j).timePoints =  [nan nan];
        elseif isnan(cropped_event_boundaries)
            good_candidate_event = 0;
            cropped_event_timeBins = [];
            cropped_event_timePoints = [];
            decoder_events(i).replayEvents(j).duration = nan;
            decoder_events(i).replayEvents(j).timeBins =  [];
            decoder_events(i).replayEvents(j).timePoints =  [nan nan];
        else
            good_candidate_event = 1;
            cropped_event_timePoints = [cropped_event_timeBins(1,1) cropped_event_timeBins(end,2)];
            decoder_events(i).replayEvents(j).duration = (cropped_event_timePoints(2)-cropped_event_timePoints(1))/params.spikeSampRate;
            decoder_events(i).replayEvents(j).timeBins =  cropped_event_timeBins;
            decoder_events(i).replayEvents(j).timePoints =  cropped_event_timePoints;
        end

        if i ==4 && (isnan(decoder_events(i).replayEvents(j).duration) || decoder_events(i).replayEvents(j).duration < durationThr)
            continue
        end

        %% Ripple power and spike density:
        if compute_ripple_power_sd_metrics==1
            ripple_power_sd_metrics = load_replay_ripple_power_spike_density(full_event_timePoints_exact,cropped_event_timePoints,zscored_LFP_ripple,zscored_LFP_15_ripple,zscored_spikeDensity_pyr,zscored_spikeDensity_pyr_smoothing_15,zscored_spikeDensity_inhibitory_smoothing_15);
            decoder_events(i).replayEvents(j).ripple_power_sd_metrics = ripple_power_sd_metrics;
        end
        %% Ripple frequency:
        if compute_ripple_frequency == 1
            ripple_frequency = load_replay_ripple_frequency(full_event_timePoints_exact,cropped_event_timePoints,zscored_LFP_ripple);
            decoder_events(i).replayEvents(j).ripple_frequency = ripple_frequency;
        end
        %% Delta power:
        if compute_delta_power==1
            ripple_window = 0.05; % amount of time in seconds on either side of the ripple to look for delta.
            if i == 4
                delta_power = load_bandpower_surrounding_ripples(zscored_delta_power,[0 0 mean(cropped_event_timePoints)],ripple_window);
            else
                delta_power = load_bandpower_surrounding_ripples(zscored_delta_power,candidateEvents.(candidate_event_parent)(j,:),ripple_window);
            end
            decoder_events(i).replayEvents(j).delta_power = delta_power;
        end
        %% Cell Participation Metrics:
        if compute_cell_participation_metrics==1
            cell_participation_metrics = load_cellParticipationInReplays(params,good_candidate_event,clusters,full_event_timeBins);
            metrics = fieldnames(cell_participation_metrics);
            for m = 1:length(metrics)
                decoder_events(i).replayEvents(j).(metrics{m}) = cell_participation_metrics.(metrics{m});
            end
        end
        %% Cell spike timing during replays/ripples
        if compute_spike_timing_metrics==1
            spike_timing_metrics = load_spikeTimingWithinEvents(params,good_candidate_event,clusters,full_event_timePoints);
            metrics = fieldnames(spike_timing_metrics);
            for m = 1:length(metrics)
                decoder_events(i).replayEvents(j).(metrics{m}) = spike_timing_metrics.(metrics{m});
            end
        end
        %% Load a segment of LFP surrounding each ripple
        if compute_LFP_surrounding_ripples==1
            ripple_window = 1;
            if i == 4
                decoder_events(i).replayEvents(j).LFP = load_LFP_surrounding_ripples(LFP_Data,[0 0 mean(cropped_event_timePoints)],ripple_window);
            else
                decoder_events(i).replayEvents(j).LFP = load_LFP_surrounding_ripples(LFP_Data,candidateEvents.(candidate_event_parent)(j,:),ripple_window);
            end
        end
        %% Behavioral information
        if compute_behavioral_metrics==1
            behaviorMetrics = load_behavioral_measures_during_replay(params,full_event_timePoints,full_event_inds,laser_state,Experiment_Information,Pass_Info,Position_Data,Position_Data_scaled,...
                segment_start_times,Pass_Transitions,Reward_Epoch_Time_Boundaries_drinking,Reward_Epoch_Time_Boundaries_speed_thresholded);
            metrics = fieldnames(behaviorMetrics);
            for m = 1:length(metrics)
                decoder_events(i).replayEvents(j).(metrics{m}) = behaviorMetrics.(metrics{m});
            end
        end
        %% Add properties related to the decoded content of the replay (i.e., properties that can be directly computed from the posterior)
        if compute_posterior_metrics == 1
            posterior_full = decoder_binDecoding.posterior(full_event_inds,:);
            posterior = decoder_binDecoding.posterior(cropped_event_inds,:);

            posterior_spread_full = decoder_binDecoding.posteriorSpread(full_event_inds,:);
            posterior_spread_cropped = decoder_binDecoding.posteriorSpread(cropped_event_inds,:);
            if isempty(posterior)
                decoder_events(i).replayEvents(j).dispersion = nan;
            else

                %%    Adding dispersion to 1d: take this part out eventually:

                best_map = filtering_map;
                num_spatial_bins = size(posterior_full,2)/2;
                left_bins = 1:num_spatial_bins;
                right_bins = num_spatial_bins+1:size(posterior_full,2);

                posterior_left = posterior(:,left_bins);
                posterior_right = posterior(:,right_bins);
                posterior_full_left = posterior_full(:,left_bins);
                posterior_full_right = posterior_full(:,right_bins);
                if best_map == 1
                    best_posterior = posterior_left;
                    best_posterior_full = posterior_full_left;
                else
                    best_posterior = posterior_right;
                    best_posterior_full = posterior_full_right;
                end
                % For each time bin, find the location corresponding to the peak posterior,
                % and the center of mass of the posterior. Also find the posterior
                % associated with each of these locations.
                [posterior_at_peak, replay_sequence_peak] = max(best_posterior,[],2);
                replay_sequence_com = nan(size(best_posterior,1),1);
                posterior_at_com = nan(size(best_posterior,1),1);
                for t = 1:size(best_posterior,1)
                    replay_sequence_com(t) = round(compute_centerOfMass_1D(best_posterior(t,:)));
                    posterior_at_com(t) = best_posterior(t,replay_sequence_com(t));
                end

                dispersion = sqrt( nanmean( (replay_sequence_com(:,1)-nanmean(replay_sequence_com(:,1))).^2 ));
                decoder_events(i).replayEvents(j).dispersion = dispersion;
            end
            %%

            %             if i == 4 || i == 5 || i == 6 % bottom up filtering methods
            %                 eventMetrics = load_candidateEventMetricsDirectionalDecoding(params,posterior_full,posterior,posterior_spread_full,posterior_spread_cropped,cropped_event_timePoints,full_event_inds,filtering_map,[],Position_Data,Position_Data_scaled,Experiment_Information);
            %             else
            %                 eventMetrics = load_candidateEventMetricsDirectionalDecoding(params,posterior_full,posterior,posterior_spread_full,posterior_spread_cropped, full_event_timePoints,full_event_inds,[],[],Position_Data,Position_Data_scaled,Experiment_Information);
            %             end
            %             metrics = fieldnames(eventMetrics);
            %             for m = 1:length(metrics)
            %                 decoder_events(i).replayEvents(j).(metrics{m}) = eventMetrics.(metrics{m});
            %             end
        end

    if compute_rank_order_correlation==1
        [rho_left,p_left,rho_right,p_right] = rank_order_correlation(clusters,full_event_timePoints_exact,x_centers);
        decoder_events(i).replayEvents(j).spearman_rho_left = rho_left;
        decoder_events(i).replayEvents(j).spearman_pval_left = p_left;
        decoder_events(i).replayEvents(j).spearman_rho_right = rho_right;
        decoder_events(i).replayEvents(j).spearman_pval_right = p_right;
    end



        %% Display progress
        disp([num2str(j) '/' num2str(length(candidate_events))])
    end
end



dt=whos('decoder_events');
if dt.bytes < 2e+09
    save('decoder_candidateEvents.mat','decoder_events')
else
    save('decoder_candidateEvents.mat','decoder_events','-v7.3')
end

