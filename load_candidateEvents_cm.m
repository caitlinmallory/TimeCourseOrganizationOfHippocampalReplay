
clearvars -except dayFiles day directory rat windows hand_clustered_only
load Experiment_Information
load Analysis_Information
load clusters
load Position_Data
load Behavior_Data
load laser_state
load spikeDensity
load candidateEvents


restrict_to_high_quality_cells = 1;

kilosort_contamPct_Thr = 20;
mClust_clusterIsolationThr = 0.2; %0.2
mClust_isolation_distanceThr = 40; %40
msort_clusterNoiseOverlapThr  = 0.03; %0.03
msort_clusterIsolationThr = 0.95; %0.95

if restrict_to_high_quality_cells == 1
    % determine which clusters meet criterion for 'good isolation'
    if isfield(clusters,'contamPct')
        well_isolated_cluster = ([clusters.contamPct] < kilosort_contamPct_Thr )';
    elseif isfield(clusters,'L_Ratio')
        well_isolated_cluster = ([clusters.L_Ratio] < mClust_clusterIsolationThr)';
    elseif isfield(clusters,'isolation')
        well_isolated_cluster = ([clusters.isolation] >= msort_clusterIsolationThr & [clusters.noise_overlap] <=msort_clusterNoiseOverlapThr)';
    else
        keyboard
    end
else
    well_isolated_cluster = ones(length(clusters),1);
end

good_cluster_inds = find(well_isolated_cluster==1 & [clusters.Excitatory]'==1);
bimodal_inds = find(well_isolated_cluster == 1 & [clusters.Modality]' == 2);
unimodal_inds = find(well_isolated_cluster == 1 & [clusters.Modality]' == 1);
excitatory_inds = find([clusters.Excitatory]' == 1);
inhibitory_inds = find([clusters.Excitatory]' == 0);

num_bins_to_check_for_local_start = 4;
bin_thr_for_local_start = 16;
num_position_bins_to_mask_out_local_posterior = 16;


if ~isnan(Experiment_Information.reward_zone_size)
    num_bins_in_reward_zone = round(Experiment_Information.reward_zone_size/binSize);
else
    num_bins_in_reward_zone = 16;
end

speedRangeThr = [0 speedThr];

Run_Times = Experiment_Information.Run_Times;
Sleep_Times = Experiment_Information.Sleep_Times;
spikeSampRate = Experiment_Information.spikeSampRate;
spatialDim = Experiment_Information.spatialDim;
Track_Type = Experiment_Information.Track_Type;
if spatialDim == 1
    directionalDecoding = 1;
    numSpatialBins = [1 numSpatialBins(2)];
    load binDecoding_02
else
    load binDecoding_08
end

times_list = cat(2, (cat(2, Run_Times{:})), Sleep_Times{:});
Times_day = [min(times_list) max(times_list)];

%positions
times = load_timeBins_cm(Times_day,decoder_binDecoding(1).shiftSizeDecoding*spikeSampRate,decoder_binDecoding(1).windowSizeDecoding*spikeSampRate);
Position_Data_full = load_positions_full(Run_Times,Sleep_Times,times,Position_Data,4);
Position_Data_full(isnan(Position_Data_full(:,5)),5) = 0;

%Lap behavior
Pass_Info_full = load_positions_full(Run_Times,Sleep_Times,times,Pass_Info,2);


% this is a cell that will hold a structure for each candidate event type
% 1 = Ripple_Events, 2 = Spike_Density_Events, 3 = Laser_Triggered_Events,
% 4 = John Widloski's filtering method


% pull out the start times for each segment of the day- this will be used
% to compute the time into the session that each replay occured.
segment_start_times = [];
for i = 1:length(Experiment_Information.Segments)
    segment_start_times = [segment_start_times; Experiment_Information.Segments(i).Times(1)];
end

% Load zscored ripple power and sd amplitude
load zscored_ripple_power
zscored_LFP = zscored_ripple_power;
load zscored_sd
zscored_spikeDensity = zscored_sd;
load zscored_sd_pyr
zscored_spikeDensity_pyr = zscored_sd;


reward_zone_entries = Pass_Transitions(Pass_Transitions(:,2) >=1,:);
reward_zone_exits = Pass_Transitions(Pass_Transitions(:,2) < 0,:);

for sessionNum_decoder = 1:size(Run_Times,1)

    timeBins = decoder_binDecoding(sessionNum_decoder).timeBins;
    decoding_time_bin_centers = mean(decoder_binDecoding(sessionNum_decoder).timeBins,2);

    laser_state_sub = compute_dataTemporalConcatenation(laser_state,Times_day);
    laser_state_sub = compute_dataInterpolation(laser_state_sub,decoding_time_bin_centers,[]);
    laser_state_sub(:,2) = round(laser_state_sub(:,2));

    Position_Data_sub = compute_dataTemporalConcatenation(Position_Data,Times_day);
    Position_Data_sub = compute_dataInterpolation(Position_Data_sub,decoding_time_bin_centers,[]);
    Position_Data_sub_scaled = compute_locsToBins(Position_Data_sub,numSpatialBins,x_edges,y_edges);

    %extract events:
    %For spike density events, ripple events, and laser events, this is
    %converting from start and and stop times to start and stop indices
    %(with resect to a vector of decoding bin times).
    event_boundaries = {};
    event_boundaries{1}(:,1) = find(hist(candidateEvents.ripple_events(:,1),decoding_time_bin_centers) == 1);
    event_boundaries{1}(:,2) = find(hist(candidateEvents.ripple_events(:,2),decoding_time_bin_centers) == 1);
    event_boundaries{2}(:,1) = find(hist(candidateEvents.spike_events(:,1),decoding_time_bin_centers) == 1);
    event_boundaries{2}(:,2) = find(hist(candidateEvents.spike_events(:,2),decoding_time_bin_centers) == 1);

    if ~isempty(candidateEvents.laser_events)
        event_boundaries{3}(:,1) = find(hist(candidateEvents.laser_events(:,1),decoding_time_bin_centers) == 1);
        event_boundaries{3}(:,2) = find(hist(candidateEvents.laser_events(:,2),decoding_time_bin_centers) == 1);
    end
    if ~isempty(candidateEvents.filtered_ripple_events)
        boundaries = load_fieldVec(candidateEvents.filtered_ripple_events,'full_event_boundaries',2);

        event_boundaries{5}(:,1) = boundaries(:,1);
        event_boundaries{5}(:,2) = boundaries(:,2);
        cropped_event_boundaries{5} = nan(size(candidateEvents.filtered_ripple_events,2),2);
        for row = 1:size(candidateEvents.filtered_ripple_events,2)
            if ~isnan(candidateEvents.filtered_ripple_events(row).boundaries)
                cropped_event_boundaries{5}(row,:) = candidateEvents.filtered_ripple_events(row).boundaries;
            end
        end
    end
    if ~isempty(candidateEvents.filtered_spike_events)
        boundaries = load_fieldVec(candidateEvents.filtered_spike_events,'full_event_boundaries',2);

        event_boundaries{6}(:,1) = boundaries(:,1);
        event_boundaries{6}(:,2) = boundaries(:,2);
        cropped_event_boundaries{6} = nan(size(candidateEvents.filtered_spike_events,2),2);
        for row = 1:size(candidateEvents.filtered_spike_events,2)
            if ~isnan(candidateEvents.filtered_spike_events(row).boundaries)
                cropped_event_boundaries{6}(row,:) = candidateEvents.filtered_spike_events(row).boundaries;
            end
        end
    end


    if isfield(candidateEvents,'filtering') == 1
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

        filtered_map_1_boundaries = candidateEvents.filtering(sessionNum_decoder).boundaries_direction_1;
        filtered_map_2_boundaries = candidateEvents.filtering(sessionNum_decoder).boundaries_direction_2;

        event_boundaries{4} = [filtered_map_1_boundaries; filtered_map_2_boundaries];
        event_boundaries{4}(:,3) = [ones(length(filtered_map_1_boundaries),1); 2*ones(length(filtered_map_2_boundaries),1)];
        %         event_boundaries{4} = sortrows(event_boundaries{4},1);
    end

    if spatialDim == 2
        %compute properties of candidate events and load struct
        for i = 1:length(event_boundaries)
            event_boundaries_sub = event_boundaries{i};
            for j = 1:size(event_boundaries_sub,1)
                indSub = event_boundaries_sub(j,1):event_boundaries_sub(j,2);
                indNaN = find(isnan(x_NaN(indSub,2)));
                replay = x(indSub,2:3);
                replay_NaN = x_NaN(indSub,2:3);
                replay_NaNremoved = replay_NaN; replay_NaNremoved(indNaN,:) = [];
                times_NaNremoved = x(indSub,1); times_NaNremoved(indNaN) = [];
                decoder_events{i}(sessionNum_decoder).replayEvents(j).indNaN = indNaN;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).indData = indSub(setdiff(1:length(indSub),indNaN))';
                decoder_events{i}(sessionNum_decoder).replayEvents(j).timeBins = timeBins(indSub,:);
                decoder_events{i}(sessionNum_decoder).replayEvents(j).timePoints = [decoding_time_bin_centers(indSub(1),1) decoding_time_bin_centers(indSub(end),1)];
                decoder_events{i}(sessionNum_decoder).replayEvents(j).duration = (decoding_time_bin_centers(indSub(end),1)-decoding_time_bin_centers(indSub(1),1))/spikeSampRate;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).maxJump_NaN = max(compute_sequenceJumps(replay_NaN));
                decoder_events{i}(sessionNum_decoder).replayEvents(j).maxJump_NaNremoved = max(compute_sequenceJumps(replay_NaNremoved));
                decoder_events{i}(sessionNum_decoder).replayEvents(j).maxJump_NaNremoved_time = max(diff(times_NaNremoved));
                decoder_events{i}(sessionNum_decoder).replayEvents(j).replay = replay;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).dispersion = compute_sequenceDispersion(replay_NaN);
                decoder_events{i}(sessionNum_decoder).replayEvents(j).ratLoc = [Position_Data_sub_scaled(indSub(1),2),Position_Data_sub_scaled(indSub(1),3)];
                decoder_events{i}(sessionNum_decoder).replayEvents(j).ratSpeed = Position_Data_sub(indSub(1),5);
                decoder_events{i}(sessionNum_decoder).replayEvents(j).ratHD = Position_Data_sub(indSub(1),4);
                decoder_events{i}(sessionNum_decoder).replayEvents(j).laser_state = laser_state_sub(indSub,2);
            end
        end

    elseif spatialDim == 1 && directionalDecoding == 1

        for i = [5 6]

            %         for i = 1:length(event_boundaries)

            event_boundaries_sub = event_boundaries{i};

            for j = 1:size(event_boundaries_sub,1)

                indSub = event_boundaries_sub(j,1):event_boundaries_sub(j,2);
                % indices of decoder_binDecoding that correspond to the
                % current candidate event

                event_timeBins = timeBins(indSub,:);
                % time bins for decoding the current event:
                % event_timeBins(:,2)-event_timeBins(:,1) is equal to your time
                % window for decoding. Diff(event_timeBins(:,1)) is equal to
                % your window slide.
                event_timepoints = mean(event_timeBins,2);
                % event_timepoints is the midpoint of each window/the mean of
                % event_timeBins

                if i == 5 % if using ripple filtered events you may also want to get the indices of the portion of the replay that was selected.

                    indSub_cropped = cropped_event_boundaries{5}(j,1):cropped_event_boundaries{5}(j,2);
                    if isnan(indSub_cropped)
                        event_timeBins_cropped = [nan nan];
                    else
                        event_timeBins_cropped = timeBins(indSub_cropped,:);
                    end
                elseif i == 6 % if using spikefiltered events you may also want to get the indices of the portion of the replay that was selected.
                    indSub_cropped = cropped_event_boundaries{6}(j,1):cropped_event_boundaries{6}(j,2);
                    if isnan(indSub_cropped)
                        event_timeBins_cropped = [nan nan];
                    else
                        event_timeBins_cropped = timeBins(indSub_cropped,:);
                    end

                else
                    indSub_cropped = event_boundaries(j,1):event_boundaries(j,2);
                    event_timeBins_cropped = timeBins(indSub_cropped,:);
                end

                if i == 4
                    % for events pulled out through John's filtering
                    % method, we need to know which map the filtered event
                    % came from. Sometimes the method will pull out events
                    % in both maps at the same times- we only want to
                    % consider an event in map1 if the majority of the
                    % posterior during that time period came from map 1,
                    % and vice versa. we will also remove the bins that
                    % John's method eliminated (due to high posterior
                    % spread or velocity).

                    filtering_map = event_boundaries_sub(j,3);
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).filtering_map = filtering_map;

                    if filtering_map == 1
                        indNaN = find(isnan(candidateEvents.filtering(sessionNum_decoder).x_NaN_direction_1(indSub,2)));
                        replay = candidateEvents.filtering(sessionNum_decoder).x_direction_1(indSub,2:3);
                        replay_NaN = candidateEvents.filtering(sessionNum_decoder).x_NaN_direction_1(indSub,2:3);
                        replay_NaNremoved = replay_NaN; replay_NaNremoved(indNaN,:) = [];
                        times_NaNremoved = candidateEvents.filtering(sessionNum_decoder).x_direction_1(indSub,1); times_NaNremoved(indNaN) = [];
                    elseif filtering_map == 2
                        indNaN = find(isnan(candidateEvents.filtering(sessionNum_decoder).x_NaN_direction_2(indSub,2)));
                        replay = candidateEvents.filtering(sessionNum_decoder).x_direction_2(indSub,2:3);
                        replay_NaN = candidateEvents.filtering(sessionNum_decoder).x_NaN_direction_2(indSub,2:3);
                        replay_NaNremoved = replay_NaN; replay_NaNremoved(indNaN,:) = [];
                        times_NaNremoved = candidateEvents.filtering(sessionNum_decoder).x_direction_2(indSub,1); times_NaNremoved(indNaN) = [];
                    end
                elseif i == 5 
                    filtering_map = candidateEvents.filtered_ripple_events(j).bestMap;
                    indNaN = find(isnan(candidateEvents.filtered_ripple_events(j).replay_NaN(:,1)));
                elseif  i == 6
                    filtering_map = candidateEvents.filtered_spike_events(j).bestMap;
                    indNaN = find(isnan(candidateEvents.filtered_spike_events(j).replay_NaN(:,1)));
                else
                    indNaN = [];
                end

                %determine the spike density and ripple power during this
                %event.
                %TODO: decide whether this should be for the entire event,
                %or just for the selected 'replay' portion.
                timePoints = [event_timeBins(1,1) event_timeBins(end,2)];
                decoder_events{i}(sessionNum_decoder).replayEvents(j).timePoints_og = timePoints;


                timePoints_cropped = [event_timeBins_cropped(1,1) event_timeBins_cropped(end,2)];


                zscored_LFP_sub = compute_dataTemporalConcatenation(zscored_LFP,timePoints);
                replay_ripple_power = nanmax(zscored_LFP_sub(:,end));

                zscored_LFP_sub_cropped = compute_dataTemporalConcatenation(zscored_LFP,timePoints_cropped);
                replay_ripple_power_cropped = nanmax(zscored_LFP_sub_cropped(:,end));

                zscored_spikeDensity_sub = compute_dataTemporalConcatenation(zscored_spikeDensity,timePoints);
                replay_spikeDensity_power = nanmax(zscored_spikeDensity_sub(:,end));

                zscored_spikeDensity_sub_cropped = compute_dataTemporalConcatenation(zscored_spikeDensity,timePoints_cropped);
                replay_spikeDensity_power_cropped = nanmax(zscored_spikeDensity_sub_cropped(:,end));

                zscored_spikeDensity_sub_pyr = compute_dataTemporalConcatenation(zscored_spikeDensity_pyr,timePoints);
                replay_spikeDensity_power_pyr = nanmax(zscored_spikeDensity_sub_pyr(:,end));

                zscored_spikeDensity_sub_pyr_cropped = compute_dataTemporalConcatenation(zscored_spikeDensity_pyr,timePoints_cropped);
                replay_spikeDensity_power_pyr_cropped = nanmax(zscored_spikeDensity_sub_pyr_cropped(:,end));

                %                 figure()
                %                 plot(zscored_LFP_sub(:,1),zscored_LFP_sub(:,2))
                %                 hold on
                %                 plot(zscored_LFP_sub(:,1),zscored_LFP_sub(:,4))


                %if the replay occured in between segments,
                %replay_ripple_power and replay_spikeDensity power will not
                %be computed- set to nan.
                if isempty(replay_ripple_power)
                    replay_ripple_power = nan;
                end
                if isempty(replay_ripple_power_cropped)
                    replay_ripple_power_cropped = nan;
                end
                if isempty(replay_spikeDensity_power)
                    replay_spikeDensity_power = nan;
                end
                if isempty(replay_spikeDensity_power_cropped)
                    replay_spikeDensity_power_cropped = nan;
                end
                if isempty(replay_spikeDensity_power_pyr)
                    replay_spikeDensity_power_pyr = nan;
                end
                if isempty(replay_spikeDensity_power_pyr_cropped)
                    replay_spikeDensity_power_pyr_cropped = nan;
                end

                %determine which pass this CE occured on:
                [~, pass_ind]= min(abs(Pass_Info_full(:,1)-event_timepoints(1)));
                [~, pos_ind] = min(abs(Position_Data_full(:,1)-event_timepoints(1)));



                %check whether the rat was truly drinking at the time of
                %this event:

                drinking = sum(event_timepoints(1) > Reward_Epoch_Time_Boundaries_drinking(:,1) ...
                    & event_timepoints(1) < Reward_Epoch_Time_Boundaries_drinking(:,2));

                decoder_events{i}(sessionNum_decoder).replayEvents(j).drinking = drinking;




                %compute the time of the replay since the start of the
                %session
                [time_in_session, session_ind] = min(event_timepoints(1) - segment_start_times(event_timepoints(1) - segment_start_times > 0));
                time_in_session = time_in_session/spikeSampRate; % sec
                decoder_events{i}(sessionNum_decoder).replayEvents(j).time_in_session = time_in_session;

                %compute the time of the replay since the last laser onset
                %and the last laser offset.
                laser_state_transitions = diff([laser_state_sub(1,2); laser_state_sub(:,2)]);
                laser_onset_times = laser_state_sub(laser_state_transitions == 1,1);
                laser_offset_times = laser_state_sub(laser_state_transitions == -1,1);
                replay_start = event_timepoints(1);

                time_since_laser_transition = event_timepoints(1)-laser_onset_times;
                [t,~] = min(time_since_laser_transition(time_since_laser_transition > 0));
                if(isempty(t))
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).time_since_last_laser_onset = nan;
                else
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).time_since_last_laser_onset = t./spikeSampRate;
                end

                time_since_laser_transition = event_timepoints(1)-laser_offset_times;
                [t,~] = min(time_since_laser_transition(time_since_laser_transition > 0));
                if(isempty(t))
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).time_since_last_laser_offset = nan;
                else
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).time_since_last_laser_offset = t./spikeSampRate;
                end


                % determine if the replay occured inside one of the reward
                % zones. If so, also label which one it was.
                ratPos = [Position_Data_sub_scaled(indSub(1),2),Position_Data_sub_scaled(indSub(1),3)];

                % distance from each reward_end:
                decoder_events{i}(sessionNum_decoder).replayEvents(j).distance_from_end_1 = abs(Experiment_Information.reward_locations(1) - ratPos(1)*binSize);
                decoder_events{i}(sessionNum_decoder).replayEvents(j).distance_from_end_2 = abs(Experiment_Information.reward_locations(2)- ratPos(1)*binSize);

                % which reward zone he's closer to:
                [~,reward_zone] = min([abs(Position_Data_sub_scaled(indSub(1),2) - num_bins_in_reward_zone),abs(Position_Data_sub_scaled(indSub(1),2) - (numSpatialBins(2) - num_bins_in_reward_zone))]);

                % binary categorization of whether he was in a reward zone
                % or not:
                if ratPos(1) <= num_bins_in_reward_zone || ratPos(1) >= numSpatialBins(2) - num_bins_in_reward_zone
                    in_reward_zone = 1;
                else
                    in_reward_zone = 0;
                end

                decoder_events{i}(sessionNum_decoder).replayEvents(j).in_reward_zone = in_reward_zone;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).reward_zone = reward_zone;

                % determine if this replay occured during a run:
                if sum(ismember(Experiment_Information.Segments(session_ind).Flags,(11)) > 0)
                    run_session = 1;
                else
                    run_session = 0;
                end

                if run_session == 1
                    % Classify event as 'engaged, entry (1)','engaged, exit' (-1), or
                    % 'disengaged' (0);

                    decoder_events{i}(sessionNum_decoder).replayEvents(j).ratPos = ratPos;


                    % time_from_transition = min([time_since_nearest_entry time_until_nearest_exit])

                    [time_from_transition, nearest_transitions_ind] = min(abs(Pass_Transitions(:,1) - event_timepoints(1)));
                    time_from_transition = time_from_transition/spikeSampRate;
                    if time_from_transition > 5
                        engagement_classification = 0;
                    elseif time_from_transition <= 5 && Pass_Transitions(nearest_transitions_ind,2) < 0
                        engagement_classification = -1;
                    elseif time_from_transition <= 5 && Pass_Transitions(nearest_transitions_ind,2) > 0
                        engagement_classification = 1;
                    end
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).engagement_classification = engagement_classification;

                    % Regardless of whether the event is 'engaged' or
                    % 'disengaged' we might want to know if the replay occured
                    % closer to a reward zone entry or closer to a reward zone
                    % exit.
                    % -1 = exit from left end, 1 = entry to left end, -2 = exit from right end, 2 = entry to right end
                    entry_exit_classification = Pass_Transitions(nearest_transitions_ind,2);
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).entry_exit_classification = entry_exit_classification;

                    % Compute the time of the replay since the last
                    % entry, and the time until the next exit.

                    reward_entry_times = Pass_Transitions(Pass_Transitions(:,2)>0,:);
                    time_since_reward_zone_entry = event_timepoints(1) - reward_entry_times(:,1);
                    time_since_reward_zone_entry(time_since_reward_zone_entry < 0) = inf;
                    [time_since_reward_zone_entry, nearest_entry_ind] = min(time_since_reward_zone_entry);
                    time_since_reward_zone_entry = time_since_reward_zone_entry./spikeSampRate;

                    drink_onsets = Reward_Epoch_Time_Boundaries_drinking(:,1);
                    time_since_drinks = event_timepoints(1) - drink_onsets;
                    time_since_drinks(time_since_drinks < 0) = inf;
                    [time_since_drink_onset,~] = min(time_since_drinks);
                    time_since_drink_onset(time_since_drink_onset == inf) = nan;
                    time_since_drink_onset = time_since_drink_onset./spikeSampRate;


                    differences = Pass_Transitions(Pass_Transitions(:,2)<0,1) - event_timepoints(1);
                    differences(differences < 0) = inf;
                    [time_till_reward_zone_exit, nearest_exit_ind] = min(differences);
                    time_till_reward_zone_exit = time_till_reward_zone_exit/spikeSampRate;


                    decoder_events{i}(sessionNum_decoder).replayEvents(j).time_since_reward_zone_entry = time_since_reward_zone_entry;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).time_till_reward_zone_exit = time_till_reward_zone_exit;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).time_since_drink_onset = time_since_drink_onset;

                    nearest_entry = reward_zone_entries(nearest_entry_ind);
                    entry_exits_diff = reward_zone_exits(:,1) - nearest_entry;
                    entry_exits_diff(entry_exits_diff < 0) = inf;
                    [val, nearest_exit_ind] = min(entry_exits_diff);
                    nearest_exit = reward_zone_exits(nearest_exit_ind);

                    laser_state_during_stopping_period = compute_dataTemporalConcatenation(laser_state_sub,[nearest_entry nearest_exit]);
                    running_speed_during_stopping_period = compute_dataTemporalConcatenation(Position_Data,[nearest_entry nearest_exit]);
                    running_speed_during_stopping_period = running_speed_during_stopping_period(:,5);

                    if ismember(1,Experiment_Information.Segments(session_ind).Flags)
                        if sum(laser_state_during_stopping_period(:,2)) > 0
                            intended_laser_state_during_stopping_period = 1;
                        else
                            intended_laser_state_during_stopping_period = 0;
                        end
                    elseif ismember(2,Experiment_Information.Segments(session_ind).Flags)
                        intended_laser_state_during_stopping_period = mode(laser_state_during_stopping_period(running_speed_during_stopping_period < speedThr,2));
                    else
                        intended_laser_state_during_stopping_period = 0;
                    end


                    if sum(laser_state_sub(indSub,2)) > 0
                        event_laser_state = 1;
                    else
                        event_laser_state = 0;
                    end

                    if event_laser_state == intended_laser_state_during_stopping_period
                        decoder_events{i}(sessionNum_decoder).replayEvents(j).laser_state_matches_intended_laser_state_for_stopping_period = 1;
                    else
                        decoder_events{i}(sessionNum_decoder).replayEvents(j).laser_state_matches_intended_laser_state_for_stopping_period = 0;
                    end


                else
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).engagement_classification = nan;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).entry_exit_classification = nan;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).time_since_drink_onset = nan;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).time_since_reward_zone_entry = nan;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).time_till_reward_zone_exit = nan;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).laser_state_matches_intended_laser_state_for_stopping_period = nan;
                end

                % Re-do the decoding on the candidate event to
                % pull out some new metrics (can add these to
                % load_binDecoding_cm in the future.



                event_timeBins_NaN_removed = event_timeBins; event_timeBins_NaN_removed(indNaN,:) = [];
                if ~isempty(event_timeBins_NaN_removed)
                    numSpks = load_numSpks_timeBins(event_timeBins_NaN_removed ,clusters,(1:length(clusters)),decoder_binDecoding(1).shiftSizeDecoding*spikeSampRate,decoder_binDecoding(1).windowSizeDecoding*spikeSampRate);


                    % find the percentage of cells that fired at least 1 spike in the replay:
                    numSpks_sum = sum(numSpks,2);
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).fraction_of_cells_participating = sum(numSpks_sum>0)/size(numSpks,1);

                    % add in the number of cells that fired at least 1
                    % spike in the replay:
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).number_of_cells_participating = sum(numSpks_sum>0);

                    % add in the average number of spikes emiited by each
                    % participating excitatory cell, and the average firing
                    % rate (spikes emitted/replay length)
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_num_spikes_all_excitatory_cells = nanmean(numSpks_sum(excitatory_inds));
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_fr_all_excitatory_cells = nanmean(numSpks_sum(excitatory_inds))/(size(numSpks,2)*decoder_binDecoding(1).shiftSizeDecoding);
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_num_spikes_all_inhibitory_cells = nanmean(numSpks_sum(inhibitory_inds));
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_fr_all_inhibitory_cells = nanmean(numSpks_sum(inhibitory_inds))/(size(numSpks,2)*decoder_binDecoding(1).shiftSizeDecoding);
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_num_spikes_all_bimodal_cells = nanmean(numSpks_sum(bimodal_inds));
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_fr_all_bimodal_cells = nanmean(numSpks_sum(bimodal_inds))/(size(numSpks,2)*decoder_binDecoding(1).shiftSizeDecoding);
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_num_spikes_all_unimodal_cells = nanmean(numSpks_sum(unimodal_inds));
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_fr_all_unimodal_cells = nanmean(numSpks_sum(unimodal_inds))/(size(numSpks,2)*decoder_binDecoding(1).shiftSizeDecoding);
                   

                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_num_spikes_participating_excitatory_cells = nanmean(numSpks_sum(intersect(find(numSpks_sum>0),excitatory_inds)));
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_fr_participating_excitatory_cells = nanmean(numSpks_sum(intersect(find(numSpks_sum>0),excitatory_inds)))/(size(numSpks,2)*decoder_binDecoding(1).shiftSizeDecoding);
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_num_spikes_participating_inhibitory_cells = nanmean(numSpks_sum(intersect(find(numSpks_sum>0),inhibitory_inds)));
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_fr_participating_inhibitory_cells = nanmean(numSpks_sum(intersect(find(numSpks_sum>0),inhibitory_inds)))/(size(numSpks,2)*decoder_binDecoding(1).shiftSizeDecoding);
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_num_spikes_participating_bimodal_cells = nanmean(numSpks_sum(intersect(find(numSpks_sum>0),bimodal_inds)));
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_fr_participating_bimodal_cells = nanmean(numSpks_sum(intersect(find(numSpks_sum>0),bimodal_inds)))/(size(numSpks,2)*decoder_binDecoding(1).shiftSizeDecoding);
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_num_spikes_participating_unimodal_cells = nanmean(numSpks_sum(intersect(find(numSpks_sum>0),unimodal_inds)));
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_fr_participating_unimodal_cells = nanmean(numSpks_sum(intersect(find(numSpks_sum>0),unimodal_inds)))/(size(numSpks,2)*decoder_binDecoding(1).shiftSizeDecoding);



                    % add in the number of cells in the session
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).total_num_cells_in_session = size(numSpks,1);

                    participating_clusters = numSpks_sum > 0;
                    participating_bimodal_cells = participating_clusters(bimodal_inds);
                    participating_unimodal_cells = participating_clusters(unimodal_inds);
                    participating_good_excitatory_cells = participating_clusters(good_cluster_inds);

                    if isempty(good_cluster_inds)
                        fraction_of_good_excitatory_cells_participating = nan;
                        number_of_good_excitatory_cells_participating = nan;
                        number_good_excitatory_cell_in_session = nan;
                    else
                        fraction_of_good_excitatory_cells_participating = sum(participating_good_excitatory_cells)/length(good_cluster_inds);
                        number_of_good_excitatory_cells_participating = sum(participating_good_excitatory_cells);
                        number_of_good_excitatory_cells_in_session = length(good_cluster_inds);
                    end
                    if isempty(bimodal_inds)
                        fraction_of_bimodal_cells_participating = nan;
                    else
                        fraction_of_bimodal_cells_participating = sum(participating_bimodal_cells)/length(bimodal_inds);
                    end

                    if isempty(unimodal_inds)
                        fraction_of_unimodal_cells_participating = nan;
                    else
                        fraction_of_unimodal_cells_participating = sum(participating_unimodal_cells)/length(unimodal_inds);
                    end



                else

                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_num_spikes_all_excitatory_cells = nan;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_fr_all_excitatory_cells = nan;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_num_spikes_all_inhibitory_cells = nan;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_fr_all_inhibitory_cells = nan;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_num_spikes_all_bimodal_cells = nan;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_fr_all_bimodal_cells = nan;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_num_spikes_all_unimodal_cells = nan;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_fr_all_unimodal_cells = nan;
                   

                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_num_spikes_participating_excitatory_cells  = nan;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_fr_participating_excitatory_cells = nan;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_num_spikes_participating_inhibitory_cells = nan;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_fr_participating_inhibitory_cells = nan;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_num_spikes_participating_bimodal_cells = nan;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_fr_participating_bimodal_cells = nan;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_num_spikes_participating_unimodal_cells = nan;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).mean_fr_participating_unimodal_cells = nan;
















                    decoder_events{i}(sessionNum_decoder).replayEvents(j).fraction_of_cells_participating = NaN;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).number_of_cells_participating = NaN;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).total_num_cells_in_session = NaN;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).number_of_good_excitatory_cells_in_session = NaN;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).number_of_good_excitatory_cells_participating = NaN;

                    fraction_of_bimodal_cells_participating = NaN;
                    fraction_of_unimodal_cells_participating = NaN;
                    fraction_of_good_excitatory_cells_participating = NaN;
                    number_of_good_excitatory_cells_participating = NaN;
                    number_of_good_excitatory_cells_in_session = NaN;


                end

                decoder_events{i}(sessionNum_decoder).replayEvents(j).fraction_of_good_excitatory_cells_participating = fraction_of_good_excitatory_cells_participating;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).number_of_good_excitatory_cells_participating = number_of_good_excitatory_cells_participating;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).number_of_good_excitatory_cells_in_session = number_of_good_excitatory_cells_in_session;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).fraction_of_bimodal_cells_participating = fraction_of_bimodal_cells_participating;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).fraction_of_unimodal_cells_participating = fraction_of_unimodal_cells_participating;
                %check if the posterior should be re-computed:
                %Update 8/3/22: I don't want to recompute the posterior-
                %it's already saved and also, computing numSpks on the
                %whole session versus just times of interest yields
                %slightly different results, making it hard to compare

                posterior = decoder_binDecoding(sessionNum_decoder).posterior(indSub,:);
                num_spatial_bins = size(posterior,2);
                left_map = posterior(:,1:num_spatial_bins/2);
                right_map = posterior(:,num_spatial_bins/2+1:num_spatial_bins);
                % Mask out the rat's current location, and calculate how
                % much remote posterior is in the left map versus the right
                % map.
                bins_to_mask = find(abs(1:num_spatial_bins-ratPos(1)) < num_position_bins_to_mask_out_local_posterior);

                left_map(:,bins_to_mask) = nan;
                right_map(:,bins_to_mask) = nan;


                decoder_events{i}(sessionNum_decoder).replayEvents(j).percent_non_local_posterior_in_left_map = nansum(nansum(left_map))/nansum(nansum(left_map+right_map));
                decoder_events{i}(sessionNum_decoder).replayEvents(j).percent_non_local_posterior_in_right_map = nansum(nansum(right_map))/nansum(nansum(left_map+right_map));





                %                 if compute_posterior == 1
                %                     [posterior,posteriorCOM,posteriorSpread,posteriorPeak] =  compute_BayesianDecoding_fast_1D(numSpks,rateMap_smoothed,numSpatialBins,decoder_binDecoding(1).windowSizeDecoding,return_fullPosterior,directionalDecoding);
                %                 else
                %                     posterior =  [decoder_events{i}(sessionNum_decoder).replayEvents(j).posterior_left decoder_events{i}(sessionNum_decoder).replayEvents(j).posterior_right];
                %                 end
                % Metrics already computed:
                % posterior, posteriorCOM, posteriorPeak
                if i == 4 || i == 5 || i == 6 % bottom up filtering methods

                    eventMetrics = load_candidateEventMetricsDirectionalDecoding(posterior,event_timepoints,filtering_map,indNaN);
                else
                    eventMetrics = load_candidateEventMetricsDirectionalDecoding(posterior,event_timepoints,[],[]);
                end


                % Finally, remove the NaN's when reporting the timeBins and
                % timePoints associated with this replay.

                indSub_orig =indSub;
                indSub(indNaN) = [];

                % duration_og and duration will only differ for filtered
                % spike or filtered ripple events. duration_og is the
                % duration of the underlying candidate event (i.e., spike
                % density event). Duration is the length of the replay
                % segment that was pulled out of the candidate event.

                if ~isempty(indSub)
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).timeBins =  timeBins(indSub,:);
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).timePoints_og =  [decoding_time_bin_centers(indSub_orig(1),1) decoding_time_bin_centers(indSub_orig(end),1)];
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).timePoints =  [decoding_time_bin_centers(indSub(1),1) decoding_time_bin_centers(indSub(end),1)];
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).duration = (decoding_time_bin_centers(indSub(end),1)-decoding_time_bin_centers(indSub(1),1))/spikeSampRate;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).duration_og = (decoding_time_bin_centers(indSub_orig(end),1)-decoding_time_bin_centers(indSub_orig(1),1))/spikeSampRate;
                    
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).ratLoc = [Position_Data_sub_scaled(indSub(1),2),Position_Data_sub_scaled(indSub(1),3)];
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).ratSpeed = Position_Data_sub(indSub(1),5);
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).ratHD = Position_Data_sub(indSub(1),4);
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).laser_state = laser_state_sub(indSub,2);
                else
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).timeBins = timeBins(indSub_orig,:); % we might still want to know the time that this candidate event happened even if there was no good content.
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).timePoints = [decoding_time_bin_centers(indSub_orig(1),1) decoding_time_bin_centers(indSub_orig(end),1)];
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).timePoints_og = [decoding_time_bin_centers(indSub_orig(1),1) decoding_time_bin_centers(indSub_orig(end),1)];
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).duration_og = (decoding_time_bin_centers(indSub_orig(end),1)-decoding_time_bin_centers(indSub_orig(1),1))/spikeSampRate;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).duration = nan;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).ratLoc = [Position_Data_sub_scaled(indSub_orig(1),2),Position_Data_sub_scaled(indSub_orig(1),3)];
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).ratSpeed = Position_Data_sub(indSub_orig(1),5);
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).ratHD = Position_Data_sub(indSub_orig(1),4);
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).laser_state = laser_state_sub(indSub_orig,2);
                end


                decoder_events{i}(sessionNum_decoder).replayEvents(j).replay_ripple_power = replay_ripple_power;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).replay_ripple_power_cropped = replay_ripple_power_cropped;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).replay_spikeDensity_power = replay_spikeDensity_power;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).replay_spikeDensity_power_cropped = replay_spikeDensity_power_cropped;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).replay_spikeDensity_power_pyr = replay_spikeDensity_power_pyr;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).replay_spikeDensity_power_pyr_cropped = replay_spikeDensity_power_pyr_cropped;
                        

                decoder_events{i}(sessionNum_decoder).replayEvents(j).posterior_left = eventMetrics.posterior_left;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).posterior_right = eventMetrics.posterior_right;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).weighted_r = eventMetrics.weighted_r;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).weighted_r = eventMetrics.weighted_r;

                %                 decoder_events{i}(sessionNum_decoder).replayEvents(j).replay_score = eventMetrics.replay_score;
                %                 decoder_events{i}(sessionNum_decoder).replayEvents(j).replay_slope = eventMetrics.replay_slope;

                decoder_events{i}(sessionNum_decoder).replayEvents(j).com = eventMetrics.replay_sequence_com;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).peak = eventMetrics.replay_sequence_peak;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).sharpness = eventMetrics.sharpness;

                decoder_events{i}(sessionNum_decoder).replayEvents(j).ave_sharpness_at_peak = eventMetrics.ave_sharpness_at_peak;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).total_posterior_left = eventMetrics.total_posterior_left;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).total_posterior_right = eventMetrics.total_posterior_right;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).best_map = eventMetrics.best_map;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).range = eventMetrics.posterior_range_peak;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).range_normalized = eventMetrics.posterior_range_normalized;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).range_over_time = eventMetrics.best_map_distance_over_time;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).coverage_wu = eventMetrics.coverage_wu;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).cum_coverage = eventMetrics.cum_coverage;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).max_jump_distance = eventMetrics.max_jump_distance_peak;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).ave_jump_distance = eventMetrics.ave_jump_distance_peak;

                [startsLocal, start_distance_from_rat, endsLocal, end_distance_from_rat] = check_replay_for_local_start_or_stop(eventMetrics.replay_sequence_peak, ratPos, spatialDim, num_bins_to_check_for_local_start, bin_thr_for_local_start);
                decoder_events{i}(sessionNum_decoder).replayEvents(j).startsLocal = startsLocal;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).endsLocal = endsLocal;
                decoder_events{i}(sessionNum_decoder).replayEvents(j).start_distance_from_rat = start_distance_from_rat; %in bins currently
                decoder_events{i}(sessionNum_decoder).replayEvents(j).end_distance_from_rat = end_distance_from_rat; %in bins currently
                decoder_events{i}(sessionNum_decoder).replayEvents(j).start_distance_from_rat = start_distance_from_rat; %in bins currently
                decoder_events{i}(sessionNum_decoder).replayEvents(j).end_distance_from_rat = end_distance_from_rat; %in bins currently






                % add forward/reverse designation
                % 1 = forward; 2 = reverse; nan means the replay was
                % static, had a weighted correlation of 0. we should remove
                % these events.
                if (eventMetrics.best_map == 1 && eventMetrics.weighted_r < 0) || (eventMetrics.best_map == 2 && eventMetrics.weighted_r > 0)
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).direction = 1;
                elseif (eventMetrics.best_map == 1 && eventMetrics.weighted_r > 0) || (eventMetrics.best_map == 2 && eventMetrics.weighted_r < 0)
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).direction = 2;
                else
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).direction = nan;
                end


                % add congruent designation
                if (eventMetrics.best_map == 1 && eventMetrics.weighted_r < 0 && reward_zone == 2) || (eventMetrics.best_map == 2 && eventMetrics.weighted_r > 0 && reward_zone == 1)
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).forward_congruent_with_rat_location = 1;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).reverse_congruent_with_rat_location = 0;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).congruent_with_rat_location=1;
                elseif (eventMetrics.best_map == 1 && eventMetrics.weighted_r > 0 && reward_zone == 1) || (eventMetrics.best_map == 2 && eventMetrics.weighted_r < 0 && reward_zone == 2)
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).forward_congruent_with_rat_location = 0;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).reverse_congruent_with_rat_location = 1;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).congruent_with_rat_location=1;
                else
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).forward_congruent_with_rat_location = 0;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).reverse_congruent_with_rat_location = 0;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).congruent_with_rat_location=0;
                end

                % add incongruent designation
                if (eventMetrics.best_map == 1 && eventMetrics.weighted_r < 0 && reward_zone == 1) || (eventMetrics.best_map == 2 && eventMetrics.weighted_r > 0 && reward_zone == 2)
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).forward_incongruent_with_rat_location = 1;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).reverse_incongruent_with_rat_location = 0;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).incongruent_with_rat_location=1;
                elseif (eventMetrics.best_map == 1 && eventMetrics.weighted_r > 0 && reward_zone == 2) || (eventMetrics.best_map == 2 && eventMetrics.weighted_r < 0 && reward_zone == 1)
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).forward_incongruent_with_rat_location = 0;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).reverse_incongruent_with_rat_location = 1;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).incongruent_with_rat_location=1;
                else
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).forward_incongruent_with_rat_location = 0;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).reverse_incongruent_with_rat_location = 0;
                    decoder_events{i}(sessionNum_decoder).replayEvents(j).incongruent_with_rat_location=0;
                end
                display(['finished ' num2str(j) '/' num2str(length(event_boundaries_sub))])
                drawnow
            end
        end

    end
end


try
    save('decoder_candidateEvents.mat','decoder_events')
catch
    save('decoder_candidateEvents.mat','decoder_events','-v7.3')
end

