
function eventMetrics = load_behavioral_measures_during_replay(params,event_timePoints,event_inds,laser_state,Experiment_Information,Pass_Info,Position_Data,Position_Data_scaled,...
    segment_start_times,Pass_Transitions,Reward_Epoch_Time_Boundaries_drinking,Reward_Epoch_Time_Boundaries_speed_thresholded)

% params.end_zone_size = 30; % cm
% params.speedThr = 5; % cm/s
% params.spikeSampRate = 30000;

% Determine which pass (2 passes=1 lap) this event occured on.
pass = compute_dataInterpolation(Pass_Info,event_timePoints(1),[]);
eventMetrics.pass_number = pass(2);


% determine whether the laser was on or off during this event
if sum(laser_state(event_inds,2)) > 0
    eventMetrics.laser_state = 1;
else
    eventMetrics.laser_state = 0;
end

%check whether the rat was truly drinking at the time of this event:
eventMetrics.drinking = sum(event_timePoints(1) > Reward_Epoch_Time_Boundaries_drinking(:,1) ...
    & event_timePoints(1) < Reward_Epoch_Time_Boundaries_drinking(:,2));

%compute the time of the replay since the start of the session
%session
[time_in_session, session_ind] = min(event_timePoints(1) - segment_start_times(event_timePoints(1) - segment_start_times >= 0));
time_in_session = time_in_session/params.spikeSampRate; % sec
eventMetrics.time_in_session = time_in_session;

%compute the time of the replay since the last laser onset
%and the last laser offset.
laser_state_transitions = diff([laser_state(1,2); laser_state(:,2)]);
laser_onset_times = laser_state(laser_state_transitions == 1,1);
laser_offset_times = laser_state(laser_state_transitions == -1,1);

time_since_laser_transition = event_timePoints(1)-laser_onset_times;
[time_since_laser_transition,~] = min(time_since_laser_transition(time_since_laser_transition > 0));
if(isempty(time_since_laser_transition))
    eventMetrics.time_since_last_laser_onset = nan;
else
    eventMetrics.time_since_last_laser_onset = time_since_laser_transition./params.spikeSampRate;
end

time_since_laser_transition = event_timePoints(1)-laser_offset_times;
[time_since_laser_transition,~] = min(time_since_laser_transition(time_since_laser_transition > 0));
if(isempty(time_since_laser_transition))
    eventMetrics.time_since_last_laser_offset = nan;
else
    eventMetrics.time_since_last_laser_offset = time_since_laser_transition./params.spikeSampRate;
end

% Rat's position, head direction, and speed during the candidate event
ratPos = [Position_Data(event_inds(1),2),Position_Data(event_inds(1),3)];
ratPos_scaled = [Position_Data_scaled(event_inds(1),2),Position_Data_scaled(event_inds(1),3)];
eventMetrics.ratPos = ratPos_scaled;
eventMetrics.ratSpeed = Position_Data(event_inds(1),5);
eventMetrics.ratHD = Position_Data(event_inds(1),4);

% Compute the average head direction for each reward zone. Then calculate
% how far off the rat's actual head direction is.

mean_hd_left_end = circ_mean(Position_Data(Position_Data(:,2)>= 0 & Position_Data(:,2)<10 & Position_Data(:,5)<5,4));
mean_hd_left_end(mean_hd_left_end<0)=mean_hd_left_end+2*pi;

mean_hd_right_end = circ_mean(Position_Data(Position_Data(:,2)<= Experiment_Information.maze_size & Position_Data(:,2)>Experiment_Information.maze_size - 10 & Position_Data(:,5)<5,4));
mean_hd_right_end(mean_hd_right_end<0)=mean_hd_right_end+2*pi;
 
eventMetrics.deviation_from_ave_hd_at_left_end = circ_dist(mean_hd_left_end,eventMetrics.ratHD);
eventMetrics.deviation_from_ave_hd_at_right_end = circ_dist(mean_hd_right_end,eventMetrics.ratHD);

% Distance between rat's position during the candidate event and
% each end zone.
reward_locations = [0 Experiment_Information.maze_size];
eventMetrics.distance_from_end_1 = abs(reward_locations(1) - ratPos(1));
eventMetrics.distance_from_end_2 = abs(reward_locations(2) - ratPos(1));
[~,eventMetrics.reward_zone] = min([abs(reward_locations(1) - ratPos(1)),abs(reward_locations(2) - ratPos(1))]);

if ratPos(1) <= reward_locations(1) + params.end_zone_size || ratPos(1) >= reward_locations(2)-params.end_zone_size
    eventMetrics.in_reward_zone = 1;
else
    eventMetrics.in_reward_zone = 0;
end

% determine if this replay occured during a run:
if sum(ismember(Experiment_Information.Segments(session_ind).Flags,(11)) > 0)
    run_session = 1;
else
    run_session = 0;
end

if run_session == 1


    % Determine whether the replay occured
    % closer to a reward zone entry or closer to a reward zone
    % exit.
    % -1 = exit from left end, 1 = entry to left end, -2 = exit from right end, 2 = entry to right end
    [~, nearest_transitions_ind] = min(abs(Pass_Transitions(:,1) - event_timePoints(1)));
    entry_exit_classification = Pass_Transitions(nearest_transitions_ind,2);
    eventMetrics.entry_exit_classification = entry_exit_classification;

    % Compute the time of the replay since the last
    % entry, and the time until the next exit.

    reward_entry_times = Reward_Epoch_Time_Boundaries_speed_thresholded(:,1);
    time_since_reward_zone_entry = event_timePoints(1) - reward_entry_times;
    time_since_reward_zone_entry(time_since_reward_zone_entry < 0) = inf;
    [time_since_reward_zone_entry, nearest_entry_ind] = min(time_since_reward_zone_entry);
    time_since_reward_zone_entry = time_since_reward_zone_entry./params.spikeSampRate;

    drink_onsets = Reward_Epoch_Time_Boundaries_drinking(:,1);
    time_since_drinks = event_timePoints(1) - drink_onsets;
    time_since_drinks(time_since_drinks < 0) = inf;
    [time_since_drink_onset,~] = min(time_since_drinks);
    time_since_drink_onset(time_since_drink_onset == inf) = nan;
    time_since_drink_onset = time_since_drink_onset./params.spikeSampRate;

    reward_exit_times = Reward_Epoch_Time_Boundaries_speed_thresholded(:,2);
    time_till_reward_zone_exit = reward_exit_times(nearest_entry_ind)-event_timePoints(1);
    time_till_reward_zone_exit(time_till_reward_zone_exit < 0) = nan; % rare case where the rat stopped in the middle of the track, so the replay is actually after the last exit.
    time_till_reward_zone_exit = time_till_reward_zone_exit/params.spikeSampRate;

    % Classify event as 'engaged, entry (1)','engaged, exit' (-1), or
    % 'disengaged' (0);

    if time_since_reward_zone_entry > 5 && time_till_reward_zone_exit > 5
        eventMetrics.engagement_classification = 0;
    elseif time_since_reward_zone_entry <= 5
        eventMetrics.engagement_classification = 1;
    elseif time_since_reward_zone_entry > 5 && time_till_reward_zone_exit < 5
        eventMetrics.engagement_classification = -1;
    else
        eventMetrics.engagement_classification = nan;
    end

    eventMetrics.time_since_reward_zone_entry = time_since_reward_zone_entry;
    eventMetrics.time_till_reward_zone_exit = time_till_reward_zone_exit;
    eventMetrics.time_since_drink_onset = time_since_drink_onset;

    laser_state_during_stopping_period = compute_dataTemporalConcatenation(laser_state,[reward_entry_times(nearest_entry_ind) reward_exit_times(nearest_entry_ind)]);
    running_speed_during_stopping_period = compute_dataTemporalConcatenation(Position_Data,[reward_entry_times(nearest_entry_ind) reward_exit_times(nearest_entry_ind)]);
    running_speed_during_stopping_period = running_speed_during_stopping_period(:,5);

    if ismember(1,Experiment_Information.Segments(session_ind).Flags)
        if sum(laser_state_during_stopping_period(:,2)) > 0
            intended_laser_state_during_stopping_period = 1;
        else
            intended_laser_state_during_stopping_period = 0;
        end
    elseif ismember(2,Experiment_Information.Segments(session_ind).Flags)
        intended_laser_state_during_stopping_period = mode(laser_state_during_stopping_period(running_speed_during_stopping_period < params.speedThr,2));
    else
        intended_laser_state_during_stopping_period = 0;
    end

    if eventMetrics.laser_state == intended_laser_state_during_stopping_period
        eventMetrics.laser_state_matches_intended_laser_state_for_stopping_period = 1;
    else
        eventMetrics.laser_state_matches_intended_laser_state_for_stopping_period = 0;
    end


else
    eventMetrics.engagement_classification = nan;
    eventMetrics.entry_exit_classification = nan;
    eventMetrics.time_since_drink_onset = nan;
    eventMetrics.time_since_reward_zone_entry = nan;
    eventMetrics.time_till_reward_zone_exit = nan;
    eventMetrics.laser_state_matches_intended_laser_state_for_stopping_period = nan;
end
