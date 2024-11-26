function t = load_replays_from_individual_session(candidate_events_to_plot,use_duration_og,override_coverage_thr,override_weighted_r_thr,override_posterior_diff_thr,sde_amplitude_criterion,ripple_power_criterion,lower_limit_time_since_reward_zone_entry,upper_limit_time_since_reward_zone_entry)

% Outputs: a structure, t, containing metrics for all the candidate events
% in the session.

% Inputs: 
% 1. candidate_events_to_plot: 'ripple_filtered' or 'spike_filtered'
% 2. use_duration_og flag. set to 1 if you want to use the duration of the
% underlying ripple or spike event from which the replay was pulled out.
% Set to 0 if you want to use the duration of the smooth replay segment
% If you are looking to pull out all qualifying spike density events or
% ripple events, set this to 1. If you're looking to pull out replays, set
% to 0.
% 3. Override_coverage. This is the coverage threshold to be considered a
% replay (from 0 to 1). If empty, the coverage threshold will be pulled
% from load_replay_criterion.
% 4. Override_weighted_r_thr. The lower threshold on weighted correlation
% required to be a replay. If empty, the weighted r threshold will be
% pulled from load_replay_criterion.
% 5. override_posterior_diff_thr. If set to zero, all candidate events will
% pass. Otherwise, this is the minimum difference between the total
% posterior in the right map and left map
% 5. sde_amplitude_criterion. The minimum amplitude required to be
% considered a SDE. Set to -inf if considering ripple events, 3 if
% considering spike density events.
% 6. ripple_power_criterion. The minimum ripple power required to be
% considered a ripple event. Set to -inf if considering spike density
% events, 3 if considering ripple events.
% 7. upper_limit_time_since_reward_zone_entry. Set to inf to include all
% events. You can also limit (for example, to the first 10 seconds of a
% stopping period).

% The following could also be made into inputs, but they rarely change:
time_into_session_thr = inf; % set to inf to including everything
end_zone_thr = 1; % set to 0 to include everthing
distance_from_end_tight_threshold = 30; % set to inf to include everthing
engagement_classifications_to_use = [-1 1 0]; % set to [-1,1,0] to include everything
% upper_limit_time_since_reward_zone_entry = inf;
% lower_limit_time_since_reward_zone_entry = 0;
long_replays_only = 1;
num_cells_participating_thr = 5;
fraction_cells_participating_thr = 0;

% posterior_difference_to_use = [{'percent_non_local_posterior_in_left_map_full'},{'percent_non_local_posterior_in_right_map_full'}];
% posterior_difference_to_use = [{'percent_posterior_ends_masked_in_left_map_full'},{'percent_posterior_ends_masked_in_right_map_full'}];
 posterior_difference_to_use = [{'percent_posterior_in_left_map'},{'percent_posterior_in_right_map'}];


load decoder_candidateEvents

[~, duration_thr, max_jump_distance_thr, coverage_thr_percent_of_track, coverage_thr, weighted_r_thr, posterior_in_map_thr] = load_replay_criterion(long_replays_only,override_coverage_thr,override_weighted_r_thr, override_posterior_diff_thr);

binSize = 2;
coverage_thr = coverage_thr/binSize; % coverage threshold needs to be in bins

only_consider_replays_with_intended_laser_state = 1;


if strcmp(candidate_events_to_plot,'laser')
    coverage_thr = 0;
    sde_amplitude_criterion = 0;
    max_jump_distance_thr = 1;
    duration_thr = 0;
    coverage_thr_percent_of_track = 0;
    weighted_r_thr = 0;
    posterior_in_map_thr = 0;
end


load Experiment_Information
load Analysis_Information
load clusters
load Position_Data
load laser_state


if strcmp(candidate_events_to_plot,'ripple')
    decoder_events_sub = decoder_events(1);
elseif strcmp(candidate_events_to_plot,'spike')
    decoder_events_sub = decoder_events(2);
elseif strcmp(candidate_events_to_plot,'laser')
    decoder_events_sub = decoder_events(3);
elseif strcmp(candidate_events_to_plot,'filtered')
    decoder_events_sub = decoder_events(4);
elseif strcmp(candidate_events_to_plot,'ripple_filtered')
    decoder_events_sub = decoder_events(5);
elseif strcmp(candidate_events_to_plot,'spike_filtered')
    decoder_events_sub = decoder_events(6);
end

sessionNum = 1;
session_flags = Experiment_Information.Segments(sessionNum).Flags;
% if this was a sleep session, ignore all the position requirements.
if sum(ismember(session_flags,12)) > 0
    end_zone_thr = 0;
end

% get the indices of the candidate events for the session of interest
t = struct2table(decoder_events_sub.replayEvents);


t.replay_spikeDensity_power = t.ripple_power_sd_metrics.zscored_sd_pyr_125.peak_power_over_fraction_of_event(:,4);
t.replay_ripple_power = t.ripple_power_sd_metrics.zscored_ripple_power.peak_power_over_fraction_of_event(:,4);
% Flag 19 means this was a saline session. change all
% 'laser states' to 0. Flag 18 means this was a CNO session.
% change all 'laser states' to 1.
if sum(ismember(session_flags,19)>0)
    t.laser_state(:) = zeros(height(t),1);
end
if sum(ismember(session_flags,18)>0)
    t.laser_state(:) = ones(height(t),1);
end

ind_session = ismember(t.timePoints_og,compute_dataTemporalConcatenation(t.timePoints_og,Experiment_Information.Segments(sessionNum).Times),'rows');

% Add columns to the table distinguishing if the replay events met
% various criterion
t.meets_engagement_requirement = zeros(height(t),1);
t.meets_engagement_requirement(ismember(t.engagement_classification,engagement_classifications_to_use),:) = 1;
t.meets_time_since_reward_zone_entry_requirement = zeros(height(t),1);
t.meets_time_since_reward_zone_entry_requirement(t.time_since_drink_onset < upper_limit_time_since_reward_zone_entry & ...
    t.time_since_drink_onset > lower_limit_time_since_reward_zone_entry) = 1;
t.meets_in_end_zone_requirement = zeros(height(t),1);
t.meets_in_end_zone_requirment(t.in_reward_zone == 1) = 1;
t.distance_from_end(t.reward_zone == 1) = t.distance_from_end_1(t.reward_zone == 1);
t.distance_from_end(t.reward_zone == 2) = t.distance_from_end_2(t.reward_zone == 2);
t.meets_tight_in_end_zone_requirement = zeros(height(t),1);
t.meets_tight_in_end_zone_requirement(t.distance_from_end < distance_from_end_tight_threshold) = 1;


if sum(ismember(session_flags,12)) >= 1 % if sleep session, ignore position requirements
    t.meets_engagement_requirement = ones(height(t),1);
    t.meets_in_end_zone_requirment = ones(height(t),1);
    t.meets_tight_in_end_zone_requirement =  ones(height(t),1);
    t.meets_time_since_reward_zone_entry_requirement = ones(height(t),1);
    t.in_reward_zone = ones(height(t),1);
    t.laser_state_matches_intended_laser_state_for_stopping_period = ones(height(t),1);
end


if only_consider_replays_with_intended_laser_state == 0
    t.laser_state_matches_intended_laser_state_for_stopping_period = ones(height(t),1);
end

t(ind_session == 0,:) = [];
% Remove events that occurred outside the segment of interest from
% decoder_events.

if use_duration_og == 1
    t.event_duration = t.duration_og;
else
    t.event_duration = t.duration;
end

t.posterior_difference = abs(t.(posterior_difference_to_use{1}) - t.(posterior_difference_to_use{2}));

% good_candidate_event = all events that meet position/timing
% criterion (but nothing about weighted correlation, posterior
% etc). replay = meets all criterion; good_replay = meets all
% criterion and is long
t.good_candidate_event = ...
    t.replay_spikeDensity_power >= sde_amplitude_criterion & ...
    t.replay_ripple_power >= ripple_power_criterion & ...
    t.event_duration > duration_thr  & ...
    t.time_in_session < time_into_session_thr & ...
    t.in_reward_zone >= end_zone_thr & ...
    t.meets_tight_in_end_zone_requirement == 1 & ...
    t.meets_engagement_requirement == 1 & ...
    t.meets_time_since_reward_zone_entry_requirement == 1 & ...
    t.laser_state_matches_intended_laser_state_for_stopping_period == 1 & ...
    t.num_of_cells_participating >= num_cells_participating_thr & ...
    t.fraction_of_cells_participating >= fraction_cells_participating_thr;

t.replay = ...
    t.good_candidate_event == 1 & ...
    abs(t.weighted_r)>weighted_r_thr & ...
    t.posterior_difference>posterior_in_map_thr & ...
    t.max_jump_distance < max_jump_distance_thr & ...
        t.range_normalized > coverage_thr_percent_of_track & ...
            t.range > coverage_thr;
%     & t.range > coverage_thr & t.range_normalized >= coverage_thr_percent_of_track;


t.good_candidate_event_laser_on = t.good_candidate_event & t.laser_state == 1;
t.good_candidate_event_laser_off = t.good_candidate_event & t.laser_state == 0;
t.replay_laser_on = t.replay == 1 & t.laser_state == 1;
t.replay_laser_off = t.replay == 1 & t.laser_state == 0;
t.forward_replay = t.replay == 1 & t.direction == 1;
t.reverse_replay = t.replay ==1 & t.direction == 2;
t.forward_replay_laser_on = t.replay == 1 & t.direction == 1 & t.laser_state == 1;
t.forward_replay_laser_off = t.replay== 1 & t.direction == 1 & t.laser_state == 0;
t.reverse_replay_laser_on = t.replay == 1 & t.direction == 2 & t.laser_state == 1;
t.reverse_replay_laser_off = t.replay == 1 & t.direction == 2 & t.laser_state == 0;
t.congruent_replay = t.replay == 1 & t.congruent_with_rat_location == 1;
t.incongruent_replay = t.replay == 1 & t.incongruent_with_rat_location == 1;
t.forward_congruent_replay = t.forward_replay==1 & t.congruent_with_rat_location==1;
t.reverse_congruent_replay = t.reverse_replay==1 & t.congruent_with_rat_location==1;
t.forward_incongruent_replay = t.forward_replay==1 & t.incongruent_with_rat_location==1;
t.reverse_incongruent_replay = t.reverse_replay==1 & t.incongruent_with_rat_location==1;

session_str = {pwd};
t.session_str(:) = repmat(session_str,[height(t),1]);