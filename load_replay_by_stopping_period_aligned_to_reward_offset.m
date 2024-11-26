clearvars -except dayFiles day directory rat windows hand_clustered_only

% ripple, spike, or laser, or filtered
plot_fig = 0;
min_time_into_stopping_period = 0;
max_time_into_stopping_period = inf;
replay_weighted_r_thr = 0.6;
replay_coverage_thr = 0.2;
replay_posterior_diff_thr = 0.33;

sessionNum=1;
posSampRate = 20; % Hz

load Experiment_Information
load Analysis_Information
load Position_Data
load Behavior_Data
load laser_state

load decoder_candidateEvents
load zscored_sd
load zscored_ripple_power

Position_Data_sub = compute_dataTemporalConcatenation(Position_Data,Experiment_Information.Segments(1).Times);
laser_state_sub = compute_dataTemporalConcatenation(laser_state,Experiment_Information.Segments(1).Times);
laser_state_sub = compute_dataInterpolation(laser_state_sub,Position_Data_sub(:,1),[]);

t_events = struct();
event_hists=struct();

% Remove any stopping periods that have been nan'ed out for some reason
nan_inds = find(isnan(Reward_Epoch_Time_Boundaries_drinking(:,1)) | isnan(Reward_Epoch_Time_Boundaries_drinking(:,2)));
Reward_Epoch_Time_Boundaries_drinking(nan_inds,:) = [];
Reward_Epoch_Time_Boundaries_speed_thresholded(nan_inds,:) = [];
%% all sdes:
events = 'spike_filtered';
sde_thr = 3;
ripple_thr = -inf;
use_duration_og = 1;
posterior_diff_thr = 0;
coverage_thr = -inf;
weighted_r_thr = -inf;

t_events.('sde') = load_replays_from_individual_session(events,use_duration_og,coverage_thr,...
    weighted_r_thr,posterior_diff_thr,sde_thr,ripple_thr,min_time_into_stopping_period,max_time_into_stopping_period);

%% all ripples:
events = 'ripple_filtered';
sde_thr = -inf;
ripple_thr = 3;
use_duration_og = 1;
posterior_diff_thr = 0;
coverage_thr = -inf;
weighted_r_thr = -inf;

t_events.('ripples') = load_replays_from_individual_session(events,use_duration_og,coverage_thr,...
    weighted_r_thr,posterior_diff_thr,sde_thr,ripple_thr,min_time_into_stopping_period,max_time_into_stopping_period);

%% spike filtered replays:
events = 'spike_filtered';
sde_thr = 3;
ripple_thr = -inf;
use_duration_og = 0;
posterior_diff_thr = replay_posterior_diff_thr;
coverage_thr = replay_coverage_thr;
weighted_r_thr = replay_weighted_r_thr;

sub_events = load_replays_from_individual_session(events,use_duration_og,coverage_thr,...
    weighted_r_thr,posterior_diff_thr,sde_thr,ripple_thr,min_time_into_stopping_period,max_time_into_stopping_period);
t_events.('reverse_replays_spike') = sub_events(sub_events.reverse_replay==1,:);
t_events.('forward_replays_spike') = sub_events(sub_events.forward_replay==1,:);
t_events.('congruent_replays_spike') = sub_events(sub_events.congruent_replay==1,:);
t_events.('incongruent_replays_spike') = sub_events(sub_events.incongruent_replay==1,:);
t_events.('forward_congruent_replays_spike') = sub_events(sub_events.forward_congruent_replay==1,:);
t_events.('forward_incongruent_replays_spike') = sub_events(sub_events.forward_incongruent_replay==1,:);
t_events.('reverse_congruent_replays_spike') = sub_events(sub_events.reverse_congruent_replay==1,:);
t_events.('reverse_incongruent_replays_spike') = sub_events(sub_events.reverse_incongruent_replay==1,:);

%% ripple filtered replays:
events = 'ripple_filtered';
sde_thr = -inf;
ripple_thr = 3;
use_duration_og = 0;
posterior_diff_thr = replay_posterior_diff_thr;
coverage_thr = replay_coverage_thr;
weighted_r_thr = replay_weighted_r_thr;

sub_events = load_replays_from_individual_session(events,use_duration_og,coverage_thr,...
    weighted_r_thr,posterior_diff_thr,sde_thr,ripple_thr,min_time_into_stopping_period,max_time_into_stopping_period);
t_events.('reverse_replays_ripple') = sub_events(sub_events.reverse_replay==1,:);
t_events.('forward_replays_ripple') = sub_events(sub_events.forward_replay==1,:);
t_events.('congruent_replays_ripple') = sub_events(sub_events.congruent_replay==1,:);
t_events.('incongruent_replays_ripple') = sub_events(sub_events.incongruent_replay==1,:);
t_events.('forward_congruent_replays_ripple') = sub_events(sub_events.forward_congruent_replay==1,:);
t_events.('forward_incongruent_replays_ripple') = sub_events(sub_events.forward_incongruent_replay==1,:);
t_events.('reverse_congruent_replays_ripple') = sub_events(sub_events.reverse_congruent_replay==1,:);
t_events.('reverse_incongruent_replays_ripple') = sub_events(sub_events.reverse_incongruent_replay==1,:);

%%

hist_length = 10; % sec
hist_bin_width = 0.5; % sec
pre_stopping_period = 10; % sec
hist_length = hist_length + pre_stopping_period; % sec

% Reward_Epoch_Time_Boundaries_drinking contains the start and stop times of each
% drinking period
drinking_period_boundaries = compute_dataTemporalConcatenation(Reward_Epoch_Time_Boundaries_drinking,Experiment_Information.Segments(sessionNum).Times);

% Reward epoch boundaries has the start and stop times of the entire
% stopping period
stopping_period_boundaries = compute_dataTemporalConcatenation(Reward_Epoch_Time_Boundaries_speed_thresholded,Experiment_Information.Segments(sessionNum).Times);



% For each visit to a reward zone, count the number of forward/reverse
% replays.
drinking_period_laser_state = nan(size(drinking_period_boundaries,1),1);
drinking_period_duration = nan(size(drinking_period_boundaries,1),1);

event_counts = nan(length(t_events),size(drinking_period_boundaries,1));
event_types = fieldnames(t_events);

reward_locations = [0 Experiment_Information.maze_size];

if plot_fig == 1
    fig1 = figure;
    plot(Position_Data_sub(laser_state_sub(:,2)==0,1)./spikeSampRate,Position_Data_sub(laser_state_sub(:,2)==0,2),'.k')
    hold on
    plot(Position_Data_sub(laser_state_sub(:,2)==1,1)./spikeSampRate,Position_Data_sub(laser_state_sub(:,2)==1,2),'.r')
    hold on

    fig2 = figure;
    plot(Position_Data_sub(:,1)./spikeSampRate,Position_Data_sub(:,2),'.','color',[0.5 0.5 0.5]);
    hold on

    fig3 = figure;
    plot(Position_Data_sub(:,1)./spikeSampRate,Position_Data_sub(:,2),'.','color',[0.5 0.5 0.5]);
    hold on
end

t_event_summary = table();
t_event_summary.event_type(:) = fieldnames(t_events);


for i = 1:size(drinking_period_boundaries,1)

    drink_start_ind = find(Position_Data_sub(:,1) == drinking_period_boundaries(i,1)) ;
    drink_stop_ind =  find(Position_Data_sub(:,1) == drinking_period_boundaries(i,2));

    stopping_period_start_ind = find(Position_Data_sub(:,1) == stopping_period_boundaries(i,1));
    stopping_period_stop_ind = find(Position_Data_sub(:,1) == stopping_period_boundaries(i,2));

    % Duration of the stopping period should only include a) stationary
    % periods and b) periods where the laser state matches the intended
    % laser state for the stopping period.
    % NEW: stopping period should only include the times when the
    % rat was very close to the reward (within 15 cm);
    drinking_period = (Position_Data_sub(:,1)>= drinking_period_boundaries(i,1) & Position_Data_sub(:,1)<= drinking_period_boundaries(i,2));
    stopping_period = (Position_Data_sub(:,1)>= stopping_period_boundaries(i,1) & Position_Data_sub(:,1)<= stopping_period_boundaries(i,2));

    if sum(drinking_period)/posSampRate > max_time_into_stopping_period  + 1
        keyboard
    end

    stationary = Position_Data_sub(:,5) < speedThr;

    position_distance_from_reward_one = abs(Position_Data_sub(:,2) - reward_locations(1));
    position_distance_from_reward_two = abs(Position_Data_sub(:,2) - reward_locations(2));
    position_distance_from_nearest_reward = min([position_distance_from_reward_one position_distance_from_reward_two],[],2);

    % determine the intended laser state of the stopping period.
    % Each stopping period can only be either laser off, or laser on.
    laser_state_sub_sub = laser_state_sub(drinking_period==1 & stationary ==1,:);

    %10/27/22: changing so that if the laser was on at all, the
    %intended laser state is ON. %10/31/22: this doesn't work well
    %for sessions when laser was on during the run.
    session_flags = Experiment_Information.Segments(sessionNum).Flags;
    if ismember(1,session_flags)
        if sum(laser_state_sub_sub(:,2)) > 0
            intended_laser_state = 1;
        else
            intended_laser_state = 0;
        end
    elseif ismember(2,session_flags)
        intended_laser_state = mode(laser_state_sub_sub(:,2));
    else
        intended_laser_state = 0;
    end

    drinking_period_laser_state(i) = intended_laser_state;

    % find the indices that meet all the criterion for this stopping
    % period:


    drinking_period_inds_correct_laser_state = find(drinking_period==1 & laser_state_sub(:,2) == intended_laser_state);
    stopping_period_inds_correct_laser_state = find(stopping_period==1 & laser_state_sub(:,2) == intended_laser_state);



    % it's okay if the laser got turned off slightly before the end
    % of the stopping period or after the beginning, but the
    % portion that is being examined needs to be one continuous
    % segment. If there were 2 segments for some reason, take the
    % last one.

    if any(diff(drinking_period_inds_correct_laser_state) ~= 1)
        last_switch_to_correct_laser_state = find(diff(drinking_period_inds_correct_laser_state) > 1,1,'last');
        drinking_period_inds_correct_laser_state = drinking_period_inds_correct_laser_state(last_switch_to_correct_laser_state+1:end);
    end

    % The start and stop times that will be used for this segment:
    t_drink_start = Position_Data_sub(drinking_period_inds_correct_laser_state(1),1);
    t_drink_stop = Position_Data_sub(drinking_period_inds_correct_laser_state(end),1);
    
    % The start and stop times that will be used for this segment:

        

    if plot_fig == 1
        figure(fig1)
        if intended_laser_state == 0
            plot(Position_Data_sub(drink_start_ind,1)./spikeSampRate,Position_Data_sub(drink_start_ind,2),'ob'); hold on;
            plot(Position_Data_sub(drink_stop_ind,1)./spikeSampRate,Position_Data_sub(drink_stop_ind,2),'ob'); hold on;
           
        elseif intended_laser_state == 1
            plot(Position_Data_sub(drink_start_ind,1)./spikeSampRate,Position_Data_sub(drink_start_ind,2),'ob'); hold on;
            plot(Position_Data_sub(drink_stop_ind,1)./spikeSampRate,Position_Data_sub(drink_stop_ind,2),'ob'); hold on;
            
        end
    end

    drinking_period_duration(i) = length(drinking_period_inds_correct_laser_state)*(1/posSampRate);
    if drinking_period_duration(i)==0
        continue
    end

    for j = 1:length(event_types)

        % find the events that occured during this stopping period
        replay_times = t_events.(event_types{j}).timePoints_og(:,1);
        replay_occured_during_stopping_period=ismember(replay_times,compute_dataTemporalConcatenation(replay_times,[t_drink_start t_drink_stop]),'rows');

        event_inds = find(t_events.(event_types{j}).laser_state_matches_intended_laser_state_for_stopping_period & replay_occured_during_stopping_period);
        event_counts(j,i) = length(event_inds);

        % The time that the rat finished drinking will be t=0. we will look
        % +/- 5 seconds out from this time.

        hist_stop_time = Reward_Epoch_Time_Boundaries_drinking(i,2);
        hist_start_time = hist_stop_time - 14*spikeSampRate;
%         hist_start_time = hist_mid_time - 5*spikeSampRate;
%         hist_stop_time = hist_mid_time + 5*spikeSampRate;


        t_bin_edges = hist_start_time:hist_bin_width*spikeSampRate:hist_stop_time;
        t_bin_centers = (t_bin_edges(1:end-1)+t_bin_edges(2:end))./2;


        event_hists.(event_types{j})(i,:) =  histcounts(replay_times(event_inds,1),t_bin_edges);
        event_hists.(event_types{j})(i,t_bin_centers < stopping_period_boundaries(i,1) | t_bin_centers > stopping_period_boundaries(i,2)) = nan;

        % Also nan out any time periods that had the wrong laser state
        laser_state_sub_sub = compute_dataInterpolation(laser_state_sub,t_bin_centers',[]);
        laser_state_sub_sub = laser_state_sub_sub(:,2);
        event_hists.(event_types{j})(i,laser_state_sub_sub ~= intended_laser_state) = nan;
    end
end
%%
t2 = table;
t2.pass_number = (1:size(drinking_period_boundaries,1))';
t2.laser_state = drinking_period_laser_state;
laser_off_pass_count = 1;
laser_on_pass_count = 1;

for i = 1:height(t2)
    if t2.laser_state(i) == 0
        t2.laser_state_pass_count(i) = laser_off_pass_count;
        laser_off_pass_count = laser_off_pass_count + 1;
    elseif t2.laser_state(i) == 1
        t2.laser_state_pass_count(i) = laser_on_pass_count;
        laser_on_pass_count = laser_on_pass_count + 1;
    end
end

t2.duration = drinking_period_duration;

for event_type = 1:size(event_counts,1)
t2.(event_types{event_type})(:) = event_counts(event_type,:);
end
t2(t2.duration == 0,:) = [];

save('replay_by_stopping_period_aligned_to_drink_offset','t2','event_hists');



