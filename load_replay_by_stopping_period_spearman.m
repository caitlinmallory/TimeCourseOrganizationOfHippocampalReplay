clearvars -except dayFiles day directory rat windows hand_clustered_only

% ripple, spike, or laser, or filtered
plot_fig = 0;
min_time_into_stopping_period = 0;
max_time_into_stopping_period = inf;
distance_from_end_tight_threshold = 30;

sessionNum=1;
posSampRate = 20; % Hz

load Experiment_Information
load Analysis_Information
load Position_Data
load Behavior_Data
load laser_state

load decoder_candidateEvents_noLFP


Position_Data_sub = compute_dataTemporalConcatenation(Position_Data,Experiment_Information.Segments(1).Times);
laser_state_sub = compute_dataTemporalConcatenation(laser_state,Experiment_Information.Segments(1).Times);
laser_state_sub = compute_dataInterpolation(laser_state_sub,Position_Data_sub(:,1),[]);

Reward_Epoch_Time_Boundaries_speed_thresholded(isnan(sum(Reward_Epoch_Time_Boundaries_speed_thresholded,2)),:) = [];
% If there are any reward epoch time boundaries with NaN's remaining,
% remove these rows.

t_events = struct();
event_hists=struct();

%% spike filtered replays:
events = 'spike_filtered';

sub_events = load_replays_from_individual_session_spearman(events,min_time_into_stopping_period,max_time_into_stopping_period);

% t_events.('reverse_replays_spike') = sub_events(sub_events.reverse_replay==1,:);
% t_events.('forward_replays_spike') = sub_events(sub_events.forward_replay==1,:);
% t_events.('congruent_replays_spike') = sub_events(sub_events.congruent_replay==1,:);
% t_events.('incongruent_replays_spike') = sub_events(sub_events.incongruent_replay==1,:);
% t_events.('forward_incongruent_replays_spike') = sub_events(sub_events.forward_incongruent_replay==1,:);
% t_events.('reverse_incongruent_replays_spike') = sub_events(sub_events.reverse_incongruent_replay==1,:);

t_events.('forward_congruent_replays_spike') = sub_events(sub_events.forward_congruent_replay==1,:);
t_events.('reverse_congruent_replays_spike') = sub_events(sub_events.reverse_congruent_replay==1,:);

num_event_types = length(fieldnames(t_events));
%%

hist_length = 14; % sec % We will smooth the end
hist_bin_width = 0.5; % sec
pre_stopping_period = 0; % sec
hist_length = hist_length + pre_stopping_period; % sec

% Reward_Epoch_Time_Boundaries contains the start and stop times of each
% visit to a reward zone.
stopping_period_boundaries = compute_dataTemporalConcatenation(Reward_Epoch_Time_Boundaries_speed_thresholded,Experiment_Information.Segments(sessionNum).Times);
%stopping_period_boundaries = compute_dataTemporalConcatenation(Reward_Epoch_Time_Boundaries_drinking,Experiment_Information.Segments(sessionNum).Times);

% For each visit to a reward zone, count the number of forward/reverse
% replays.
stopping_period_laser_state = nan(size(stopping_period_boundaries,1),1);
stopping_period_duration_stationary = nan(size(stopping_period_boundaries,1),1);
stopping_period_duration = (Reward_Epoch_Time_Boundaries_speed_thresholded(:,2)-Reward_Epoch_Time_Boundaries_speed_thresholded(:,1))./spikeSampRate;
stopping_period_reward_zone = nan(size(stopping_period_boundaries,1),1);

event_counts = nan(num_event_types,size(stopping_period_boundaries,1));
first_event_times = nan(num_event_types,size(stopping_period_boundaries,1));
last_event_times = nan(num_event_types,size(stopping_period_boundaries,1));

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

% If the stopping period has been limited to the first X-seconds,
% change stopping period exits to be whichever is first, the actual end
% time of the stopping period, or the start time of the stopping period +
% x-seconds.
stopping_period_boundaries(:,2) =  min(stopping_period_boundaries(:,2),stopping_period_boundaries(:,1) + max_time_into_stopping_period*spikeSampRate);


for i = 1:size(stopping_period_boundaries,1)

    start_ind = find(Position_Data_sub(:,1) == stopping_period_boundaries(i,1)) ;
    stop_ind =  find(Position_Data_sub(:,1) == stopping_period_boundaries(i,2));

    % Which end did this stopping period occur on?
    stopping_period_reward_zone(i) = mode(Pass_Info(Pass_Info(:,1)>stopping_period_boundaries(i,1) & Pass_Info(:,1)<stopping_period_boundaries(i,2),3));


    % Duration of the stopping period should only include a) stationary
    % periods and b) periods where the laser state matches the intended
    % laser state for the stopping period.
    % NEW: stopping period should only include the times when the
    % rat was very close to the reward (within 15 cm);
    stopping_period = (Position_Data_sub(:,1)>= stopping_period_boundaries(i,1) & Position_Data_sub(:,1)<= stopping_period_boundaries(i,2));


    % Limit the stopping period to the first X-seconds (set by
    % upper_lmit_time_since_reward_zone_entry). This will exclude
    % portions of the stopping period exceeding this value- meant
    % to eliminate times when the animal spent a very long time at
    % one end.
    stopping_period (Position_Data_sub(:,1) > stopping_period_boundaries(i,1) +  spikeSampRate*max_time_into_stopping_period) = 0;

    if sum(stopping_period)/posSampRate > max_time_into_stopping_period  + 1
        keyboard
    end

    stationary = Position_Data_sub(:,5) < speedThr;

    position_distance_from_reward_one = abs(Position_Data_sub(:,2) - reward_locations(1));
    position_distance_from_reward_two = abs(Position_Data_sub(:,2) - reward_locations(2));
    position_distance_from_nearest_reward = min([position_distance_from_reward_one position_distance_from_reward_two],[],2);
    in_endzone_tight = position_distance_from_nearest_reward <= distance_from_end_tight_threshold;

    % determine the intended laser state of the stopping period.
    % Each stopping period can only be either laser off, or laser on.
    laser_state_sub_sub = laser_state_sub(stopping_period==1 & stationary ==1,:);

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

    stopping_period_laser_state(i) = intended_laser_state;

    % find the indices that meet all the criterion for this stopping
    % period:

    stopping_period_inds = find(stopping_period==1 & stationary==1 & in_endzone_tight==1 & laser_state_sub(:,2) == intended_laser_state);

    if plot_fig == 1
        figure(fig1)
        if intended_laser_state == 0
            plot(Position_Data_sub(start_ind,1)./spikeSampRate,Position_Data_sub(start_ind,2),'ok','MarkerFaceColor','k'); hold on;
            plot(Position_Data_sub(stop_ind,1)./spikeSampRate,Position_Data_sub(stop_ind,2),'ok','MarkerFaceColor','k'); hold on;
            plot(Position_Data_sub(stopping_period_inds,1)./spikeSampRate,Position_Data_sub(stopping_period_inds,2),'ok'); hold on;
        elseif intended_laser_state == 1
            plot(Position_Data_sub(start_ind,1)./spikeSampRate,Position_Data_sub(start_ind,2),'or','MarkerFaceColor','r'); hold on;
            plot(Position_Data_sub(stop_ind,1)./spikeSampRate,Position_Data_sub(stop_ind,2),'or','MarkerFaceColor','r'); hold on;
            plot(Position_Data_sub(stopping_period_inds,1)./spikeSampRate,Position_Data_sub(stopping_period_inds,2),'or'); hold on;
        end

        figure(fig2)
        if intended_laser_state == 0
            plot(Position_Data_sub(start_ind,1)./spikeSampRate,Position_Data_sub(start_ind,2),'ok','MarkerFaceColor','k'); hold on;
            plot(Position_Data_sub(stop_ind,1)./spikeSampRate,Position_Data_sub(stop_ind,2),'ok','MarkerFaceColor','k'); hold on;
            plot(Position_Data_sub(stopping_period_inds,1)./spikeSampRate,Position_Data_sub(stopping_period_inds,2),'.k'); hold on;
        elseif intended_laser_state == 1
            plot(Position_Data_sub(start_ind,1)./spikeSampRate,Position_Data_sub(start_ind,2),'or','MarkerFaceColor','r'); hold on;
            plot(Position_Data_sub(stop_ind,1)./spikeSampRate,Position_Data_sub(stop_ind,2),'or','MarkerFaceColor','r'); hold on;
            plot(Position_Data_sub(stopping_period_inds,1)./spikeSampRate,Position_Data_sub(stopping_period_inds,2),'.r'); hold on;
        end
    end


    stopping_period_duration_stationary(i) = length(stopping_period_inds)*(1/posSampRate);

    if stopping_period_duration_stationary(i)==0
        continue
    end


    for j = 1:length(event_types)

        % find the events that occured during this stopping period
        replay_times = t_events.(event_types{j}).timePoints_og(:,1);
        replay_occured_during_stopping_period=ismember(replay_times,compute_dataTemporalConcatenation(replay_times,[stopping_period_boundaries(i,1) stopping_period_boundaries(i,2)]),'rows');

        event_inds = find(t_events.(event_types{j}).laser_state_matches_intended_laser_state_for_stopping_period & replay_occured_during_stopping_period);
        event_counts(j,i) = length(event_inds);
        if ~isempty(event_inds)
            first_event_times(j,i) = t_events.(event_types{j}){event_inds(1),'time_since_reward_zone_entry'};
            last_event_times(j,i) = t_events.(event_types{j}){event_inds(end),'time_till_reward_zone_exit'};
        end
        % to plot replay rate over time, we need to be the stopping period.
        % time needs to be continous. Use the first and last index that met all
        % the criterion above to compute the start and stop times for the
        % stopping period.

        t_start = Position_Data_sub(stopping_period_inds(1),1) - pre_stopping_period*spikeSampRate;
        t_stop = Position_Data_sub(stopping_period_inds(end),1);

        t_bin_edges = t_start:hist_bin_width*spikeSampRate:(t_start+hist_length*spikeSampRate);
        t_bin_centers = (t_bin_edges(1:end-1)+t_bin_edges(2:end))./2;

        event_hists.(event_types{j})(i,:) =  histcounts(replay_times(event_inds,1),t_bin_edges);
        event_hists.(event_types{j})(i,t_bin_centers > t_stop) = nan;
    end
end
%%
t2 = table;
t2.pass_number = (1:size(stopping_period_boundaries,1))';
t2.laser_state = stopping_period_laser_state;
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

t2.duration = stopping_period_duration;
t2.duration_stationary = stopping_period_duration_stationary;
t2.reward_zone = stopping_period_reward_zone;
first_event_time_table = table();
last_event_time_table = table();
for event_type = 1:size(event_counts,1)
    t2.(event_types{event_type})(:) = event_counts(event_type,:);
    first_event_time_table.(event_types{event_type})(:) = first_event_times(event_type,:)';
    last_event_time_table.(event_types{event_type})(:) = last_event_times(event_type,:)';
end
rows_to_remove = find(t2.duration==0);


t2(rows_to_remove,:) = [];
first_event_time_table(rows_to_remove,:) = [];
last_event_time_table(rows_to_remove,:) = [];

save('replay_by_stopping_period_pearsons','t2','first_event_time_table','last_event_time_table','event_hists');



