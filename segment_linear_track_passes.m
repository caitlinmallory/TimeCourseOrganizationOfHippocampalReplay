clearvars -except dayFiles day directory rat windows hand_clustered_only
load Experiment_Information
load Position_Data
load laser_state

plot_figure = 1;

speedThr = 5;
drinking_speedThr = 5;
end_zone_size = 30;
min_drink_duration_thr = 2;
deltThr = 0;
reward_zone_tolerance = 5;
Experiment_Information.spikeSampRate = 30000;
% Main output being used is Reward_Epoch_Time_Boundaries_speed_thresholded
% Entry times are the first times on each pass when the rat is within reward_zone_tolerance of the
% reward location and moving below speedThr. Exit times are
% the last time the rat is within end_zone_size and moving below speedThr.


run_segments = [];
for i = 1:length(Experiment_Information.Segments)
    if sum(ismember(Experiment_Information.Segments(i).Flags,11)) > 0
        run_segments = [run_segments; i];
    end
end


reward_locations = [0 Experiment_Information.maze_size];

Pass_Info = [];
Pass_Transitions = [];
Run_Epoch_Time_Boundaries = [];
Reward_Epoch_Time_Boundaries = [];
Reward_Epoch_Time_Boundaries_matched = [];
Reward_Epoch_Time_Boundaries_speed_thresholded = [];

for i = 1:length(run_segments)


    session_times = Experiment_Information.Segments(run_segments(i)).Times;

    session_start_time = session_times(1);

    Position_Data_sub = compute_dataTemporalConcatenation(Position_Data,session_times);
    laser_state_sub = compute_dataTemporalConcatenation(laser_state,session_times);
    laser_state_sub = compute_dataInterpolation(laser_state_sub,Position_Data_sub(:,1),[]);
    
    end_one = find(Position_Data_sub(:,2) < reward_locations(1) + end_zone_size);
    end_two = find(Position_Data_sub(:,2) > reward_locations(2) - end_zone_size);


    % Segment into laps and find the average time spent per lap at each end,
    % and each end stationary
    lap_info = zeros(size(Position_Data_sub(:,6)));
    lap_info(end_one) = 1;
    lap_info(end_two) = 2;

    diff_lap_info = diff(lap_info);
    diff_lap_info = [0; diff_lap_info];

    end_transitions = find(abs(diff_lap_info) > 0);
    diff_end_transitions = diff(diff_lap_info(end_transitions));
    diff_end_transitions = [0; diff_end_transitions];

    % exits from left side are -1
    % entries to right side are 2
    % exits from right side are -2
    % entries to left side are 1

    % difference in transitions should be 3 -4 3 -2 ...
    % an exit from the right followed by an entry to the right is 2 - (-2) = 4
    % an exit from the left followed by a re-entry to the left is 1 - (-1) = 2
    % these are "bad transitions", don't consider them to be separate laps.

    bad_transitions = find(diff_end_transitions == 2 | diff_end_transitions == 4);

    lap_info_clean = lap_info;
    for n = 1:length(bad_transitions)
        lap_info_clean(end_transitions(bad_transitions(n)-1):end_transitions(bad_transitions(n))) = lap_info_clean(end_transitions(bad_transitions(n)-1)-1);
    end

    diff_lap_info_clean = diff(lap_info_clean);
    diff_lap_info_clean = [0; diff_lap_info_clean];

    end_transitions = find(abs(diff_lap_info_clean) > 0);
    diff_end_transitions = diff(diff_lap_info_clean(end_transitions));
    diff_end_transitions = [0; diff_end_transitions];

    % left_exits = find(diff_lap_info_clean == -1);
    % right_exits = find(diff_lap_info_clean == -2);
    exits = find(diff_end_transitions < 0);

    % If the animal started in and end, add in the first exit that wasn't
    % identified in the manner above.
    if lap_info_clean(1) > 0
        exits = [1; exits];
    end

    for n = 1:length(exits)
        if exits(n) == length(end_transitions)
            if lap_info_clean(end_transitions(exits(n))-1) == 1
                lap_info_clean(end_transitions(exits(n)):end) = -2;
            elseif lap_info_clean(end_transitions(exits(n))-1) == 2
                lap_info_clean(end_transitions(exits(n)):end) = -1;
            end
        else
            if lap_info_clean(end_transitions(exits(n))-1) == 1
                lap_info_clean(end_transitions(exits(n)):end_transitions(exits(n)+1)-1) = -2;
            elseif lap_info_clean(end_transitions(exits(n))-1) == 2
                lap_info_clean(end_transitions(exits(n)):end_transitions(exits(n)+1)-1) = -1;
            end
        end
    end

    %find lap data if the rat started in the middle of the track heading left:
    if lap_info_clean(end_transitions(1)) == 1
        lap_info_clean(1:end_transitions(1)-1) = -1;
    end
    %find lap data if the rat started in the middle of the track heading right:
    if lap_info_clean(end_transitions(1)) == 2
        lap_info_clean(1:end_transitions(1)-1) = -2;
    end

    % Assign a pass number to each pass
    Pass_Info_sub = zeros(size(Position_Data_sub(:,6)));
    pass = 1;
    for n = 1:length(exits)
        if n == length(exits)
            Pass_Info_sub(end_transitions(exits(n)):end) = pass;
            pass = pass;
        else
            Pass_Info_sub(end_transitions(exits(n)):end_transitions(exits(n+1))-1) = pass;
            pass = pass + 1;
        end
    end
    % if the animal stopped in the middle of the track, include that last
    % segment in the last pass
    if lap_info_clean(end) < 0
        Pass_Info_sub(end_transitions(end):length(Pass_Info_sub)) = pass;
    end

    Pass_Info_sub = Pass_Info_sub + 1;


    if plot_figure == 1

        figure()
        p1 = plot(Position_Data_sub(:,1)./30000 - session_start_time/30000,Position_Data_sub(:,2),'LineWidth',2);
        xlabel('time (s)')
        ylabel('X-Position (cm)')
        hold on;
        p2 = plot(Position_Data_sub(end_one,1)./30000 - session_start_time/30000,Position_Data_sub(end_one,2),'.k','MarkerSize',10);
        hold on;
        p3 = plot(Position_Data_sub(end_two,1)./30000 - session_start_time/30000,Position_Data_sub(end_two,2),'.r','MarkerSize',10);
        legend([p2 p3],{'Left', 'Right'},'Location','northoutside')

        figure();
        plot(Position_Data_sub(lap_info_clean == 1,1)./30000 - session_start_time/30000,Position_Data_sub(lap_info_clean == 1,2),'.k','MarkerSize',10);
        hold on;
        plot(Position_Data_sub(lap_info_clean == 2,1)./30000 - session_start_time/30000,Position_Data_sub(lap_info_clean == 2,2),'.r','MarkerSize',10);
        hold on
        plot(Position_Data_sub(lap_info_clean == -1,1)./30000 - session_start_time/30000,Position_Data_sub(lap_info_clean == -1,2),'.b','MarkerSize',10);
        hold on
        plot(Position_Data_sub(lap_info_clean == -2,1)./30000 - session_start_time/30000,Position_Data_sub(lap_info_clean == -2,2),'.g','MarkerSize',10);
        hold on
        xlabel('time (s)')
        ylabel('X-Position (cm)')
        legend({'Left-end', 'Right-end', 'Moving left', 'Moving right'},'Location','northoutside')


        figure();
        plot(Position_Data_sub(:,1)./30000 - session_start_time/30000,Position_Data_sub(:,2),'color',[0.5 0.5 0.5]);
        hold on
        for p = 1:length(unique(Pass_Info_sub(~isnan(Pass_Info_sub))))
            plot(Position_Data_sub(Pass_Info_sub == p,1)./30000 - session_start_time/30000,Position_Data_sub(Pass_Info_sub == p,2));
            hold on
        end
    end


    direction_of_transition = nan(length(end_transitions),1);
    direction_of_transition(Position_Data_sub(end_transitions,2)-Position_Data_sub(end_transitions-1,2)>0) = 1;
    direction_of_transition(Position_Data_sub(end_transitions,2)-Position_Data_sub(end_transitions-1,2)<0) = -1;
    location_of_transition = ones(length(end_transitions),1);

    location_of_transition(abs(Position_Data_sub(end_transitions,2)-end_zone_size) > abs(Position_Data_sub(end_transitions,2)-(Experiment_Information.maze_size-end_zone_size))) = 2;


    % transitions has to columns:
    % col 1 = the timestamp of the transition into or out of an end zone.
    % col 2 = the key corresponding to the exit/entry type

    % Entry/Exsit key:
    % Entry to end 1 (left) = 1
    % Exit from end 1 (left) = -1
    % Entry to end 2 (right) = 2
    % Exit from end 2 (right) = -2

    Pass_Transitions_sub = [Position_Data_sub(end_transitions,1) nan(length(end_transitions),1)];
    Pass_Transitions_sub(direction_of_transition == 1 & location_of_transition == 2,2) = 2;
    Pass_Transitions_sub(direction_of_transition == 1 & location_of_transition == 1,2) = -1;
    Pass_Transitions_sub(direction_of_transition == -1 & location_of_transition == 2,2) = -2;
    Pass_Transitions_sub(direction_of_transition == -1 & location_of_transition == 1,2) = 1;

    if plot_figure == 1
        figure()
        plot(Position_Data_sub(laser_state_sub(:,2)==0,1)./30000 - session_start_time/30000,Position_Data_sub(laser_state_sub(:,2)==0,2),'.k');
        hold on
        plot(Position_Data_sub(laser_state_sub(:,2)==1,1)./30000 - session_start_time/30000,Position_Data_sub(laser_state_sub(:,2)==1,2),'.r');
        for pass = 1:length(Pass_Transitions_sub)
            ind = find(Position_Data_sub(:,1) == Pass_Transitions_sub(pass,1));
            plot(Position_Data_sub(ind,1)./30000 - session_start_time/30000,Position_Data_sub(ind,2),'.k'); hold on;
        end
    end

    top_boundary_one = reward_locations(1) + end_zone_size;
    bottom_boundary_two = reward_locations(2) - end_zone_size;
    location_rat = ones(length(Position_Data_sub),1);
    location_rat(abs(Position_Data_sub(:,2)-top_boundary_one) > abs(Position_Data_sub(:,2)-bottom_boundary_two)) = 2;

    Pass_Info_sub = [Position_Data_sub(:,1) Pass_Info_sub location_rat];
    % Also add which end the rat is closer to at end given time (1 = left, 2 =
    % right)

    % Find the start and stop times of each reward zone entry/exit.
    % the first transition is the rat leaving the end zone; But we want to look
    % at replays that occured before that.

    Reward_Epoch_Time_Boundaries_sub = nan(length(unique(Pass_Info_sub(:,2))),2);
    Reward_Epoch_Time_Boundaries_sub_speed_thresholded = nan(length(unique(Pass_Info_sub(:,2))),2);

    for n = 1:length(unique(Pass_Info_sub(:,2)))
        correct_pass = Pass_Info_sub(:,2)== n;
        % in_end_zone = Position_Data_sub(:,2) < end_zone_size | Position_Data_sub(:,2) > (Experiment_Information.maze_size-end_zone_size);
        in_end_zone = Position_Data_sub(:,2) < top_boundary_one | Position_Data_sub(:,2) > bottom_boundary_two;
        stopping_period_inds = find(correct_pass == 1 & in_end_zone == 1);

        if length(stopping_period_inds) > 2
            % if ~isempty(stopping_period_inds)
            % if the animal moved in and out of the end zone, take the
            % first chunk where he was in end zone and stationary
            stopping_period_inds_diff = diff(stopping_period_inds);
            stopping_period_inds_diff = [stopping_period_inds_diff(1); stopping_period_inds_diff];
            stopping_period_inds_diff(stopping_period_inds_diff~=1) = nan;
            [boundaries,lengths] = compute_allSequences_NaNseparated(stopping_period_inds_diff);
            [boundaries_merged,lengths_merged] = compute_allSequences_NaNseparated_merge_stopping_periods(Position_Data_sub,stopping_period_inds,boundaries,end_zone_size,inf,1,Experiment_Information.spikeSampRate);

            % instead of taking the longest segment, I think I should take
            % the first.
            %[~,best_segment_ind] = max(lengths_merged);
            boundaries_merged = boundaries_merged(1,:);

            stopping_period_start_ind = stopping_period_inds(boundaries_merged(1));
            stopping_period_end_ind = stopping_period_inds(boundaries_merged(2));

            Reward_Epoch_Time_Boundaries_sub(n,1) = Position_Data_sub(stopping_period_start_ind,1);
            Reward_Epoch_Time_Boundaries_sub(n,2) = Position_Data_sub(stopping_period_end_ind,1);


            if plot_figure == 1
                plot(Position_Data_sub(stopping_period_start_ind,1)./30000 - session_start_time/30000,Position_Data_sub(stopping_period_start_ind,2),'ok','MarkerFaceColor','k'); hold on;
                plot(Position_Data_sub(stopping_period_end_ind,1)./30000 - session_start_time/30000,Position_Data_sub(stopping_period_end_ind,2),'ok','MarkerFaceColor','k'); hold on;
                plot(Position_Data_sub(:,1)./30000 - session_start_time/30000,Position_Data_sub(:,5),'m'); hold on;
            end
            sub_position_segment = Position_Data_sub(stopping_period_start_ind:stopping_period_end_ind,:);
            laser_state_sub_sub = laser_state_sub(stopping_period_start_ind:stopping_period_end_ind,:);
   
            % find end-zone start time: the first time the rat is in the
            % end zone and velocity < 5 cm/s
            % find end-zone stop time: the last time the rate is in the end
            % zone and velocity < 5 cm/s
            distance_from_reward = min(abs(sub_position_segment(:,2)-reward_locations(1)),abs(sub_position_segment(:,2)-reward_locations(2)));
            sub_ind_start = find(sub_position_segment(:,5)<speedThr & distance_from_reward < reward_zone_tolerance,1,'first');
            sub_ind_end = find(sub_position_segment(:,5)<speedThr,1,'last');

            if isempty(sub_ind_start)
                Reward_Epoch_Time_Boundaries_sub_speed_thresholded(n,1) = nan;
                Reward_Epoch_Time_Boundaries_sub_speed_thresholded(n,2) = nan;
            else
                Reward_Epoch_Time_Boundaries_sub_speed_thresholded(n,1) = sub_position_segment(sub_ind_start,1);
                Reward_Epoch_Time_Boundaries_sub_speed_thresholded(n,2) = sub_position_segment(sub_ind_end,1);

                if plot_figure == 1
                    plot(sub_position_segment(sub_ind_start,1)./30000 - session_start_time/30000,sub_position_segment(sub_ind_start,2),'og'); hold on;
                    plot(sub_position_segment(sub_ind_end,1)./30000 - session_start_time/30000,sub_position_segment(sub_ind_end,2),'or'); hold on;                    
                title('green=end zone entry < 5 cm/s, red=last exit moving > 5 cm/s')
                end
            end

        else
            Reward_Epoch_Time_Boundaries_sub(n,1) = nan;
            Reward_Epoch_Time_Boundaries_sub(n,2) = nan;
            Reward_Epoch_Time_Boundaries_sub_speed_thresholded(n,1) = nan;
            Reward_Epoch_Time_Boundaries_sub_speed_thresholded(n,2) = nan;
        end
    end


    % Minor fix to pass transitions: if the first reward zone transition is an
    % exit, make the start of the Position Data be the 'entry' (i.e., the rat
    % started in the end zone). If the last transition is an entry, make the
    % end of the Position data be the last 'exit' (i.e., the rat ended in the
    % reward zone).
    if Pass_Transitions_sub(1,2) < 0
        Pass_Transitions_sub = [Position_Data_sub(1,1) -1*Pass_Transitions_sub(1,2); Pass_Transitions_sub];
    end
    if Pass_Transitions_sub(end,2) > 0
        Pass_Transitions_sub = [Pass_Transitions_sub; Position_Data_sub(end,1) -1*Pass_Transitions_sub(end,2)];
    end

    bad_inds = find(isnan(Reward_Epoch_Time_Boundaries_sub(:,1)));

    Reward_Epoch_Time_Boundaries_sub(bad_inds,:) = [];
    Reward_Epoch_Time_Boundaries_sub_speed_thresholded(bad_inds,:) = [];

    Pass_Info = [Pass_Info; Pass_Info_sub];
    Pass_Transitions = [Pass_Transitions; Pass_Transitions_sub];
    Reward_Epoch_Time_Boundaries = [Reward_Epoch_Time_Boundaries; Reward_Epoch_Time_Boundaries_sub];
    Reward_Epoch_Time_Boundaries_speed_thresholded = [Reward_Epoch_Time_Boundaries_speed_thresholded; Reward_Epoch_Time_Boundaries_sub_speed_thresholded];


    % Laps should start with an end zone exit and end with and end zone entry:
    lap_starts = Pass_Transitions(Pass_Transitions(:,2)<0,1);
    lap_stops = Pass_Transitions(Pass_Transitions(:,2)>0,1);
    if lap_stops(1) < lap_starts(1)
        lap_stops(1) = [];
    end
    if lap_starts(end) > lap_stops(end)
        lap_starts(end) = [];
    end
    Run_Epoch_Time_Boundaries = [lap_starts lap_stops];


    % New: limit reward_epoch_time_boundaries to match behavior from laser
    % on and off
    % For each end zone, let's further restrict to compareable parts on the
    % maze (determined above).
    Reward_Epoch_Time_Boundaries_matched_sub = nan(size(Reward_Epoch_Time_Boundaries_sub));

    if ismember(11,Experiment_Information.Segments(i).Flags)
        laser_state_sub = compute_dataTemporalConcatenation(laser_state,[Position_Data_sub(1,1) Position_Data_sub(end,1)]);
        laser_state_sub = compute_dataInterpolation(laser_state_sub,Position_Data_sub(:,1),[]);


        if exist('LaserTimestampData.mat')==2
            load LaserTimestampData.mat
            laser_events = laser_timestamps_cat(laser_timestamps_cat(:,2)==1,1); % laser on timestamps
        else
            disp('No laser timestamp data found! Assuming this was a control session')
            laser_events = [];
            laser_timestamps_cat = [];
        end
        if ~isempty(laser_events)


            laser_timestamps_sub = compute_dataTemporalConcatenation(laser_timestamps_cat,Experiment_Information.Segments(i).Times);
            [laser_on_windows, laser_off_windows, long_laser_off_windows] = get_laser_chunk_boundaries(laser_timestamps_sub,Experiment_Information.spikeSampRate, []);

            %Hacky: but for now, if the laser on windows were short (less than
            %110 ms), this was ripple-triggered laser session, so don't analyze
            %in chunks
            if ~isempty(laser_on_windows)
                if mode((laser_on_windows(:,2)-laser_on_windows(:,1))./Experiment_Information.spikeSampRate) < 0.11
                    laser_on_windows = [];
                    laser_off_windows = [];
                    long_laser_off_windows = [];
                end
            end

        else
            laser_on_windows = [];
            laser_off_windows = [];
            long_laser_off_windows = [];
        end

        if ~isempty(laser_on_windows)
        laser_on_positions = compute_dataInterpolation(Position_Data_sub, laser_on_windows(:,1), []);
        laser_on_positions = laser_on_positions(:,2);

        end_zone_1_entry = quantile(laser_on_positions(laser_on_positions < 75),0.75);
        end_zone_2_entry = quantile(laser_on_positions(laser_on_positions > 75),.25);

        laser_off_positions = compute_dataInterpolation(Position_Data_sub, laser_off_windows(:,1), []);
        laser_off_positions = laser_off_positions(:,2);

        end_zone_1_exit= quantile(laser_off_positions(laser_off_positions < 75),0.25);
        end_zone_2_exit = quantile(laser_off_positions(laser_off_positions > 75),0.75);

        % You want to look at reward epoch time boundaries, or run laps,
        % depending on the question.


        all_epoch_inds_matched = [];
        for j = 1:length(Reward_Epoch_Time_Boundaries_sub)
            [~, start_ind ]= min(abs(Position_Data_sub(:,1) - Reward_Epoch_Time_Boundaries_sub(j,1)));
            [~, stop_ind ]= min(abs(Position_Data_sub(:,1) - Reward_Epoch_Time_Boundaries_sub(j,2)));
            epoch_inds = start_ind:stop_ind;
            Position_Data_epoch = Position_Data_sub(epoch_inds,:);
            if nanmean(Position_Data_epoch(:,2)) < 75 % hack to find end 1's
                start_ind_sub = find(Position_Data_epoch(:,2) < end_zone_1_entry,1,'first');
                Position_Data_epoch(1:start_ind_sub,:) = nan;
                stop_ind_sub = find(Position_Data_epoch(:,2) > end_zone_1_exit,1,'first');
                if isempty(stop_ind_sub) % rat never left the end zone
                    stop_ind_sub = find(~isnan(Position_Data_epoch(:,2)),1,'last');
                end
            else
                start_ind_sub = find(Position_Data_epoch(:,2) > end_zone_2_entry,1,'first');
                Position_Data_epoch(1:start_ind_sub,:) = nan;
                stop_ind_sub = find(Position_Data_epoch(:,2) < end_zone_2_exit ,1,'first');
                if isempty(stop_ind_sub) % rat never left the end zone
                    stop_ind_sub = find(~isnan(Position_Data_epoch(:,2)),1,'last');
                end
            end
            epoch_inds_matched = epoch_inds(start_ind_sub):epoch_inds(stop_ind_sub);
            all_epoch_inds_matched = [all_epoch_inds_matched; epoch_inds_matched'];
            if ~isempty(epoch_inds_matched) % occasionally the rat did not enter the end zone at all.
                Reward_Epoch_Time_Boundaries_matched_sub(j,1) = Position_Data_sub(epoch_inds_matched(1),1);
                Reward_Epoch_Time_Boundaries_matched_sub(j,2) = Position_Data_sub(epoch_inds_matched(end),1);
            end
        end

        laser_on_inds = all_epoch_inds_matched(laser_state_sub(all_epoch_inds_matched,2)==1);
        laser_off_inds =  all_epoch_inds_matched(laser_state_sub(all_epoch_inds_matched,2)==0);

%         if plot_figure
%             figure()
% 
%             plot(Position_Data_sub(:,1)./30000 - session_start_time/30000,Position_Data_sub(:,2))
%             hold on
%             plot(Position_Data_sub(laser_on_inds,1)./30000 - session_start_time/30000,Position_Data_sub(laser_on_inds,2),'.r')
%             hold on
%             plot(Position_Data_sub(laser_off_inds,1)./30000 - session_start_time/30000,Position_Data_sub(laser_off_inds,2),'.k')
%         end
    end

    Reward_Epoch_Time_Boundaries_matched = [Reward_Epoch_Time_Boundaries_matched; Reward_Epoch_Time_Boundaries_matched_sub];
    end
end

%%
% Go through and select actual drinking periods

Reward_Epoch_Time_Boundaries_drinking = nan(length(Reward_Epoch_Time_Boundaries_speed_thresholded),2);
Reward_Epoch_Time_Boundaries_endzone = nan(length(Reward_Epoch_Time_Boundaries_speed_thresholded),1);

if plot_figure == 1
figure()
plot(Position_Data(:,1)./30000 - session_start_time/30000, Position_Data(:,2));
hold on
plot(Position_Data(:,1)./30000 - session_start_time/30000, Position_Data(:,5),'color',[0.5 0.5 0.5]);
end
for i = 1:length(Reward_Epoch_Time_Boundaries_speed_thresholded)
    Position_Data_sub = compute_dataTemporalConcatenation(Position_Data,[Reward_Epoch_Time_Boundaries_speed_thresholded(i,1) Reward_Epoch_Time_Boundaries_speed_thresholded(i,2)]);

    if isempty(Position_Data_sub)
        Reward_Epoch_Time_Boundaries_drinking(i,:) = [NaN NaN];
        continue
    end

    if abs(nanmean(Position_Data_sub(:,2))-reward_locations(1)) < abs(nanmean(Position_Data_sub(:,2))-reward_locations(2))
        Reward_Epoch_Time_Boundaries_endzone(i) = 1;
    else
        Reward_Epoch_Time_Boundaries_endzone(i) = 2;
    end

    drinking = (Position_Data_sub(:,2) > reward_locations(1) - reward_zone_tolerance & ...
        Position_Data_sub(:,2) < reward_locations(1) + reward_zone_tolerance) | ...
        (Position_Data_sub(:,2) > reward_locations(2) - reward_zone_tolerance & ...
        Position_Data_sub(:,2) < reward_locations(2) + reward_zone_tolerance);
    stationary = Position_Data_sub(:,5) < drinking_speedThr;

    good_inds = drinking & stationary;

    x_NaN = nan(length(Position_Data_sub),1);
    x_NaN(good_inds) = 1;
    % find the longest stretch of stationary periods near the reward zone
    [boundaries,lengths] = compute_allSequences_NaNseparated_cm(x_NaN);

    % merge periods that are very close in time
    [boundaries,lengths] = compute_allSequences_NaNseparated_merge_time_only(Position_Data_sub(:,1),boundaries,deltThr*30000);

    % remove drinks less than min_drink_duration_thr
    boundaries(lengths*0.05 < min_drink_duration_thr,:) = [];
    lengths(lengths*0.05 < min_drink_duration_thr,:) = [];


    if isempty(boundaries)
        Reward_Epoch_Time_Boundaries_drinking(i,:) = [NaN NaN];
        continue
    end

    % Take the first drink period only
    boundaries = boundaries(1,:);
    drink_duration = lengths(1,:)*0.05;

% if drink_duration > 20
%     keyboard
% end
    Reward_Epoch_Time_Boundaries_drinking(i,:) = [Position_Data_sub(boundaries(1),1) Position_Data_sub(boundaries(2),1)];


    drinking_inds = find(Position_Data(:,1) >= Reward_Epoch_Time_Boundaries_drinking(i,1) & Position_Data(:,1) <= Reward_Epoch_Time_Boundaries_drinking(i,2));
    
    if plot_figure ==1
    plot(Position_Data(drinking_inds,1)./30000 - session_start_time/30000,Position_Data(drinking_inds,2),'.k')
    end
end





save('Behavior_Data','Pass_Info','Pass_Transitions','Reward_Epoch_Time_Boundaries_endzone','Reward_Epoch_Time_Boundaries','Reward_Epoch_Time_Boundaries_speed_thresholded','Reward_Epoch_Time_Boundaries_matched','Run_Epoch_Time_Boundaries','Reward_Epoch_Time_Boundaries_drinking')