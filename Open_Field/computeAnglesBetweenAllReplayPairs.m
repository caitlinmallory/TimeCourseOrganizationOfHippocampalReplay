% This code computes the angles between all pairs of replays within a
% stopping period.

clearvars -except dayFiles day directory rat windows hand_clustered_only

load Experiment_Information.mat
load Analysis_Information.mat
load Position_Data.mat
load replayEvents

% Load the true drinking periods, which mark the time intervals of interest
load true_drink_periods.mat

% Scale the position data into spatial bins based on the given number of bins and edges
Position_Data_scaled = compute_locsToBins(Position_Data,numSpatialBins,x_edges,y_edges);

% Pull out all unique run sessions in this day:
run_segments = [];
for i = 1:length(Experiment_Information.Segments)
    if ismember(11,Experiment_Information.Segments(i).Flags)
        run_segments = [run_segments; i]; % Add segment index to run_segments
    end
end

% One segment = one run with a consistent home well. Bolt had up to 2 homed
% wells, 4 sessions per day. CM1 had 2 home well/2 sessions per day.
for segment_count = 1:length(run_segments)
    rat_position_during_stops = nan(length(true_drink_periods_summary(segment_count).true_drink_periods),2);

    % For each stopping period, extract the average location of the rat
    for i = 1:length(true_drink_periods_summary(segment_count).true_drink_periods)
        ratPos = compute_dataTemporalConcatenation(Position_Data_scaled,[true_drink_periods_summary(segment_count).true_drink_periods(i,1) true_drink_periods_summary(segment_count).true_drink_periods(i,2)]);
        if isempty(ratPos)
            continue
        end
        rat_position_during_stops(i,:) = ratPos(1,2:3);
    end

    % Get the segment details: segment ID, times, and decoder object
    segment_id = run_segments(segment_count);
    segment_times = Experiment_Information.Segments(segment_id).Times;
    decoder = Experiment_Information.Segments(segment_id).Decoder;

    % For pull out the drink times
    true_drink_periods = true_drink_periods_summary(segment_count).true_drink_periods; % times

    % Load and filter the replay events for this segment (home/away events)
    replay_events = struct2table(decoder_replay(decoder).replayEvents);
    replay_events = replay_events(replay_events.away_event==1 | replay_events.home_event==1,:);

    % Initialize arrays to store cleaned and upsampled replays
    all_replays_NaNremoved = cell(height(replay_events),1);
    all_replays_upsampled = cell(height(replay_events),1);

    % pull out all replays at the start, upsample to be 10000 points long.
    for i = 1:height(replay_events)
        % Remove NaN values from the replay data and store cleaned replay
        replay_NaNremoved = [mean(replay_events.timeBins{i},2) replay_events.replay{i}];
        replay_NaNremoved(replay_events.indNaN{i},:) = [];
        all_replays_NaNremoved{i} = replay_NaNremoved;
        if size(replay_NaNremoved,1) == 1
            all_replays_upsampled{i} = nan;
        else
            all_replays_upsampled{i} = interp1(1:size(replay_NaNremoved,1),[replay_NaNremoved(:,2) replay_NaNremoved(:,3)], linspace(1,size(replay_NaNremoved,1),10000));
        end
    end

    % Process each replay event for intersection with concentric rings around the rat
    for i = 1:height(replay_events)
        tic
        replay = all_replays_upsampled{i}; % Get the current replay trajectory
        replay_timePoint = replay_events.timePoints(i); % Get the timestamp of the replay event

        % Find the true drink period that the replay occurred during
        current_drink_period_ind = find(replay_timePoint>=true_drink_periods(:,1) & replay_timePoint<=true_drink_periods(:,2));
        if isempty(current_drink_period_ind)  % If no exact match, find closest
            [~,current_drink_period_ind] = min(abs(replay_timePoint - true_drink_periods(:,1)));
        end

        % Store the rat's position during the current stopping period
        replay_events.ratPos(i,:) = rat_position_during_stops(current_drink_period_ind,:);


        % Compute the intersection points of the replay trajectory with concentric circles
        replay_events.intersections{i} =  computePointsOfIntersectionWithConcentricCircles(replay,replay_events.ratPos(i,:));
        disp([num2str(i) '/' num2str(height(replay_events))]) % display progress
        toc
    end % end looping through detected replays

    %%

    % Initialize variables to store data for angle comparisons
    all_angles=[];
    all_comparison_inds = [];
    all_time_differences = [];
    all_times_into_stopping_period = [];
    all_dispersion = [];
    all_start_distance_from_rat = [];

    % For each stopping period, comparing the angle of all replays
    for i = 1:length(true_drink_periods_summary(segment_count).true_drink_periods)
        % Get replay events that occurred during the current stopping period
        replays_sub = replay_events(replay_events.timePoints(:,1)>=true_drink_periods_summary(segment_count).true_drink_periods(i,1) & replay_events.timePoints(:,1)<=true_drink_periods_summary(segment_count).true_drink_periods(i,2),:);
        replays_sub_ind = find(replay_events.timePoints(:,1)>=true_drink_periods_summary(segment_count).true_drink_periods(i,1) & replay_events.timePoints(:,1)<=true_drink_periods_summary(segment_count).true_drink_periods(i,2));

        % Skip if fewer than 2 replays exist during this period
        if (isempty(replays_sub_ind) | length(replays_sub_ind)<2)
            continue
        end

        % Create pairwise comparisons between replays
        comparisons = nchoosek(1:height(replays_sub),2);
        comparison_inds = [replays_sub_ind(comparisons(:,1)) replays_sub_ind(comparisons(:,2))];


        % Calculate angles between replay pairs
        ratPos = replays_sub.ratPos(1,:); % Get the rat's position during the stop
        angles = nan(size(comparison_inds,1),90); % Initialize matrix to store angles

        % Loop through each pair of replays and calculate angles
        for j = 1:size(comparison_inds,1)
            ind1 = comparison_inds(j,1); ind2 = comparison_inds(j,2);

            % Get the replay trajectory vectors relative to the rat's position
            u = [replay_events.intersections{ind1}(:,1)-ratPos(1) replay_events.intersections{ind1}(:,2)-ratPos(2)];
            v = [replay_events.intersections{ind2}(:,1)-ratPos(1) replay_events.intersections{ind2}(:,2)-ratPos(2)];

            % Compute the angle between the two replay vectors
            for k = 1:90
                angles(j,k) = rad2deg(atan2(u(k,1)*v(k,2)-u(k,2)*v(k,1),u(k,1)*v(k,1)+u(k,2)*v(k,2)));
            end
        end

        all_angles = [all_angles; angles];
        all_comparison_inds = [all_comparison_inds; comparison_inds];
        all_time_differences = [all_time_differences; (abs(replay_events.timePoints(comparison_inds(:,1),1)-replay_events.timePoints(comparison_inds(:,2),1)))./30000];
        all_times_into_stopping_period = [all_times_into_stopping_period; [replay_events.time_since_real_drink_onset(comparison_inds(:,1)) replay_events.time_since_real_drink_onset(comparison_inds(:,2))]];
        all_dispersion = [all_dispersion; [replay_events.dispersion(comparison_inds(:,1)) replay_events.dispersion(comparison_inds(:,2))]];
        all_start_distance_from_rat = [all_start_distance_from_rat; [replay_events.startDistFromRat(comparison_inds(:,1)) replay_events.startDistFromRat(comparison_inds(:,2))]];
    end

    % store the results for this segment
    angles_between_replays(segment_count).replays = replay_events;
    angles_between_replays(segment_count).all_angles = all_angles;
    angles_between_replays(segment_count).all_comparison_inds = all_comparison_inds;
    angles_between_replays(segment_count).all_time_differences = all_time_differences;
    angles_between_replays(segment_count).all_times_into_stopping_period = all_times_into_stopping_period;
    angles_between_replays(segment_count).all_dispersion = all_dispersion;
    angles_between_replays(segment_count).all_start_distance_from_rat = all_start_distance_from_rat;

end

% Save the computed angles between replays to a file for later analysis
save('angles_between_replays','angles_between_replays')

