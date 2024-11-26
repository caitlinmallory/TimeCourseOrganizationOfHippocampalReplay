% Figure 3J

% Must have run do_combine_replayPropertiesAcrossSessions first to load and
% filter data.

% Flag to control whether to correct for multiple comparisons using a correction factor
correct_for_multiple_comparisons = 0;

% Set the default axes colors to black
set(groot, {'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'}, {'k','k','k'})

% Specify the path to save the figure for the manuscript
fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures';

% Define restricted stop times for the analysis (start and end times)
restricted_stop_times = [0 inf];

% % Look at all forward/reverse congruent replays within a stopping period where laser is off
% Filter data where laser is off (laser_state == 0) and time_since_reward_zone_entry is within the restricted stop times
data = t_congruent_replay_og(t_congruent_replay_og.laser_state == 0 & ...
    t_congruent_replay_og.time_since_reward_zone_entry >= restricted_stop_times(1) & t_congruent_replay_og.time_since_reward_zone_entry <= restricted_stop_times(2),:);

% Calculate the ratio of forward direction replays
ratio = (sum(data.direction == 1) / height(data));

% Create a unique list of sessions and pass numbers
[C, ia, ic] = unique(data(:, {'session_str', 'pass_number'}));

% Initialize variables for time differences, replay types, and session information
all_time_differences = [];
all_replay_is_different_type = [];
all_session_id = {};
all_lap = [];
all_comparison_inds = [];

% Loop through each session and pass number combination (lap)
for lap = 1:height(C)
    % Extract session and lap information
    this_session = C.session_str{lap};
    this_lap = C.pass_number(lap);

    % Subset data for the current session and lap
    data_sub = data(strcmp(data.session_str, C.session_str{lap}) == 1 & data.pass_number == C.pass_number(lap), :);
    replays_sub_ind = find(strcmp(data.session_str, C.session_str{lap}) == 1 & data.pass_number == C.pass_number(lap));

    % Skip if there are fewer than 2 replays in the current subset
    if height(data_sub) < 2
        continue
    end

    % Generate all pairwise comparisons for replays in the stopping period
    comparisons = nchoosek(1:height(data_sub), 2);
    all_comparison_inds = [all_comparison_inds; [replays_sub_ind(comparisons(:,1)) replays_sub_ind(comparisons(:,2))]];

    % Initialize arrays to store time differences and replay type differences
    lap_time_differences = nan(size(comparisons, 1), 1);
    lap_replay_is_different_type = nan(size(comparisons, 1), 1);

    % Loop through each pairwise comparison
    for i = 1:size(comparisons, 1)
        % Calculate the time difference between the two replays
        lap_time_differences(i) = abs(data_sub.timePoints(comparisons(i, 1), 1) - data_sub.timePoints(comparisons(i, 2), 1)) / 30000;

        % Check if the two replays are of different types (forward vs reverse direction)
        if (data_sub.direction(comparisons(i, 1)) == 1 && data_sub.direction(comparisons(i, 2)) == 2) || ...
                (data_sub.direction(comparisons(i, 1)) == 2 && data_sub.direction(comparisons(i, 2)) == 1)
            lap_replay_is_different_type(i) = 1; % Different map types
        else
            lap_replay_is_different_type(i) = 0; % Same map type
        end
    end

    % Accumulate the time differences and replay type differences for all laps
    all_time_differences = [all_time_differences; lap_time_differences];
    all_replay_is_different_type = [all_replay_is_different_type; lap_replay_is_different_type];
    all_session_id = [all_session_id; repmat({this_session}, [size(comparisons, 1), 1])];
    all_lap = [all_lap; repmat(this_lap, [size(comparisons, 1), 1])];
end

% Create a table to store the results of comparisons
switches = table();
switches.time_diff = all_time_differences;
switches.session_id = all_session_id;
switches.lap = all_lap;
switches.time_into_stopping_period = [data.time_since_reward_zone_entry(all_comparison_inds(:,1)) data.time_since_reward_zone_entry(all_comparison_inds(:,2))] / 60; % Time in minutes
switches.time_into_session = [data.time_in_session(all_comparison_inds(:,1)) data.time_in_session(all_comparison_inds(:,2))] / 60; % Time in minutes
switches.direction = [data.direction(all_comparison_inds(:,1)) data.direction(all_comparison_inds(:,2))];
switches.map = [data.best_map(all_comparison_inds(:,1)) data.best_map(all_comparison_inds(:,2))];
switches.different = all_replay_is_different_type;

% Find comparisons where the time difference is less than 1.5 seconds (for visualization)
example_inds = find(all_time_differences < 1.5);
examples = switches(example_inds, :);

%% Plot preparation
windowSize = 0.4;   % Size of the time window (in minutes)
windowShift = 0.1;  % Shift between consecutive time windows
start_time = 0;     % Starting time for binning
end_time = 4;       % Ending time for binning

% Create time bins
bin_start = (start_time:windowShift:(end_time - windowSize))';
bin_end = bin_start + windowSize;
bin_edges = [bin_start bin_end];
time_bin_centers = mean(bin_edges, 2);

% Initialize arrays to store data for each time bin
data_fraction = nan(length(time_bin_centers), 1);
data_n = nan(length(time_bin_centers), 1);
data_n_forward = nan(length(time_bin_centers), 1);
data_n_reverse = nan(length(time_bin_centers), 1);
data_n_switches = nan(length(time_bin_centers), 1);

% Number of bootstrap iterations
num_boots = 1000;
bootstrapped_data_fraction = nan(length(bin_edges), num_boots);

% Loop through each time bin and calculate various statistics
for bin = 1:length(bin_edges)
    % Filter data to include only comparisons within the current time bin
    switches_sub = switches(switches.time_diff >= bin_edges(bin, 1) & switches.time_diff < bin_edges(bin, 2), :);

    % Count forward and reverse events
    [~, ib, ~] = unique(switches_sub.time_into_stopping_period);
    data_n_forward(bin) = sum(switches_sub.direction(ib) == 1);
    data_n_reverse(bin) = sum(switches_sub.direction(ib) == 2);

    % Calculate the fraction of "different" replay types in this bin
    data_sub = switches_sub.different;
    data_n(bin) = height(data_sub);
    data_n_switches(bin) = sum(data_sub);
    data_fraction(bin) = sum(data_sub) / height(switches_sub);

    % Bootstrap resampling to estimate confidence intervals
    if ~isempty(data_sub)
        for j = 1:num_boots
            rand_inds = randi(length(data_sub), [length(data_sub), 1]);
            bootstrapped_data_fraction(bin, j) = sum(data_sub(rand_inds)) / length(data_sub);
        end
    end
end

% Number of valid data points in the "unbinned" data
data_n_unbinned = height(switches(switches.time_diff > 0 & switches.time_diff <= 4, :));
min_data_n = min(data_n);
max_data_n = max(data_n);

%% Shuffle the data to create a null distribution
num_shuffles = 500;
shuffled_data_fraction = nan(length(bin_edges), num_shuffles);
switches_shuffle = switches;
for shuffle = 1:num_shuffles
    % Shuffle the time differences between comparisons
    switches_shuffle.time_diff = switches.time_diff(randperm(height(switches), height(switches)));

    % Loop through each time bin and calculate shuffled data fractions
    for bin = 1:length(bin_edges)
        switches_sub = switches_shuffle(switches_shuffle.time_diff >= bin_edges(bin, 1) & switches_shuffle.time_diff < bin_edges(bin, 2), :);

        % Count forward and reverse events for shuffled data
        [~, ib, ~] = unique(switches_sub.time_into_stopping_period);
        shuffled_data_n_forward(bin) = sum(switches_sub.direction(ib) == 1);
        shuffled_data_n_reverse(bin) = sum(switches_sub.direction(ib) == 2);

        % Calculate the fraction of "different" replay types in the shuffled data
        shuffled_data_sub = switches_sub.different;
        shuffled_data_n(bin, shuffle) = height(shuffled_data_sub);
        shuffled_data_n_switches(bin, shuffle) = sum(shuffled_data_sub);
        shuffled_data_fraction(bin, shuffle) = sum(shuffled_data_sub) / height(switches_sub);
    end
end

% Calculate p-values for the observed data vs. shuffled data
pvals = nan(length(time_bin_centers), 1);
for i = 1:length(time_bin_centers)
    right_side = sum(shuffled_data_fraction(i, :) > data_fraction(i));
    left_side = sum(shuffled_data_fraction(i, :) < data_fraction(i));
    count = min(right_side, left_side);
    pvals(i) = (count + 1) / (num_shuffles + 1);
end

% Apply multiple comparisons correction if specified
alpha = 0.05;
if correct_for_multiple_comparisons == 1
    alpha = alpha / size(shuffled_data_fraction, 1);
end

% Calculate shuffled data quantiles for significance testing
shuffled_low = quantile(shuffled_data_fraction, alpha / 2, 2);
shuffled_high = quantile(shuffled_data_fraction, 1 - alpha / 2, 2);

% Identify significant bins based on data fraction being outside the shuffled confidence interval
sig_bins = find(data_fraction > shuffled_high | data_fraction < shuffled_low);

%% Calculate bootstrapped confidence intervals for the observed data
low = quantile(bootstrapped_data_fraction, alpha / 2, 2);
high = quantile(bootstrapped_data_fraction, 1 - alpha / 2, 2);
err = [abs(data_fraction - low) abs(data_fraction - high)];

% Create the plot
figure('Position', [2313 876 175 125])  % Set figure size
h = shadedErrorBar(time_bin_centers, data_fraction .* 100, err .* 100); % Plot with shaded error bars
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';

% Customize plot labels and appearance
xlabel('Time between replays (s)');
ylabel('% Opposite');
ylim([0 100]);
xticks(0:1:end_time);
xtickangle(0);
hold on;

% Plot shuffled confidence intervals and significant bins
plot(time_bin_centers, shuffled_low .* 100, '--k');
plot(time_bin_centers, shuffled_high .* 100, '--k');
plot(time_bin_centers(sig_bins), 100 * ones(length(sig_bins), 1), '.k');

% Make the plot square
axis square;
