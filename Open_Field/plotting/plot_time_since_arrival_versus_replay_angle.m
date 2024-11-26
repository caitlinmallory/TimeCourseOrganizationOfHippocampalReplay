% Define file path where figures will be saved
fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures';

% Flag for plotting shuffled data (currently set to no plot of shuffled data)
plot_shuffles = 0;

% Uncomment the following lines to define custom color maps (currently not used)
% green_colormap = customcolormap(linspace(0,1,100), {'#014419','#1c7735','#5aae60','#a6db9d','#d7f1d6'});
% purple_colormap = customcolormap(linspace(0,1,100), {'#410149','#762a84','#9b6fac','#c1a5cd','#e7d4e8'});

% Set flag to normalize plots (currently no normalization)
normalize_plots = 0;

% Define color limits for plotting (main limits and normalized limits)
color_limits_normalized = [2 11];
color_limits = [0 200]; % Main limits
% color_limits = [0 40]; % Alternative limits (for supplementary data)

% Flag for plotting mean angular histograms (currently set to not plot)
plot_mean_angular_hist = 0;

% Define time bins for analysis (two time intervals for comparison)
time_bin_1 = [0 2];  % First time window (0 to 2)
time_bin_2 = [4 10]; % Second time window (4 to 10)

% Flags to plot angular histograms for the defined time bins (both are set to plot)
plot_time_bin_1_angular_hist = 1;
plot_time_bin_2_angular_hist = 1;

% Define bin sizes and edges for the time analysis
bin_start = (start_time:windowShift:(end_time-windowSize))'; % Start of time bins
bin_end = bin_start + windowSize; % End of time bins

% Define time bin edges and calculate bin centers
bin_edges = [bin_start bin_end];
time_bin_centers = mean(bin_edges,2);

% Define angle bins for angular displacement analysis (currently 19 bins)
angle_bins = linspace(0, 180, 19);  % 19 bins between 0 and 180 degrees
angle_bin_centers = mean([angle_bins(1:end-1)' angle_bins(2:end)'],2); % Center of each angle bin

% Initialize matrices to store angle distributions for different conditions
replay_past_path_angle_map = nan(length(time_bin_centers), length(angle_bin_centers));
replay_future_path_angle_map = nan(length(time_bin_centers), length(angle_bin_centers));

% If shuffling data is to be plotted, initialize shuffled angle matrices
replay_shuffled_past_path_angle_map = nan(length(time_bin_centers), length(angle_bin_centers));
replay_shuffled_future_path_angle_map = nan(length(time_bin_centers), length(angle_bin_centers));
replay_shuffled_path_angle_map = nan(length(time_bin_centers), length(angle_bin_centers));
replay_shuffled_path_all_angle_map = nan(length(time_bin_centers), length(angle_bin_centers));

% Loop through time bins to calculate angular displacement for each time point
for i = 1:length(time_bin_centers)
    replay_past_angles = binned_properties.meanAngDisplacement_pastPath(i,:).data{:};
    replay_future_angles = binned_properties.meanAngDisplacement_futPath(i,:).data{:};
    
    % If shuffling is enabled, calculate shuffled angles for comparison
    if plot_shuffles == 1
        replay_shuffled_past_angles = binned_properties.mean_shuffle_past(i,:).data{:};
        replay_shuffled_future_angles = binned_properties.mean_shuffle_future(i,:).data{:};
        replay_shuffled_angles = binned_properties.mean_shuffle_all(i,:).data{:};
        replay_shuffled_all_angles = binned_properties.all_shuffled_paths(i,:).data{:};
    end

    % Loop through angle bins to count the occurrences of angles within each bin
    for j = 1:length(angle_bins)-1
        replay_past_path_angle_map(i,j) = sum(replay_past_angles >= angle_bins(j) & replay_past_angles < angle_bins(j+1));
        replay_future_path_angle_map(i,j) = sum(replay_future_angles >= angle_bins(j) & replay_future_angles < angle_bins(j+1));

        % If shuffling is enabled, count occurrences for shuffled data
        if plot_shuffles == 1
            replay_shuffled_past_path_angle_map(i,j) = sum(replay_shuffled_past_angles >= angle_bins(j) & replay_shuffled_past_angles < angle_bins(j+1));
            replay_shuffled_future_path_angle_map(i,j) = sum(replay_shuffled_future_angles >= angle_bins(j) & replay_shuffled_future_angles < angle_bins(j+1));
            replay_shuffled_path_angle_map(i,j) = sum(replay_shuffled_angles >= angle_bins(j) & replay_shuffled_angles < angle_bins(j+1));
            replay_shuffled_path_all_angle_map(i,j) = sum(replay_shuffled_all_angles >= angle_bins(j) & replay_shuffled_all_angles < angle_bins(j+1));
        end
    end
end

% Calculate the total number of events and percentage of angles for past and future
total_number_events_past = sum(replay_past_path_angle_map,2);
total_number_angles_past = sum(replay_past_path_angle_map,1);
pcnt_angles_past = total_number_angles_past ./ sum(total_number_events_past);

total_number_events_future = sum(replay_future_path_angle_map,2);
total_number_angles_future = sum(replay_future_path_angle_map,1);
pcnt_angles_future = total_number_angles_future ./ sum(total_number_events_future);

% If shuffling is enabled, calculate for shuffled data as well
if plot_shuffles == 1
    total_number_events_shuffled_past = sum(replay_shuffled_past_path_angle_map,2);
    total_number_angles_shuffled_past = sum(replay_shuffled_past_path_angle_map,1);
    pcnt_angles_shuffled_past = total_number_angles_shuffled_past ./ sum(total_number_events_shuffled_past);

    total_number_events_shuffled_future = sum(replay_shuffled_future_path_angle_map,2);
    total_number_angles_shuffled_future = sum(replay_shuffled_future_path_angle_map,1);
    pcnt_angles_shuffled_future = total_number_angles_shuffled_future ./ sum(total_number_events_shuffled_future);

    total_number_events_shuffled_path = sum(replay_shuffled_path_angle_map,2);
    total_number_angles_shuffled = sum(replay_shuffled_path_angle_map,1);
    pcnt_angles_shuffled_path = total_number_angles_shuffled ./ sum(total_number_events_shuffled_path);

    total_number_events_shuffled_all_path = sum(replay_shuffled_path_all_angle_map,2);
    total_number_angles_shuffled_all = sum(replay_shuffled_path_all_angle_map,1);
    pcnt_angles_shuffled_all_path = total_number_angles_shuffled_all ./ sum(total_number_events_shuffled_all_path);
end

% Calculate total number of angles for past and future in specific time bins
total_number_angles_past_time_bin_1 = sum(replay_past_path_angle_map(time_bin_centers >= time_bin_1(1) & time_bin_centers <= time_bin_1(2), :), 1);
total_number_angles_past_time_bin_2 = sum(replay_past_path_angle_map(time_bin_centers >= time_bin_2(1) & time_bin_centers <= time_bin_2(2), :), 1);

% Normalize the angle distribution for each time bin
pcnt_angles_past_time_bin_1 = total_number_angles_past_time_bin_1 ./ sum(total_number_angles_past_time_bin_1);
pcnt_angles_past_time_bin_2 = total_number_angles_past_time_bin_2 ./ sum(total_number_angles_past_time_bin_2);

% Repeat the process for future angles in the specified time bins
total_number_angles_future_time_bin_1 = sum(replay_future_path_angle_map(time_bin_centers <= 3, :), 1) ./ sum(sum(replay_future_path_angle_map(time_bin_centers <= 3, :), 1));
total_number_angles_future_time_bin_2 = sum(replay_future_path_angle_map(time_bin_centers >= time_bin_2(1) & time_bin_centers <= time_bin_2(2), :), 1);
pcnt_angles_future_time_bin_1 = total_number_angles_future_time_bin_1 ./ sum(total_number_angles_future_time_bin_1);
pcnt_angles_future_time_bin_2 = total_number_angles_future_time_bin_2 ./ sum(total_number_angles_future_time_bin_2);

% If shuffling is enabled, calculate for shuffled data in the time bins as well
if plot_shuffles == 1
    % Past angles in specific time bins for shuffled data
    total_number_angles_shuffled_past_time_bin_1 = sum(replay_shuffled_past_path_angle_map(time_bin_centers >= time_bin_1(1) & time_bin_centers <= time_bin_1(2), :), 1);
    total_number_angles_shuffled_past_time_bin_2 = sum(replay_shuffled_past_path_angle_map(time_bin_centers >= time_bin_2(1) & time_bin_centers <= time_bin_2(2), :), 1);
    pcnt_angles_shuffled_past_time_bin_1 = total_number_angles_shuffled_past_time_bin_1 ./ sum(total_number_angles_shuffled_past_time_bin_1);
    pcnt_angles_shuffled_past_time_bin_2 = total_number_angles_shuffled_past_time_bin_2 ./ sum(total_number_angles_shuffled_past_time_bin_2);

    % Future angles in specific time bins for shuffled data
    total_number_angles_shuffled_future_time_bin_1 = sum(replay_shuffled_future_path_angle_map(time_bin_centers >= time_bin_1(1) & time_bin_centers <= time_bin_1(2), :), 1);
    total_number_angles_shuffled_future_time_bin_2 = sum(replay_shuffled_future_path_angle_map(time_bin_centers >= time_bin_2(1) & time_bin_centers <= time_bin_2(2), :), 1);
    pcnt_angles_shuffled_future_time_bin_1 = total_number_angles_shuffled_future_time_bin_1 ./ sum(total_number_angles_shuffled_future_time_bin_1);
    pcnt_angles_shuffled_future_time_bin_2 = total_number_angles_shuffled_future_time_bin_2 ./ sum(total_number_angles_shuffled_future_time_bin_2);

    % Path angles in specific time bins for shuffled data
    total_number_angles_shuffled_path_time_bin_1 = sum(replay_shuffled_path_angle_map(time_bin_centers >= time_bin_1(1) & time_bin_centers <= time_bin_1(2), :), 1);
    total_number_angles_shuffled_path_time_bin_2 = sum(replay_shuffled_path_angle_map(time_bin_centers >= time_bin_2(1) & time_bin_centers <= time_bin_2(2), :), 1);
    pcnt_angles_shuffled_path_time_bin_1 = total_number_angles_shuffled_path_time_bin_1 ./ sum(total_number_angles_shuffled_path_time_bin_1);
    pcnt_angles_shuffled_path_time_bin_2 = total_number_angles_shuffled_path_time_bin_2 ./ sum(total_number_angles_shuffled_path_time_bin_2);

    % All path angles in specific time bins for shuffled data
    total_number_angles_shuffled_path_all_time_bin_1 = sum(replay_shuffled_path_all_angle_map(time_bin_centers >= time_bin_1(1) & time_bin_centers <= time_bin_1(2), :), 1);
    total_number_angles_shuffled_path_all_time_bin_2 = sum(replay_shuffled_path_all_angle_map(time_bin_centers >= time_bin_2(1) & time_bin_centers <= time_bin_2(2), :), 1);
    pcnt_angles_shuffled_path_all_time_bin_1 = total_number_angles_shuffled_path_all_time_bin_1 ./ sum(total_number_angles_shuffled_path_all_time_bin_1);
    pcnt_angles_shuffled_path_all_time_bin_2 = total_number_angles_shuffled_path_all_time_bin_2 ./ sum(total_number_angles_shuffled_path_all_time_bin_2);
end

%%
plot_figure2_angle_distributions

