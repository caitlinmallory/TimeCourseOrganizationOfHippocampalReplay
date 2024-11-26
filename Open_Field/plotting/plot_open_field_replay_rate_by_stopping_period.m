%% This script generates Fig. 2J-K: plots the rate of prospective or retrospective replays over time since stopping.


% Parameters for smoothing and window size
smooth_rate_plots = 1;           % If 1, smooth the replay rate plots
smoothing_sigma = 1;             % Smoothing window size in bins
start_time = 0;                  % Start time for the analysis (in seconds)
end_time = 14;                   % End time for the analysis (in seconds)
windowSize = 0.5;                % Bin size (in seconds)
include_trials_with_no_replays = 1;  % If 1, include trials with no replays; if 0, exclude them

% Handle trials with no replays based on the user's choice
if include_trials_with_no_replays == 1
    add_pseudo_replays_for_rate_analysis  % Add pseudo-replays for rate analysis if needed
    % Get trials with pseudo replays included
    pseudo_trials = (replay_with_pseudo_replays_for_trial(replay_with_pseudo_replays_for_trial.fake_replay_included_just_for_trial_entry == 1, :));
else
    replay_with_pseudo_replays_for_trial = replay;  % Use original replay data if no pseudo-replays are added
end

% Remove rows with missing drink period numbers from the data
replay_with_pseudo_replays_for_trial(isnan(replay_with_pseudo_replays_for_trial.drink_period_number), :) = [];

% Flags to analyze only home or away trials
home_only = 0;  % Set to 1 to look only at home trials
away_only = 0;  % Set to 1 to look only at away trials

% Flags for statistical comparisons
compare_means = 0;    % Set to 1 to compare means
compare_medians = 1;  % Set to 1 to compare medians

% Set the y-limits for rate and rate difference plots
ylim_rates = [0 0.25];         % Y-axis limits for rate plots
ylim_rate_diffs = [-0.1 0.20]; % Y-axis limits for rate differences

%% For figure 1N, plot rates of prospective and retrospective replay when the animal was at the home well or away well.
%% Data preparation for Figure 1N: analyze prospective and retrospective replays for home and away conditions
category1 = 'home';  % Label for home trials
category2 = 'away';  % Label for away trials
% Unique trial identifiers to categorize data
[C, ia, ic] = unique(replay_with_pseudo_replays_for_trial(:, {'unique_session_id', 'drink_period_number', 'home_event', 'rat_label', 'drink_period_time', 'laser_state_binary'}));
inds1 = find(C.home_event == 1);  % Indices for home trials
inds2 = find(C.home_event == 0);  % Indices for away trials

%% Data preparation for Figure 4: Analyze prospective and retrospective replays during laser on and off conditions (commented out)
% category1 = 'laser off';  % Label for laser off trials
% category2 = 'laser on';   % Label for laser on trials
% [C, ia, ic] = unique(replay_with_pseudo_replays_for_trial(:, {'unique_session_id', 'drink_period_number', 'home_event', 'rat_label', 'drink_period_time', 'laser_state_binary'}), 'stable');
% if home_only == 1
%     C = C(C.home_event == 1, :);  % Keep only home trials
% end
% if away_only == 1
%     C = C(C.home_event == 0, :);  % Keep only away trials
% end
% inds1 = find(C.laser_state_binary == 0);  % Indices for laser off
% inds2 = find(C.laser_state_binary == 1);  % Indices for laser on

%% Data preparation for Figure 3E: Analyze prospective and retrospective replays when the angle between future and past paths are tight (similar) or wide (diverging)
% category1 = 'tight';  % Label for tight angles
% category2 = 'wide';   % Label for wide angles
% [C, ia, ic] = unique(replay_with_pseudo_replays_for_trial(:, {'unique_session_id', 'drink_period_number', 'home_event', 'rat_label', 'drink_period_time', 'laser_state_binary', 'angle_between_past_future_trajectory'}));
% if home_only == 1
%     C = C(C.home_event == 1, :);  % Keep only home trials
% end
% if away_only == 1
%     C = C(C.home_event == 0, :);  % Keep only away trials
% end
% inds1 = find(C.angle_between_past_future_trajectory >= 0 & C.angle_between_past_future_trajectory <= 60);  % Indices for tight angles
% inds2 = find(C.angle_between_past_future_trajectory >= 120 & C.angle_between_past_future_trajectory <= 180);  % Indices for wide angles
%% Data preparation for Figure S3G: Analyze prospective and retrospective replays in the first or second half of the session
% category1 = 'first_half';  % Label for the first half of the session
% category2 = 'second_half'; % Label for the second half of the session
% [C, ia, ic] = unique(replay_with_pseudo_replays_for_trial(:, {'unique_session_id', 'drink_period_number', 'home_event', 'rat_label', 'drink_period_time', 'laser_state_binary'}));
% C.session_half = nan(height(C), 1);  % Initialize session_half column
% unique_sessions = unique(C.unique_session_id);  % Get unique sessions
% for i = 1:length(unique_sessions)
%     all_trials = find(C.unique_session_id == unique_sessions(i));  % Get all trials for this session
%     first_half_inds = all_trials(1:round(length(all_trials)/2));  % First half trials
%     second_half_inds = all_trials(round(length(all_trials)/2)+1:end);  % Second half trials
%     C.session_half(first_half_inds) = 1 * ones(length(first_half_inds), 1);  % Mark first half trials
%     C.session_half(second_half_inds) = 2 * ones(length(second_half_inds), 1);  % Mark second half trials
% end
% if home_only == 1
%     C = C(C.home_event == 1, :);  % Keep only home trials
% end
% if away_only == 1
%     C = C(C.home_event == 0, :);  % Keep only away trials
% end
% inds1 = find(C.session_half == 1);  % Indices for first half of session
% inds2 = find(C.session_half == 2);  % Indices for second half of session
%% Counts the number of each replay type (prospective/future, retrospective/past, or all replays) in 0.5 s bins
num_stopping_periods = height(C);  % Number of unique stopping periods
event_types = {'future', 'past', 'replay'};  % Event types to consider
bin_edges = start_time:windowSize:end_time;  % Define the time bins
bin_centers = (bin_edges(1:end-1)' + bin_edges(2:end)') / 2;  % Get the center of each time bin
bin_centers_plot = bin_centers;  % Store the bin centers for plotting

% Adjust bin centers if aligning to drink offset
if align_to_drink_offset == 1
    bin_centers_plot = -1.*bin_centers;  % Reverse time for alignment to drink offset
end

% Initialize structures to store results
combined_hists = struct();
combined_hists_rates = struct();
for i =1:length(event_types)
    combined_hists.(event_types{i}) = nan(num_stopping_periods,length(bin_edges)-1); % Initialize histogram arrays
    combined_hists_rates.(event_types{i}) = nan(num_stopping_periods,length(bin_edges)-1); % Initialize rate arrays
end

binned_data_by_stopping_period.replay = nan(num_stopping_periods,length(bin_edges)-1); % Initialize replay data array

% Loop through each stopping period and calculate event histograms
for i = 1:num_stopping_periods
    % find the data in this unique stopping period:
    replay_sub = replay_with_pseudo_replays_for_trial(replay_with_pseudo_replays_for_trial.unique_session_id == C.unique_session_id(i) & replay_with_pseudo_replays_for_trial.drink_period_number == C.drink_period_number(i),:);
    for n = 1:length(event_types)
        event_inds = find(replay_sub.(event_types{n})==1);
        if align_to_drink_onset==1
            % Align event times to drink onset
            event_hists =  histcounts(replay_sub.time_since_real_drink_onset(event_inds),bin_edges);
        elseif align_to_drink_offset==1
            % Align event times to drink offset
            event_hists =  histcounts(replay_sub.time_till_real_drink_offset(event_inds),bin_edges);
        end
        event_hists(replay_sub.drink_period_time(1) < bin_centers(:,1)) = nan;

        % length of stopping period:
        combined_hists.(event_types{n})(i,:) = event_hists;
        combined_hists_rates.(event_types{n})(i, :) = event_hists / windowSize;  % Convert to rates
    end
end

% smooth the rates on each lap, if requested
if smooth_rate_plots==1
    filter_length = smoothing_sigma*6;
    if(mod(filter_length,2)==0) % iseven
        filter_length = filter_length+1;
    end
    w = setUp_gaussFilt([1 filter_length ],smoothing_sigma);
    for i = 1:length(event_types)
        for j = 1:size(combined_hists_rates.future,1)

            rate_sub = combined_hists_rates.(event_types{i})(j,:);
            no_nan_inds = find(~isnan(rate_sub));
            rate_sub(isnan(rate_sub)) = [];
            combined_hists_rates.(event_types{i})(j,no_nan_inds) = conv(rate_sub,w,'same');
        end
    end
end

%% Generates figure 2J-K
colorlimits = [0 8]; % Defines the color limits for the colormap (minimum and maximum values).
num_laps_to_plot = 60; % Number of laps to plot (i.e., how many trials or laps).
num_bins_to_plot = 20; % Number of bins to use for time data (each bin represents a 0.5 second interval).

% Initialize arrays to store time spent in stopping periods for past and future
past_time_in_stopping_period_by_lap = nan(num_laps_to_plot,num_bins_to_plot); 
future_time_in_stopping_period_by_lap = nan(num_laps_to_plot,num_bins_to_plot);

% Define custom colormaps for visualizing the data
green_colormap = customcolormap(linspace(0,1,5), {'#014419','#1c7735','#5aae60','#a6db9d','#d7f1d6'});
purple_colormap = customcolormap(linspace(0,1,5), {'#410149','#762a84','#9b6fac','#c1a5cd','#e7d4e8'});

% Loop through the laps (trials) and extract data for future and past time spent in stopping periods
for lap = 1:num_laps_to_plot
    % Find the indices of the data for the current lap
    inds_0 = find(C.drink_period_number==lap);

    % Sum the time spent in the stopping period for the current lap (future and past)
    future_time_in_stopping_period_by_lap(lap,:) = nansum(combined_hists.('future')(inds_0,1:num_bins_to_plot),1);
    past_time_in_stopping_period_by_lap(lap,:) = nansum(combined_hists.('past')(inds_0,1:num_bins_to_plot),1);
end

% Create a new figure for plotting
fig1 = figure();
fig1.Position = [600 600 300 250]; % Set the position of the figure window
tiledlayout(2,2); % Create a 2x2 grid layout for subplots

% Plot the future time spent in stopping period
ax1 = nexttile(1); % Select the first tile in the 2x2 grid
if align_to_drink_offset == 1
    imagesc(fliplr(future_time_in_stopping_period_by_lap)); % Flip the data if aligning to drink offset
    xticks(linspace(0.5,20.5,6)); % Set x-axis ticks
    xticklabels({}) % Remove x-axis tick labels
else
    imagesc(future_time_in_stopping_period_by_lap); % Display the image for future time
    xticks(linspace(0.5,20.5,6)); % Set x-axis ticks
    xticklabels({}) % Remove x-axis tick labels
end
box off; % Remove the box around the plot
colormap(ax1, green_colormap); % Set the colormap for this plot
caxis(colorlimits); % Set the color axis limits
xtickangle(0); % Set the angle of the x-axis ticks
yticks([1,10,20,30,40,50,60]); % Set the y-axis ticks
ylabel('Trial'); % Label the y-axis
ax1.FontSize = 8; % Set font size for the axes


% Plot the past time spent in stopping period
ax3 = nexttile(3); % Select the third tile in the grid
if align_to_drink_offset == 1
    imagesc(fliplr(past_time_in_stopping_period_by_lap)); % Flip the data if aligning to drink offset
    xticks(linspace(0.5,20.5,6)); % Set x-axis ticks
    xticklabels(num2cell(0:2:10)); % Label x-axis with time intervals
    xlabel('Time till departure (s)'); % Label x-axis
else
    imagesc(past_time_in_stopping_period_by_lap); % Display the image for past time
    xticks(linspace(0.5,20.5,6)); % Set x-axis ticks
    xticklabels(num2cell(0:2:10)); % Label x-axis with time intervals
    xlabel('Time since arrival (s)'); % Label x-axis
end
box off; % Remove the box around the plot
colormap(ax3, purple_colormap); % Set the colormap for this plot
caxis(colorlimits); % Set the color axis limits
xtickangle(0); % Set the angle of the x-axis ticks
yticks([1,10,20,30,40,50,60]); % Set the y-axis ticks
ylabel('Trial'); % Label the y-axis
ax3.FontSize = 8; % Set font size for the axes

% Plot rate data (future and past time spent in stopping period)
ax2 = nexttile(2); % Select the second tile in the grid
x_mean = nanmean(combined_hists_rates.future); % Calculate mean rate for future data
x_err = nanstd(combined_hists_rates.future./sqrt((sum(~isnan(combined_hists_rates.future))))); % Calculate standard error for future data
if align_to_drink_offset == 1
    % Plot the future data, with optional smoothing
    if smooth_rate_plots == 1
        h = shadedErrorBar(bin_centers_plot((filter_length-1)/2:end), x_mean((filter_length-1)/2:end), x_err((filter_length-1)/2:end), 'lineprops', '-k'); hold on;
    else
        h = shadedErrorBar(bin_centers_plot, x_mean, x_err, 'lineprops', '-k'); hold on;
    end
    xticks(-10:2:0); % Set x-axis ticks for time till departure
    xticklabels([0 2 4 6 8 10]); % Label x-axis with appropriate time intervals
    xlim([-10 0]); % Set x-axis limits
else
    h = shadedErrorBar(bin_centers_plot, x_mean, x_err, 'lineprops', '-k'); hold on;
    xticks(0:2:10); % Set x-axis ticks for time since arrival
    xlim([0 10]); % Set x-axis limits
end
h.patch.FaceColor = color_future; % Set color for the future data plot
h.mainLine.Color = color_future; % Set the line color for the future plot
h.mainLine.LineWidth = 1; % Set the line width for the plot
h.edge(1).Color = 'none'; h.edge(2).Color = 'none'; % Remove edge lines around the shaded error
% Plot the past data
% Plot the past data
x_mean = nanmean(combined_hists_rates.past); % Calculate mean rate for past data
x_err = nanstd(combined_hists_rates.past./sqrt((sum(~isnan(combined_hists_rates.past))))); % Calculate standard error for past data
hold on;
if align_to_drink_offset == 1
    % Plot the past data, with optional smoothing
    if smooth_rate_plots == 1
        h = shadedErrorBar(bin_centers_plot((filter_length-1)/2:end), x_mean((filter_length-1)/2:end), x_err((filter_length-1)/2:end), 'lineprops', '-k'); hold on;
    else
        h = shadedErrorBar(bin_centers_plot, x_mean, x_err, 'lineprops', '-k'); hold on;
    end
    xticks(-10:2:0); % Set x-axis ticks for time till departure
    xlim([-10 0]); % Set x-axis limits
else
    h = shadedErrorBar(bin_centers_plot, x_mean, x_err, 'lineprops', '-k'); hold on;
    xticks(0:2:10); % Set x-axis ticks for time since arrival
    xlim([0 10]); % Set x-axis limits
end
h.patch.FaceColor = color_past; % Set color for the past data plot
h.mainLine.Color = color_past; % Set the line color for the past plot
h.mainLine.LineWidth = 1; % Set the line width for the plot
h.edge(1).Color = 'none'; h.edge(2).Color = 'none'; % Remove edge lines around the shaded error
xticklabels({''}); % Remove x-axis tick labels
ylabel('Events/s'); % Label y-axis
ylim(ylim_rates); % Set y-axis limits for rate data

% Plot the difference between future and past rates
ax4 = nexttile(4); % Select the fourth tile in the grid
x_mean = nanmean(combined_hists_rates.future - combined_hists_rates.past); % Calculate mean rate difference
x_err = nanstd(combined_hists_rates.future - combined_hists_rates.past)./sqrt((sum(~isnan((combined_hists_rates.future - combined_hists_rates.past))))); % Calculate standard error for rate difference
if align_to_drink_offset == 1
    % Plot the rate difference, with optional smoothing
    if smooth_rate_plots == 1
        h = shadedErrorBar(bin_centers_plot((filter_length-1)/2:end), x_mean((filter_length-1)/2:end), x_err((filter_length-1)/2:end), 'lineprops', '-k'); hold on;
    else
        h = shadedErrorBar(bin_centers_plot, x_mean, x_err, 'lineprops', '-k'); hold on;
    end
    xticks(-10:2:0); % Set x-axis ticks for time till departure
    xticklabels(num2cell(0:2:10)); % Label x-axis with time intervals
    xlim([-10 0]); % Set x-axis limits
    xlabel('Time till departure (s)'); % Label x-axis
else
    h = shadedErrorBar(bin_centers_plot, x_mean, x_err, 'lineprops', '-k'); hold on;
    xticks(0:2:10); % Set x-axis ticks for time since arrival
    xlim([0 10]); % Set x-axis limits
    xlabel('Time since arrival (s)'); % Label x-axis
end
ylim(ylim_rate_diffs); % Set y-axis limits for rate difference
h.patch.FaceColor = [0 0 0]; % Set color for the rate difference plot
h.mainLine.Color = [0 0 0]; % Set line color for the rate difference plot
h.mainLine.LineWidth = 1; % Set line width for the plot
h.edge(1).Color = 'none'; h.edge(2).Color = 'none'; % Remove edge lines around the shaded error
xtickangle(0); % Set the angle of the x-axis ticks
yline(0); % Add a horizontal line at y=0

% Perform statistical comparison (e.g., t-test or sign rank) to identify significant differences
if compare_means == 1
    sig_bins = ttest2(combined_hists_rates.future, combined_hists_rates.past); % T-test comparison
    sig_bins(isnan(sig_bins)) = 0; % Remove NaNs from the significance result
    sig_bins = logical(sig_bins); % Convert the result to a logical array
elseif compare_medians == 1
    sig_bins = zeros(length(bin_centers), 1); % Initialize array for significance results
    for bin = 1:length(bin_centers)
        [~, sig_bins(bin), ~] = signrank(combined_hists_rates.future(:, bin), combined_hists_rates.past(:, bin)); % Sign rank test for each bin
    end
    sig_bins = logical(sig_bins); % Convert the result to a logical array
end

% Plot significant bins on the figure
ylimit = gca().YLim; % Get the current y-axis limits
plot(bin_centers_plot(sig_bins), ylimit(2) * ones(sum(sig_bins)), '.k'); % Plot significant bins
ylabel('Events/s'); % Label y-axis


% Save the figure
set(gcf, 'Color', 'white', 'Renderer', 'painters'); % Set figure background color and renderer
set(gcf, 'PaperPositionMode', 'auto'); % Set automatic paper position mode
fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures'; % Define path to save the figure
saveas(gcf, fullfile(fig_path, 'Main_open_field_rate_summary'), 'pdf'); % Save the figure as a PDF

% Calculate the number of non-NaN values in the future data
all_n = sum(~isnan(combined_hists_rates.future)); 
max_n = max(all_n(bin_centers <= 10)); % Maximum number of non-NaN values for bins <= 10
min_n = min(all_n(bin_centers <= 10)); % Minimum number of non-NaN values for bins <= 10

%% 
% Generate a 2 by 2 bar graph, collapsed across time (Fig 2O, Fig 4F).
plot_fig2_rate_summary 
% Plot home prospective versus home retrospective (Generates Fig 2N, Fig3E,
% Fig 4E).
plot_figure2_prospective_retrospective_and_diff

