%% Plots angle distributions over time relative to past and future paths, as in Fig 2G

% Flag for normalizing the plots (set to 0 for no normalization)
normalize_plots = 0;

% Create a new figure with a specified size and position on screen
figure('Position', [552 561 200 250])

% Create a tiled layout with 2 rows and 3 columns for the subplots, and tight spacing between tiles
tiledlayout(2, 3, 'TileSpacing', 'tight')

% Plot the first subplot (spanning 2 columns)
nexttile(1, [1 2])

% Check if normalization is enabled
if normalize_plots == 1
    % Plot normalized future path angle map (percentage of each angle category)
    imagesc(100 .* ((replay_future_path_angle_map ./ sum(replay_future_path_angle_map, 2)))');
    caxis(color_limits_normalized)  % Apply the color axis limits for normalized data
else
    % Plot the raw future path angle map without normalization
    imagesc(replay_future_path_angle_map')
end

% Apply the 'magma' colormap for the plot
colormap(magma)

% Set the color axis for the plot (this controls the color scaling)
caxis(color_limits)

% Set the Y-axis to display with the normal direction (top to bottom)
set(gca, 'YDir', 'normal')

% Customize the Y-axis ticks (display every 6th angle bin)
yticks(0.5:6:(length(angle_bin_centers) + 1))
yticklabels(num2str(angle_bins(1:6:end)'))

% Customize the X-axis ticks (based on time, adjusted for window shift)
xticks(0.5:(2 / windowShift):length(time_bin_centers) + 1)

% Set X-axis limits to display time from the start to end
xlim([0 + windowShift, (10 / windowShift + windowShift)])

% Label the X-axis with time labels (every 2 seconds)
xticklabels(num2str([start_time:2:end_time]'))
xtickangle(0)  % Set the angle of the X-axis tick labels

% Label the axes
xlabel('Time since arrival (s)')
ylabel('|Displacement|')

% Draw a horizontal line at y=2.5 for reference
hold on
h = hline(2.5); 
h.Color = [1 1 1];  % White color for the line

% Plot the third subplot: Mean angular histogram for future path angles
nexttile(3)
hold on
% If the flag is set, plot the mean angular histogram for the entire future path
if plot_mean_angular_hist == 1
    plot(pcnt_angles_future, angle_bin_centers .* 100, 'b', 'linewidth', 2); 
    hold on;
end

% If the flag for time bin 1 is set, plot the angular histogram for future path in the first time bin
if plot_time_bin_1_angular_hist == 1
    plot(pcnt_angles_future_time_bin_1 .* 100, angle_bin_centers, 'k', 'linewidth', 2)
end

% If the flag for time bin 2 is set, plot the angular histogram for future path in the second time bin
if plot_time_bin_2_angular_hist == 1
    plot(pcnt_angles_future_time_bin_2 .* 100, angle_bin_centers, 'color', [0.6 0.6 0.6], 'linewidth', 2); 
    hold on
end

% Turn off the box around the plot
box off

% Set the X-axis limit for the angular histogram
xlim([0 10])

% Plot the fourth subplot (spanning 2 columns) for past path angles
nexttile(4, [1 2])
if normalize_plots == 1
    % Plot normalized past path angle map
    imagesc((100 .* (replay_past_path_angle_map ./ sum(replay_past_path_angle_map, 2)))');
    caxis(color_limits_normalized)  % Apply color limits for normalized data
else
    % Plot the raw past path angle map
    imagesc(replay_past_path_angle_map')
end

% Apply the 'magma' colormap for the plot
colormap(magma)
caxis(color_limits)  % Apply color axis limits
set(gca, 'YDir', 'normal')  % Set Y-axis direction
yticks(0.5:6:(length(angle_bin_centers) + 1))  % Set Y-ticks (every 6th angle bin)
yticklabels(num2str(angle_bins(1:6:end)'))  % Display angle bin labels
xticks(0.5:(2 / windowShift):length(time_bin_centers) + 1)  % Set X-ticks for time bins
xlim([0 + windowShift, (10 / windowShift + windowShift)])  % Set X-axis limits
xticklabels(num2str([start_time:2:end_time]'))  % Label the X-axis with time
xtickangle(0)  % Set X-axis tick angle
xlabel('Time since arrival (s)')  % X-axis label
ylabel('|Displacement|')  % Y-axis label

% Draw a horizontal line at y=2.5 for reference
hold on
h = hline(2.5); 
h.Color = [1 1 1];  % White color for the line

% Plot the sixth subplot: Mean angular histogram for past path angles
nexttile(6)
% If the flag is set, plot the mean angular histogram for the entire past path
if plot_mean_angular_hist == 1
    plot(pcnt_angles_past, angle_bin_centers .* 100, 'b', 'linewidth', 2); 
    hold on;
end

% If the flag for time bin 1 is set, plot the angular histogram for past path in the first time bin
if plot_time_bin_1_angular_hist == 1
    plot(pcnt_angles_past_time_bin_1 .* 100, angle_bin_centers, 'k', 'linewidth', 2); 
    hold on;
end

% If the flag for time bin 2 is set, plot the angular histogram for past path in the second time bin
if plot_time_bin_2_angular_hist == 1
    plot(pcnt_angles_past_time_bin_2 .* 100, angle_bin_centers, 'color', [0.6 0.6 0.6], 'linewidth', 2); 
    hold on
end

% Turn off the box around the plot
box off

% Set the X and Y axis limits for the angular histogram
xlim([0 10])
ylim([min(angle_bins) max(angle_bins)])

% Set the Y-ticks for the angle bins
yticks(angle_bins(1:6:end))

% Label the X-axis as the percentage of events
xlabel('% events')

% Set figure properties for saving
set(gcf, 'Color', 'white', 'Renderer', 'painters', 'PaperPositionMode', 'auto');

% Save the figure as both JPG and PDF in the specified file path
saveas(gcf, fullfile(fig_path, 'time_versus_past_future_angle'), 'jpg')
saveas(gcf, fullfile(fig_path, 'time_versus_past_future_angle'), 'pdf')
