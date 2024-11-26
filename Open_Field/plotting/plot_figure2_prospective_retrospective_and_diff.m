% Set axis limits for the rate plots and rate differences
% These are commented out but may be used for setting axis limits in the plot
% ylim_rates = [0 0.35];
% ylim_rate_diffs = [-0.2 0.3];
bin_centers_plot = bin_centers(bin_centers<=14);

sig_test = 'means';

% Significance level (alpha)
alpha = 0.05;

% Title based on the type of trials (home, away, or all trials)
if home_trials_only == 1
    title_ending = 'home_trials';
elseif away_trials_only == 1
    title_ending = 'away_trials';
else
    title_ending = 'all_trials';
end

% Transparency percentage for colors
transparency_pcnt = 1;

% Adjust colors based on transparency for plotting
colors_2 = [[1 - transparency_pcnt * (1 - colors(1, 1)), 1 - transparency_pcnt * (1 - colors(1, 2)), 1 - transparency_pcnt * (1 - colors(1, 3))]; ...
    [1 - transparency_pcnt * (1 - colors(2, 1)), 1 - transparency_pcnt * (1 - colors(2, 2)), 1 - transparency_pcnt * (1 - colors(2, 3))]];

% Calculate the number of valid (non-NaN) entries for each group (inds1 and inds2)
data1_n = [min(sum(~isnan(combined_hists_rates.future(inds1, bin_centers <= 10)))) max(sum(~isnan(combined_hists_rates.future(inds1, bin_centers <= 10))))];
data2_n = [min(sum(~isnan(combined_hists_rates.future(inds2, bin_centers <= 10)))) max(sum(~isnan(combined_hists_rates.future(inds2, bin_centers <= 10))))];

%% Version 1: Plot Group1 v. Group2 Propsective, Plot Group1 v. Group2 Retrospective, Plot Group1 Difference v. Group2 Difference)

% Set up figure for plotting
% Creating a figure with a specific size and position
figure('Position',[1921 560 350 125])
% Create a 1x3 grid layout for the plots with tight spacing
tiledlayout(1,3,'TileSpacing','tight')

ax1 = nexttile(1);
data = combined_hists_rates.future(inds1, :); % Data for group 1 (home)
x_mean = nanmean(data); % Mean value (ignoring NaNs)
x_err = nanstd(data) ./ sqrt(sum(~isnan(data))); % Standard error of the mean

% Plot the shaded error bar
if align_to_drink_offset == 1
    if smooth_rate_plots == 1
        h = shadedErrorBar(bin_centers_plot((filter_length - 1) / 2:end), x_mean((filter_length - 1) / 2:end), ...
            x_err((filter_length - 1) / 2:end), 'lineprops', 'k'); hold on;
    else
        h = shadedErrorBar(bin_centers_plot, x_mean, x_err, 'lineprops', 'k'); hold on;
    end
else
    h = shadedErrorBar(bin_centers_plot, x_mean, x_err, 'lineprops', 'k'); hold on;
end
% Customize plot appearance (color, line width, etc.)
h.patch.FaceColor = colors(1, :);
h.mainLine.Color = colors(1, :);
h.mainLine.LineWidth = 1;
h.patch.FaceColor = colors_2(1, :);
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';

hold on
% Plot data for group 2 (away)
data = combined_hists_rates.future(inds2, :);
x_mean = nanmean(data);
x_err = nanstd(data) ./ sqrt(sum(~isnan(data)));
if align_to_drink_offset == 1
    if smooth_rate_plots == 1
        h = shadedErrorBar(bin_centers_plot((filter_length - 1) / 2:end), x_mean((filter_length - 1) / 2:end), ...
            x_err((filter_length - 1) / 2:end), 'lineprops', '--k'); hold on;
    else
        h = shadedErrorBar(bin_centers_plot, x_mean, x_err, 'lineprops', '--k'); hold on;
    end
    xlim([-10 0]);
else
    h = shadedErrorBar(bin_centers_plot, x_mean, x_err, 'lineprops', '--k'); hold on;
    xlim([0 10]);
end
% Customize plot appearance for the second group
h.patch.FaceColor = colors(1, :);
h.mainLine.Color = colors(1, :);
h.mainLine.LineWidth = 1;
h.patch.FaceColor = colors_2(1, :);
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';

% Set labels and axis limits
ylabel('Events/s')
ylim(ylim_rates)
yticks(ylim_rates(1):0.1:ylim_rates(2))
xticks(0:2:10)
xtickangle(0);

clear pvals
% Perform statistical significance tests (t-tests or ranksum tests based on user input)
if strcmp(sig_test, 'means')
    for i = 1:size(combined_hists_rates.future, 2)
        a = combined_hists_rates.future(inds1, i);
        b = combined_hists_rates.future(inds2, i);
        [~, pvals(i)] = ttest2(a, b); % Perform two-sample t-test
    end
    sig_times = pvals < alpha;
elseif strcmp(sig_test, 'medians')
    for i = 1:size(combined_hists_rates.future, 2)
        a = combined_hists_rates.future(inds1, i);
        b = combined_hists_rates.future(inds2, i);
        [pvals(i)] = ranksum(a, b); % Perform ranksum test
    end
    sig_times = pvals < alpha;
end

% Mark significant times (with a black dot at the top of the plot)
sig_times = logical(sig_times);
plot(bin_centers_plot(sig_times), ylimit(2) * ones(sum(sig_times),1), '.k')

ax2 = nexttile(2);
data = combined_hists_rates.past(inds1, :); % Data for group 1 (home)
x_mean = nanmean(data);
x_err = nanstd(data) ./ sqrt(sum(~isnan(data)));
if align_to_drink_offset == 1
    if smooth_rate_plots == 1
        h = shadedErrorBar(bin_centers_plot((filter_length - 1) / 2:end), x_mean((filter_length - 1) / 2:end), ...
            x_err((filter_length - 1) / 2:end), 'lineprops', '-k'); hold on;
    else
        h = shadedErrorBar(bin_centers_plot, x_mean, x_err, 'lineprops', '-m'); hold on;
    end
    xlim([-10 0]);
else
    h = shadedErrorBar(bin_centers_plot, x_mean, x_err, 'lineprops', '-m'); hold on;
    xlim([0 10]);
end
h.patch.FaceColor = colors(2, :);
h.mainLine.Color = colors(2, :);
h.mainLine.LineWidth = 1;
h.patch.FaceColor = colors_2(2, :);
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';

hold on
% Plot data for group 2 (away)
data = combined_hists_rates.past(inds2, :);
x_mean = nanmean(data);
x_err = nanstd(data) ./ sqrt(sum(~isnan(data)));
if align_to_drink_offset == 1
    if smooth_rate_plots == 1
        h = shadedErrorBar(bin_centers_plot((filter_length - 1) / 2:end), x_mean((filter_length - 1) / 2:end), ...
            x_err((filter_length - 1) / 2:end), 'lineprops', '--k'); hold on;
    else
        h = shadedErrorBar(bin_centers_plot, x_mean, x_err, 'lineprops', '--m'); hold on;
    end
    xlim([-10 0]);
else
    h = shadedErrorBar(bin_centers_plot, x_mean, x_err, 'lineprops', '--k'); hold on;
    xlim([0 10]);
end
h.patch.FaceColor = colors(2, :);
h.mainLine.Color = colors(2, :);
h.mainLine.LineWidth = 1;
h.patch.FaceColor = colors_2(2, :);
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';

% Set labels and axis limits
ylabel('Events/s')
ylim(ylim_rates)
yticks(ylim_rates(1):0.1:ylim_rates(2))
xticks(0:2:10)
xtickangle(0);
xlabel('Time since stopping (s)')

nexttile(3)
hold on
% Plot rate difference for group 1
data1 = combined_hists_rates.future(inds1,:) - combined_hists_rates.past(inds1,:);
x_mean = nanmean(data1);
x_err = nanstd(data1)./sqrt(sum(~isnan(data1)));

if align_to_drink_offset==1
    if smooth_rate_plots==1
        h = shadedErrorBar(bin_centers_plot((filter_length-1)/2:end),x_mean((filter_length-1)/2:end),...
            x_err((filter_length-1)/2:end),'lineprops','-k'); hold on;
    else
        h = shadedErrorBar(bin_centers_plot,x_mean,...
            x_err,'lineprops','-k'); hold on;
    end
    xlim([-10 0]);
else
    h = shadedErrorBar(bin_centers_plot,x_mean,...
        x_err,'lineprops','-k'); hold on;
    xlim([0 10]);
end
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';
h.mainLine.LineWidth = 1;

% plot rate difference for group 2
data2 = combined_hists_rates.future(inds2,:) - combined_hists_rates.past(inds2,:);
x_mean = nanmean(data2);
x_err = nanstd(data2)./sqrt(sum(~isnan(data2)));
hold on;

if align_to_drink_offset==1
    if smooth_rate_plots==1
        h = shadedErrorBar(bin_centers_plot((filter_length-1)/2:end),x_mean((filter_length-1)/2:end),...
            x_err((filter_length-1)/2:end),'lineprops','--k'); hold on;
    else
        h = shadedErrorBar(bin_centers_plot,x_mean,...
            x_err,'lineprops','--k'); hold on;
    end
    xlim([-10 0]);
else
    h = shadedErrorBar(bin_centers_plot,x_mean,...
        x_err,'lineprops','--k'); hold on;
    xlim([0 10]);
end
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';
h.mainLine.LineWidth = 1;
h.mainLine.LineWidth = 1;

hold on
yline(0)
 box off
ylim(ylim_rate_diffs)
yticks(ylim_rate_diffs(1):0.1:ylim_rate_diffs(2))
xtickangle(0);
xticks(0:2:10)

% Statistics:
% Calculate the absolute difference in means between the two groups (data1 and data2)
diff_of_diffs = abs(nanmean(data1) - nanmean(data2));

% Initialize an array to store p-values (for statistical testing)
pvals = nan(size(combined_hists_rates.future, 2), 1);

% Check which statistical test is selected (shuffle, means, or medians)
if strcmp(sig_test, 'shuffle') == 1
    % If the 'shuffle' test is selected, initialize a matrix to store shuffled differences
    shuffle_diff_of_diffs = nan(1000, size(diff_of_diffs, 2));

    % Perform 1000 shuffle iterations to generate a distribution of differences under the null hypothesis
    for i = 1:1000
        % Shuffle the indices for group 1 (inds1) and group 2 (inds2)
        shuffled_inds1 = randperm(size(combined_hists_rates.future, 1), size(inds1, 1));
        shuffled_inds2 = setdiff(1:size(combined_hists_rates.future, 1), shuffled_inds1);
        
        % Calculate the difference between future and past rates for shuffled group 1 and 2
        shuffle1_diff = combined_hists_rates.future(shuffled_inds1, :) - combined_hists_rates.past(shuffled_inds1, :);
        shuffle2_diff = combined_hists_rates.future(shuffled_inds2, :) - combined_hists_rates.past(shuffled_inds2, :);
        
        % Store the absolute difference of means for the shuffled data
        shuffle_diff_of_diffs(i, :) = abs(nanmean(shuffle1_diff) - nanmean(shuffle2_diff));
    end

    % Compute the significance times based on whether the observed difference exceeds the quantile from the shuffled distribution
    sig_times = diff_of_diffs > quantile(shuffle_diff_of_diffs, 1 - alpha / 2);

elseif strcmp(sig_test, 'means')
    % If the 'means' test is selected, perform a t-test between the differences in future and past rates for both groups
    for i = 1:size(combined_hists_rates.future, 2)
        % Calculate the difference between future and past rates for group 1 and group 2
        a = combined_hists_rates.future(inds1, i) - combined_hists_rates.past(inds1, i);
        b = combined_hists_rates.future(inds2, i) - combined_hists_rates.past(inds2, i);

        % Perform a two-sample t-test between the two groups
        [~, pvals(i)] = ttest2(a, b);
    end

    % Determine the significant times based on the p-values
    sig_times = pvals < alpha;

elseif strcmp(sig_test, 'medians')
    % If the 'medians' test is selected, perform a ranksum test (non-parametric test) between the differences in future and past rates for both groups
    for i = 1:size(combined_hists_rates.future, 2)
        % Calculate the difference between future and past rates for group 1 and group 2
        a = combined_hists_rates.future(inds1, i) - combined_hists_rates.past(inds1, i);
        b = combined_hists_rates.future(inds2, i) - combined_hists_rates.past(inds2, i);

        % Perform a ranksum test between the two groups
        [pvals(i)] = ranksum(a, b);
    end

    % Determine the significant times based on the p-values
    sig_times = pvals < alpha;
end

% Convert the logical vector of significant times to boolean
sig_times = logical(sig_times);

% Plot the results for significant times (highlighted as black dots on the plot)
hold on
ylim(ylim_rate_diffs) % Set the y-axis limits for rate differences
yticks(ylim_rate_diffs(1):0.1:ylim_rate_diffs(2)) % Set y-ticks for the plot
ylimit = gca().YLim; % Get the current y-axis limits
plot(bin_centers_plot(sig_times), ylimit(2) * ones(sum(sig_times), 1), '.k') % Plot significant times as black dots

% Set figure properties (background color, renderer, etc.)
set(gcf, 'Color', 'white', 'Renderer', 'painters', 'PaperPositionMode', 'auto');

% Save the figure as both a JPG and PDF in the specified figure path with a dynamic file name
saveas(gcf, fullfile(fig_path, ['replay_rates_' category1 '_versus_' category2 '_DegThr_' num2str(thrForCategorization_include) ...
    'sdeThr_' num2str(sd_thr) '_disp' num2str(replay_dispersionThr) '_' title_ending]), 'jpg');
saveas(gcf, fullfile(fig_path, ['replay_rates_' category1 '_versus_' category2 '_DegThr_' num2str(thrForCategorization_include) ...
    'sdeThr_' num2str(sd_thr) '_disp' num2str(replay_dispersionThr) '_' title_ending]), 'pdf');


%% Version 2: (Plot Group1 Propsective v Retrospective, Plot Group2 Prospective v Retrospective, Plot Group1 Difference v. Group2 Difference)

% Create a figure window with a specific size
figure('Position', [1921 560 350 125])

% Set up a tiled layout with 1 row and 3 columns for plotting
tiledlayout(1, 3, 'TileSpacing', 'tight')


% Plot on the first tile of the layout
ax1 = nexttile(1);

% Extract and process the data for home prospective (group 1)
data = combined_hists_rates.future(inds1, :);
x_mean = nanmean(data);  % Compute the mean of the data (ignoring NaNs)
x_err = nanstd(data) ./ sqrt(sum(~isnan(data)));  % Compute the standard error

% Conditional alignment and smoothing logic
if align_to_drink_offset == 1
    if smooth_rate_plots == 1
        h = shadedErrorBar(bin_centers_plot((filter_length - 1) / 2:end), x_mean((filter_length - 1) / 2:end), ...
            x_err((filter_length - 1) / 2:end), 'lineprops', '-k'); hold on;
    else
        h = shadedErrorBar(bin_centers_plot, x_mean, x_err, 'lineprops', '-m'); hold on;
    end
    xlim([-10 0]);  % Set x-axis limits if aligning to drink offset
else
    h = shadedErrorBar(bin_centers_plot, x_mean, x_err, 'lineprops', '-k'); hold on;
    xlim([0 10]);  % Set x-axis limits if not aligning to drink offset
end

% Customize the appearance of the plot
h.patch.FaceColor = colors(1, :);
h.mainLine.Color = colors(1, :);
h.mainLine.LineWidth = 1;
h.patch.FaceColor = colors_2(1, :);
h.edge(1).Color = 'none';  % Hide edges of the shaded area
h.edge(2).Color = 'none';  % Hide edges of the shaded area

% Plot past data for home prospective (group 1)
hold on
data = combined_hists_rates.past(inds1, :);
x_mean = nanmean(data);
x_err = nanstd(data) ./ sqrt(sum(~isnan(data)));

if align_to_drink_offset == 1
    if smooth_rate_plots == 1
        h = shadedErrorBar(bin_centers_plot((filter_length - 1) / 2:end), x_mean((filter_length - 1) / 2:end), ...
            x_err((filter_length - 1) / 2:end), 'lineprops', '-k'); hold on;
    else
        h = shadedErrorBar(bin_centers_plot, x_mean, x_err, 'lineprops', '-k'); hold on;
    end
    xlim([-10 0]);
else
    h = shadedErrorBar(bin_centers_plot, x_mean, x_err, 'lineprops', '-k'); hold on;
    xlim([0 10]);
end

% Customize the appearance of the past data plot
h.patch.FaceColor = colors(2, :);
h.mainLine.Color = colors(2, :);
h.mainLine.LineWidth = 1;
h.patch.FaceColor = colors_2(2, :);
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';

% Label the y-axis and set axis properties
ylabel('Events/s')
ylim(ylim_rates)  % Set y-axis limits
yticks(ylim_rates(1):0.1:ylim_rates(2))  % Set y-ticks
xtickangle(0);
xticks(0:2:10)

% Perform statistical significance testing based on the selected method
if strcmp(sig_test, 'means')
    for i = 1:size(combined_hists_rates.future, 2)
        a = combined_hists_rates.future(inds1, i);
        b = combined_hists_rates.past(inds1, i);
        [~, pvals(i)] = ttest2(a, b);  % Perform t-test for means
    end
    sig_times = pvals < alpha;  % Identify significant times based on p-values
elseif strcmp(sig_test, 'medians')
    for i = 1:size(combined_hists_rates.future, 2)
        a = combined_hists_rates.future(inds1, i);
        b = combined_hists_rates.past(inds1, i);
        [pvals(i)] = ranksum(a, b);  % Perform ranksum test for medians
    end
    sig_times = pvals < alpha;  % Identify significant times based on p-values
end

% Plot significant times as black dots
plot(bin_centers_plot(sig_times), ylimit(2) * ones(sum(sig_times), 1), '.k')

% Plot on the second tile of the layout
ax2 = nexttile(2);

% Extract and process data for away prospective (group 2)
data = combined_hists_rates.future(inds2, :);
x_mean = nanmean(data);
x_err = nanstd(data) ./ sqrt(sum(~isnan(data)));

if align_to_drink_offset == 1
    if smooth_rate_plots == 1
        h = shadedErrorBar(bin_centers_plot((filter_length - 1) / 2:end), x_mean((filter_length - 1) / 2:end), ...
            x_err((filter_length - 1) / 2:end), 'lineprops', '--k'); hold on;
    else
        h = shadedErrorBar(bin_centers_plot, x_mean, x_err, 'lineprops', '--k'); hold on;
    end
    xlim([-10 0]);
else
    h = shadedErrorBar(bin_centers_plot, x_mean, x_err, 'lineprops', '--k'); hold on;
    xlim([0 10]);
end

% Customize the appearance of the plot
h.patch.FaceColor = colors(1, :);
h.mainLine.Color = colors(1, :);
h.mainLine.LineWidth = 1;
h.patch.FaceColor = colors_2(1, :);
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';

% Plot past data for away prospective (group 2)
hold on
data = combined_hists_rates.past(inds2, :);
x_mean = nanmean(data);
x_err = nanstd(data) ./ sqrt(sum(~isnan(data)));

if align_to_drink_offset == 1
    if smooth_rate_plots == 1
        h = shadedErrorBar(bin_centers_plot((filter_length - 1) / 2:end), x_mean((filter_length - 1) / 2:end), ...
            x_err((filter_length - 1) / 2:end), 'lineprops', '--k'); hold on;
    else
        h = shadedErrorBar(bin_centers_plot, x_mean, x_err, 'lineprops', '--k'); hold on;
    end
    xlim([-10 0]);
else
    h = shadedErrorBar(bin_centers_plot, x_mean, x_err, 'lineprops', '--k'); hold on;
    xlim([0 10]);
end

% Customize the appearance of the plot
h.patch.FaceColor = colors(2, :);
h.mainLine.Color = colors(2, :);
h.mainLine.LineWidth = 1;
h.patch.FaceColor = colors_2(2, :);
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';

% Label the x-axis and set axis properties
xlabel('Time since arrival (s)')
ylim(ylim_rates)
yticks(ylim_rates(1):0.1:ylim_rates(2))
xtickangle(0);
xticks(0:2:10)

% Perform statistical significance testing based on the selected method
if strcmp(sig_test, 'means')
    for i = 1:size(combined_hists_rates.future, 2)
        a = combined_hists_rates.future(inds2, i);
        b = combined_hists_rates.past(inds2, i);
        [~, pvals(i)] = ttest2(a, b);
    end
    sig_times = pvals < alpha;
elseif strcmp(sig_test, 'medians')
    for i = 1:size(combined_hists_rates.future, 2)
        a = combined_hists_rates.future(inds2, i);
        b = combined_hists_rates.past(inds2, i);
        [pvals(i)] = ranksum(a, b);
    end
    sig_times = pvals < alpha;
end

% Plot significant times as black dots
plot(bin_centers_plot(sig_times), ylimit(2) * ones(sum(sig_times), 1), '.k')

nexttile(3)
hold on
% Calculate the difference of future vs. past for both groups
data1 = combined_hists_rates.future(inds1, :) - combined_hists_rates.past(inds1, :);
x_mean = nanmean(data1);
x_err = nanstd(data1) ./ sqrt(sum(~isnan(data1)));

if align_to_drink_offset == 1
    if smooth_rate_plots == 1
        h = shadedErrorBar(bin_centers_plot((filter_length - 1) / 2:end), x_mean((filter_length - 1) / 2:end), ...
            x_err((filter_length - 1) / 2:end), 'lineprops', '-k'); hold on;
    else
        h = shadedErrorBar(bin_centers_plot, x_mean, x_err, 'lineprops', '-k'); hold on;
    end
    xlim([-10 0]);
else
    h = shadedErrorBar(bin_centers_plot, x_mean, x_err, 'lineprops', '-k'); hold on;
    xlim([0 10]);
end

h.edge(1).Color = 'none';
h.edge(2).Color = 'none';
h.mainLine.LineWidth = 1;

% Calculate the difference for the second group (away prospective)
data2 = combined_hists_rates.future(inds2, :) - combined_hists_rates.past(inds2, :);
x_mean = nanmean(data2);
x_err = nanstd(data2) ./ sqrt(sum(~isnan(data2)));

hold on
if align_to_drink_offset == 1
    if smooth_rate_plots == 1
        h = shadedErrorBar(bin_centers_plot((filter_length - 1) / 2:end), x_mean((filter_length - 1) / 2:end), ...
            x_err((filter_length - 1) / 2:end), 'lineprops', '--k'); hold on;
    else
        h = shadedErrorBar(bin_centers_plot, x_mean, x_err, 'lineprops', '--k'); hold on;
    end
    xlim([-10 0]);
else
    h = shadedErrorBar(bin_centers_plot, x_mean, x_err, 'lineprops', '--k'); hold on;
    xlim([0 10]);
end

h.edge(1).Color = 'none';
h.edge(2).Color = 'none';
h.mainLine.LineWidth = 1;
h.mainLine.LineWidth = 1;

hold on
yline(0)  % Add a horizontal line at zero

% Customize axis properties
box off
ylim(ylim_rate_diffs)  % Set y-axis limits for differences
yticks(ylim_rate_diffs(1):0.1:ylim_rate_diffs(2))
xtickangle(0);
xticks(0:2:10)

% Perform statistical significance testing for the differences
pvals = nan(size(combined_hists_rates.future, 2), 1);
if strcmp(sig_test, 'means')
    for i = 1:size(combined_hists_rates.future, 2)
        a = combined_hists_rates.future(inds1, i) - combined_hists_rates.past(inds1, i);
        b = combined_hists_rates.future(inds2, i) - combined_hists_rates.past(inds2, i);
        [~, pvals(i)] = ttest2(a, b);
    end
    sig_times = pvals < alpha;
elseif strcmp(sig_test, 'medians')
    for i = 1:size(combined_hists_rates.future, 2)
        a = combined_hists_rates.future(inds1, i) - combined_hists_rates.past(inds1, i);
        b = combined_hists_rates.future(inds2, i) - combined_hists_rates.past(inds2, i);
        [pvals(i)] = ranksum(a, b);
    end
    sig_times = pvals < alpha;
end

% Plot significant times for the difference of differences as black dots
hold on
ylim(ylim_rate_diffs)
yticks(ylim_rate_diffs(1):0.1:ylim_rate_diffs(2))
ylimit = gca().YLim;
plot(bin_centers_plot(sig_times), ylimit(2) * ones(sum(sig_times), 1), '.k')

% Save the figure as both a .jpg and a .pdf file
set(gcf, 'Color', 'white', 'Renderer', 'painters', 'PaperPositionMode', 'auto');
saveas(gcf, fullfile(fig_path, ['replay_rates_v2_' category1 '_versus_' category2 '_DegThr_' num2str(thrForCategorization_include) 'sdeThr_' num2str(sd_thr) '_disp' num2str(replay_dispersionThr) '_' title_ending]), 'jpg')
saveas(gcf, fullfile(fig_path, ['replay_rates_v2_' category1 '_versus_' category2 '_DegThr_' num2str(thrForCategorization_include) 'sdeThr_' num2str(sd_thr) '_disp' num2str(replay_dispersionThr) '_' title_ending]), 'pdf')


%% Version 3: (Plots Prospective and Retrospective Rates for Group1 and Group2 on top of each other)
% Creating a figure with a specific size and position
figure('Position',[1921 560 350 125])
% Create a 1x3 grid layout for the plots with tight spacing
tiledlayout(1,3,'TileSpacing','tight')

% Plot prospective, group1
ax1 = nexttile(1);
data = combined_hists_rates.future(inds1,:);
x_mean = nanmean(data);
x_err = nanstd(data)./sqrt(sum(~isnan(data)));
if align_to_drink_offset==1
    if smooth_rate_plots==1
h = shadedErrorBar(bin_centers_plot((filter_length-1)/2:end),x_mean((filter_length-1)/2:end),...
    x_err((filter_length-1)/2:end),'lineprops','-k'); hold on;
    else
         h = shadedErrorBar(bin_centers_plot,x_mean,...
    x_err,'lineprops','-k'); hold on;
    end
    xlim([-10 0]);
else
 h = shadedErrorBar(bin_centers_plot,x_mean,...
    x_err,'lineprops','-k'); hold on;
    xlim([0 10]);
end
% Set line and patch colors
h.patch.FaceColor = colors(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.patch.FaceColor = colors_2(1,:);
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';
hold on

% Plot retrospective, group1
data = combined_hists_rates.past(inds1,:);
x_mean = nanmean(data);
x_err = nanstd(data)./sqrt(sum(~isnan(data)));
if align_to_drink_offset==1
    if smooth_rate_plots==1
h = shadedErrorBar(bin_centers_plot((filter_length-1)/2:end),x_mean((filter_length-1)/2:end),...
    x_err((filter_length-1)/2:end),'lineprops','-k'); hold on;
    else
         h = shadedErrorBar(bin_centers_plot,x_mean,...
    x_err,'lineprops','-k'); hold on;
    end
    xlim([-10 0]);
else
 h = shadedErrorBar(bin_centers_plot,x_mean,...
    x_err,'lineprops','-k'); hold on;
    xlim([0 10]);
end
% Set line and patch colors
h.patch.FaceColor = colors(2,:);
h.mainLine.Color = colors(2,:);
h.mainLine.LineWidth = 1;
h.patch.FaceColor = colors_2(2,:);
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';

% Plot prospective, group2
hold on
data = combined_hists_rates.future(inds2,:);
x_mean = nanmean(data);
x_err = nanstd(data)./sqrt(sum(~isnan(data)));
if align_to_drink_offset==1
    if smooth_rate_plots==1
h = shadedErrorBar(bin_centers_plot((filter_length-1)/2:end),x_mean((filter_length-1)/2:end),...
    x_err((filter_length-1)/2:end),'lineprops','--k'); hold on;
    else
         h = shadedErrorBar(bin_centers_plot,x_mean,...
    x_err,'lineprops','--k'); hold on;
    end
    xlim([-10 0]);
else
 h = shadedErrorBar(bin_centers_plot,x_mean,...
    x_err,'lineprops','--k'); hold on;
    xlim([0 10]);
end
% Set line and patch colors
h.patch.FaceColor = colors(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.patch.FaceColor = colors_2(1,:);
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';
hold on

% Plot retrospective, group2
data = combined_hists_rates.past(inds2,:);
x_mean = nanmean(data);
x_err = nanstd(data)./sqrt(sum(~isnan(data)));
if align_to_drink_offset==1
    if smooth_rate_plots==1
h = shadedErrorBar(bin_centers_plot((filter_length-1)/2:end),x_mean((filter_length-1)/2:end),...
    x_err((filter_length-1)/2:end),'lineprops','--k'); hold on;
    else
         h = shadedErrorBar(bin_centers_plot,x_mean,...
    x_err,'lineprops','--k'); hold on;
    end
    xlim([-10 0]);
else
 h = shadedErrorBar(bin_centers_plot,x_mean,...
    x_err,'lineprops','--k'); hold on;
    xlim([0 10]);
end
% Set line and patch colors
h.patch.FaceColor = colors(2,:);
h.mainLine.Color = colors(2,:);
h.mainLine.LineWidth = 1;
h.patch.FaceColor = colors_2(2,:);
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';

%legend({'Future, off', 'Past, off', 'Future, on', 'Past, on'})
ylim(ylim_rates)
yticks(ylim_rates(1):0.1:ylim_rates(2))
xtickangle(0);
xticks(0:2:10)
ylabel('Events/s')


% Second plot: Copy of the first plot in the second tile, just for
% formating purposes.
ax2 = nexttile(2);
ax1Chil = ax1.Children;
copyobj(ax1Chil,ax2);
ylim(ylim_rates)
yticks(ylim_rates(1):0.1:ylim_rates(2))
xtickangle(0);
xticks(0:2:10)
xlabel('Time since arrival (s)')
if align_to_drink_offset
    xlim([-10 0]);
else
    xlim([0 10]);
end

nexttile(3)
% Third plot: Difference of differences between future and past for both
% groups

% Differences group2
data1 = combined_hists_rates.future(inds1,:) - combined_hists_rates.past(inds1,:);
x_mean = nanmean(data1);
x_err = nanstd(data1)./sqrt(sum(~isnan(data1)));
if align_to_drink_offset==1
    if smooth_rate_plots==1
h = shadedErrorBar(bin_centers_plot((filter_length-1)/2:end),x_mean((filter_length-1)/2:end),...
    x_err((filter_length-1)/2:end),'lineprops','-k'); hold on;
    else
         h = shadedErrorBar(bin_centers_plot,x_mean,...
    x_err,'lineprops','-k'); hold on;
    end
    xlim([-10 0]);
else
 h = shadedErrorBar(bin_centers_plot,x_mean,...
    x_err,'lineprops','-k'); hold on;
    xlim([0 10]);
end
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';
h.mainLine.LineWidth = 1;


% Differences group2
data2 = combined_hists_rates.future(inds2,:) - combined_hists_rates.past(inds2,:);
x_mean = nanmean(data2);
x_err = nanstd(data2)./sqrt(sum(~isnan(data2)));
hold on;
if align_to_drink_offset==1
    if smooth_rate_plots==1
h = shadedErrorBar(bin_centers_plot((filter_length-1)/2:end),x_mean((filter_length-1)/2:end),...
    x_err((filter_length-1)/2:end),'lineprops','--k'); hold on;
    else
         h = shadedErrorBar(bin_centers_plot,x_mean,...
    x_err,'lineprops','--k'); hold on;
    end
    xlim([-10 0]);
else
 h = shadedErrorBar(bin_centers_plot,x_mean,...
    x_err,'lineprops','--k'); hold on;
    xlim([0 10]);
end
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';
h.mainLine.LineWidth = 1;
h.mainLine.LineWidth = 1;

hold on
yline(0)
box off
ylim([-0.2 0.2])
yticks(ylim_rates(1):0.1:ylim_rates(2))
xtickangle(0);
xticks(0:2:10)

% Statistical testing and significance marking (means, or medians)
% Apply shuffle or statistical tests based on the chosen method
pvals = nan(size(combined_hists_rates.future,2),1);
if strcmp(sig_test,'means')
    for i = 1:size(combined_hists_rates.future,2)
        a = combined_hists_rates.future(inds1,i)-combined_hists_rates.past(inds1,i);
        b = combined_hists_rates.future(inds2,i)-combined_hists_rates.past(inds2,i);
        [~,pvals(i)] = ttest2(a,b);
    end
    sig_times = pvals < alpha;
elseif strcmp(sig_test,'medians')
    for i = 1:size(combined_hists_rates.future,2)
        a = combined_hists_rates.future(inds1,i)-combined_hists_rates.past(inds1,i);
        b = combined_hists_rates.future(inds2,i)-combined_hists_rates.past(inds2,i);
        [pvals(i)] = ranksum(a,b);
    end
    sig_times = pvals < alpha;
end
sig_times = logical(sig_times);

hold on
ylim(ylim_rate_diffs)
yticks(ylim_rate_diffs(1):0.1:ylim_rate_diffs(2))
ylimit = gca().YLim;
plot(bin_centers_plot(sig_times),ylimit(2)*ones(sum(sig_times),1),'.k')

set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
saveas(gcf,fullfile(fig_path,['replay_rates_v3_' category1 '_versus_' category2 '_DegThr_' num2str(thrForCategorization_include) 'sdeThr_' num2str(sd_thr) '_disp' num2str(replay_dispersionThr) '_' title_ending]),'jpg')
saveas(gcf,fullfile(fig_path,['replay_rates_v3_' category1 '_versus_' category2 '_DegThr_' num2str(thrForCategorization_include) 'sdeThr_' num2str(sd_thr) '_disp' num2str(replay_dispersionThr) '_' title_ending]),'pdf')
