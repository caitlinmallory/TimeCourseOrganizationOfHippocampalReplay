%% Look at the overall rate of retrospective or prospective replay on each trial type.

% Define the bin range for analysis
bins = 1:20;

% Calculate the rate of prospective replay for group 1 (i.e., future-related behavior)
rate_prospective_inds1 = nansum(combined_hists.future(inds1,bins),2)./(min(10,C.drink_period_time(inds1,:)));

% Calculate the rate of retrospective replay for group 1 (i.e., past-related behavior)
rate_retrospective_inds1 = nansum(combined_hists.past(inds1,bins),2)./(min(10,C.drink_period_time(inds1)));

% Calculate the rate of prospective replay for group 2 (i.e., future-related behavior)
rate_prospective_inds2 = nansum(combined_hists.future(inds2,bins),2)./(min(10,C.drink_period_time(inds2,:)));

% Calculate the rate of retrospective replay for group 2 (i.e., past-related behavior)
rate_retrospective_inds2 = nansum(combined_hists.past(inds2,bins),2)./(min(10,C.drink_period_time(inds2,:)));

% Perform a statistical comparison using the ranksum (Mann-Whitney U test) for the future replay rate between group 1 and group 2
[pval_future_inds1_inds2, ~, zval_future_inds1_inds2] = ranksum(rate_prospective_inds1, rate_prospective_inds2);

% Perform a statistical comparison for the past replay rate between group 1 and group 2
[pval_past_inds1_inds2, h, zval_past_inds1_inds2] = ranksum(rate_retrospective_inds1, rate_retrospective_inds2);

%% Shuffled p-values and shuffled interaction:

% Number of shuffles to run for generating null distributions
num_shuffles = 10000;

% Calculate the real difference in means for future and past replay rates
real_diff_future = mean(rate_prospective_inds1) - mean(rate_prospective_inds2);
real_diff_past = mean(rate_retrospective_inds1) - mean(rate_retrospective_inds2);

% Calculate the real interaction effect (difference between past and future differences)
real_interaction = real_diff_past - real_diff_future;

% Initialize arrays to store shuffled differences
shuffled_diff_future = nan(num_shuffles, 1);
shuffled_diff_past = nan(num_shuffles, 1);

% Perform shuffling to calculate null distributions for future and past replay rates
for n = 1:num_shuffles
    % Shuffle the indices for group 1
    shuffled_inds1 = randperm(height(combined_hists.future), length(inds1));
    % Get the indices for group 2 as the complement of shuffled_inds1
    shuffled_inds2 = setdiff(1:height(combined_hists.future), shuffled_inds1);

    % Calculate the difference in future replay rates for the shuffled data
    shuffled_diff_future(n) = mean(nansum(combined_hists.future(shuffled_inds1, bins), 2) ./ (min(10, C.drink_period_time(shuffled_inds1, :)))) - ...
        mean(nansum(combined_hists.future(shuffled_inds2, bins), 2) ./ (min(10, C.drink_period_time(shuffled_inds2, :))));

    % Calculate the difference in past replay rates for the shuffled data
    shuffled_diff_past(n) = mean(nansum(combined_hists.past(shuffled_inds1, bins), 2) ./ (min(10, C.drink_period_time(shuffled_inds1, :)))) - ...
        mean(nansum(combined_hists.past(shuffled_inds2, bins), 2) ./ (min(10, C.drink_period_time(shuffled_inds2, :))));
end

% Calculate the interaction effect for the shuffled data
shuffled_interaction = shuffled_diff_past - shuffled_diff_future;

% Calculate the significance of the real future, past, and interaction effects
sig_future = (sum(shuffled_diff_future > abs(real_diff_future)) + sum(shuffled_diff_future < -1 * abs(real_diff_future)) + 1) / (num_shuffles + 1);
sig_past = (sum(shuffled_diff_past > abs(real_diff_past)) + sum(shuffled_diff_past < -1 * abs(real_diff_past)) + 1) / (num_shuffles + 1);
sig_interaction = (sum(shuffled_interaction > abs(real_interaction)) + sum(shuffled_interaction < -1 * abs(real_interaction)) + 1) / (num_shuffles + 1);

% Create a table to store the rates and group information for later analysis
data_table_sub = table();
data_table_sub.rates = [rate_prospective_inds1; rate_retrospective_inds1; rate_prospective_inds2; rate_retrospective_inds2];
data_table_sub.group = [zeros(size(rate_prospective_inds1)); zeros(size(rate_retrospective_inds1)); ones(size(rate_prospective_inds2)); ones(size(rate_retrospective_inds2))];
data_table_sub.replay_class = [ones(size(rate_prospective_inds1)); 2 * ones(size(rate_retrospective_inds1)); ones(size(rate_prospective_inds2)); 2 * ones(size(rate_retrospective_inds2))];

% Convert groups and replay classes into categorical variables
groupOrder = ["1", "2"];
namedGroup = categorical(data_table_sub.group, 0:1, groupOrder);
directionOrder = ["Reverse", "Forward"];
namedDirection = categorical(data_table_sub.replay_class, [2, 1], directionOrder);

% Data for statistical tests (ANOVA)
data_sub = data_table_sub.rates;

%% (Optional) ANOVA to look at the interaction effects between group and replay type
% The ANOVA code is commented out, but it calculates the effects of group and replay class
% [p, tbl, stats] = anovan(data_sub', {namedGroup, namedDirection}, 'model', 'interaction', 'varnames', {'namedGroup', 'namedDirection'});
% [results, ~, ~, gnames] = multcompare(stats, "Dimension", [1 2]);

% Define colors for plotting
colors = [.4660 0.6740 0.1880; 0.4940 0.1840 0.5560;];
transparency_pcnt = 0.5;
colors_2 = [[1 - transparency_pcnt * (1 - colors(1, 1)), 1 - transparency_pcnt * (1 - colors(1, 2)), 1 - transparency_pcnt * (1 - colors(1, 3))]; ...
    [1 - transparency_pcnt * (1 - colors(2, 1)), 1 - transparency_pcnt * (1 - colors(2, 2)), 1 - transparency_pcnt * (1 - colors(2, 3))]];

% Prepare data and error values for plotting
data = [nanmean(rate_prospective_inds1), nanmean(rate_prospective_inds2), nanmean(rate_retrospective_inds1), nanmean(rate_retrospective_inds2)];
err = [nanstd(rate_prospective_inds1) / sqrt(nansum(~isnan(rate_prospective_inds1))), ...
    nanstd(rate_prospective_inds2) / sqrt(nansum(~isnan(rate_prospective_inds2))), ...
    nanstd(rate_retrospective_inds1) / sqrt(nansum(~isnan(rate_retrospective_inds1))), ...
    nanstd(rate_retrospective_inds2) / sqrt(nansum(~isnan(rate_retrospective_inds2)))];

% Create a bar plot to visualize the replays with error bars
figure('Position', [1986 1051 100 100]);
b1 = bar(1, data(1)); hold on;
e1 = errorbar(1, data(1), err(1), 'k', 'linestyle', 'none');
e1.CapSize = 4;
hold on;
b2 = bar(2, data(2)); hold on;
e2 = errorbar(2, data(2), err(2), 'k', 'linestyle', 'none');
e2.CapSize = 4;
b3 = bar(4, data(3)); hold on;
e3 = errorbar(4, data(3), err(3), 'k', 'linestyle', 'none');
e3.CapSize = 4;
b4 = bar(5, data(4)); hold on;
e4 = errorbar(5, data(4), err(4), 'k', 'linestyle', 'none');
e4.CapSize = 4;

% Set the colors for the bars
b1.FaceColor = colors(1, :);
b1.EdgeColor = 'none';
b2.FaceColor = colors_2(1, :);
b2.EdgeColor = 'none';
b3.FaceColor = colors(2, :);
b3.EdgeColor = 'none';
b4.FaceColor = colors_2(2, :);
b4.EdgeColor = 'none';

% Clean up plot appearance
box off;
xticks([1.5 4.5]);
xticklabels({});
ylabel('Events/s');
ylim([0 0.3])

% Set figure properties for saving
set(gcf, 'Color', 'white', 'Renderer', 'painters', 'PaperPositionMode', 'auto');
fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures';

if home_trials_only==1
    title_ending = 'home_trials';
elseif away_trials_only==1
    title_ending = 'away_trials';
else
    title_ending = 'all_trials';
end
saveas(gcf,fullfile(fig_path,['average_stopping_period_rate_' category1 '_' category2 'trials_' title_ending]),'jpg')
saveas(gcf,fullfile(fig_path,['average_stopping_period_rate_' category1 '_' category2 'trials_' title_ending]),'pdf')
