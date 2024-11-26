% Set the path for saving figures
fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures';

% Initialization for figure 1: Home and away conditions
home_only = 0; % Flag to plot only home condition
away_only = 0; % Flag to plot only away condition

% Set color schemes for plots
colors = [.4660 0.6740 0.1880; 0.4940 0.1840 0.5560;];

% Define time window for drink-related events
drink_times = [0 10]; % From 0 to 10 units of time
fig_ylim = [60 100]; % Set y-axis limits for final plot

% Split data based on whether it's home or away condition
group1 = 'home'; % Label for home group
group2 = 'away'; % Label for away group
replay_group1 = replay(replay.home_event == 1,:); % Home event group
replay_group2 = replay(replay.home_event == 0,:); % Away event group


% For Fig X?
% group1 = 'off';
% group2 = 'on';
% replay_group1 = replay(replay.laser_state_binary==0,:);
% replay_group2 = replay(replay.laser_state_binary==1,:);
% if home_only==1
% replay_group1 = replay(replay.laser_state_binary==0 & replay.home_event==1,:);
% replay_group2 = replay(replay.laser_state_binary==1 & replay.home_event==1,:);
% end
% if away_only==1
% replay_group1 = replay(replay.laser_state_binary==0 & replay.home_event==0,:);
% replay_group2 = replay(replay.laser_state_binary==1 & replay.home_event==0,:);
% end

%%
figure('Position',[803 408 400 150])
tiledlayout(1,3)
plot_ind_animals = 1; % Flag for individual animal plotting
if plot_ind_animals == 1
    % Loop through each rat to plot individual data
    for rat = 1:length(rats)
        replay_rat = replay(replay.rat_label == rats(rat),:); % Filter by rat ID
        angle_bins = linspace(0, 180, 19); % Define angle bins for path displacement
        angle_bin_centers = mean([angle_bins(1:end-1)' angle_bins(2:end)'],2); % Calculate bin centers
        replay_past_path_angle_map = nan(length(angle_bin_centers),1); % Placeholder for past angle map
        replay_future_path_angle_map = nan(length(angle_bin_centers),1); % Placeholder for future angle map

        % Compute the counts for past and future angle distributions
        for j = 1:length(angle_bins)-1
            replay_past_path_angle_map(j) = sum(replay_rat.meanAngDisplacement_pastPath >= angle_bins(j) & replay_rat.meanAngDisplacement_pastPath < angle_bins(j+1));
            replay_future_path_angle_map(j) = sum(replay_rat.meanAngDisplacement_futPath >= angle_bins(j) & replay_rat.meanAngDisplacement_futPath < angle_bins(j+1));
        end

        % Plot past and future path angle distributions for the current rat
        ax1 = nexttile(1); % Next plot in the grid
        plot(angle_bin_centers, replay_future_path_angle_map ./ sum(replay_future_path_angle_map) * 100, 'LineWidth', 1); % Future path angle
        xlim([0 180]); ylim([0 16]); hold on;

        ax2 = nexttile(2); % Next plot in the grid
        plot(angle_bin_centers, replay_past_path_angle_map ./ sum(replay_past_path_angle_map) * 100, 'LineWidth', 1); % Past path angle
        xlim([0 180]); ylim([0 16]); hold on;
    end
    hold on
end

% Recalculate the angle maps for all rats combined (summary statistics)
angle_bins = linspace(0, 180, 19);
angle_bin_centers = mean([angle_bins(1:end-1)' angle_bins(2:end)'],2);
replay_past_path_angle_map = nan(length(angle_bin_centers),1);
replay_future_path_angle_map = nan(length(angle_bin_centers),1);

% Plot the summary statistics for future and past path angle distributions
for j = 1:length(angle_bins)-1
    replay_past_path_angle_map(j) = sum(replay.meanAngDisplacement_pastPath>=angle_bins(j) & replay.meanAngDisplacement_pastPath<angle_bins(j+1));
    replay_future_path_angle_map(j) = sum(replay.meanAngDisplacement_futPath>=angle_bins(j) & replay.meanAngDisplacement_futPath<angle_bins(j+1));
end

nexttile(1)
plot(angle_bin_centers,replay_future_path_angle_map./sum(replay_future_path_angle_map).*100,'color',colors(1,:),'LineWidth',3)
xlim([0 180])
xticks([0 90 180])
xtickangle(0)
box off
ylabel('% Events')
xlabel('|Displacement|')

nexttile(2)
plot(angle_bin_centers,replay_past_path_angle_map./sum(replay_past_path_angle_map).*100,'color',colors(2,:),'LineWidth',3)
xlim([0 180])
xticks([0 90 180])
xtickangle(0)
xlabel('|Displacement|')
box off

% Calculate the mean and standard error of the mean for past and future angles
mean_abs_angle_past_rat = nan(length(rats),1);
mean_abs_angle_future_rat = nan(length(rats),1);
for rat = 1:length(rats)
    mean_abs_angle_past_rat(rat) = nanmean(replay.meanAngDisplacement_pastPath(replay.rat_label==rats(rat)));
    mean_abs_angle_future_rat(rat) = nanmean(replay.meanAngDisplacement_futPath(replay.rat_label==rats(rat)));
end

% Compute grand means
grand_mean_abs_angle_past = mean(mean_abs_angle_past_rat);
grand_mean_abs_angle_future = mean(mean_abs_angle_future_rat);

% Perform paired t-test between future and past path angles
[h,p] = ttest( ...
    mean_abs_angle_future_rat,mean_abs_angle_past_rat)

% Calculate means and standard errors for all subjects
mean_abs_angle_past = nanmean(replay.meanAngDisplacement_pastPath);
mean_abs_angle_future = nanmean(replay.meanAngDisplacement_futPath);
sem_abs_angle_past = nanstd(replay.meanAngDisplacement_pastPath)/sqrt(sum(~isnan(replay.meanAngDisplacement_pastPath)));
sem_abs_angle_future = nanstd(replay.meanAngDisplacement_futPath)/sqrt(sum(~isnan(replay.meanAngDisplacement_futPath)));

% Plot bar chart with error bars for the grand means
ax3 = nexttile(3);
[p,h,z] = ranksum(replay.meanAngDisplacement_futPath,replay.meanAngDisplacement_pastPath);
[p,h,z] = signrank(replay.meanAngDisplacement_futPath,replay.meanAngDisplacement_pastPath);

RF_n = sum(~isnan(replay.meanAngDisplacement_futPath));
RP_n = sum(~isnan(replay.meanAngDisplacement_pastPath));

% Bar chart showing the mean displacement for future and past events
x = [1 2];
data = [mean_abs_angle_future mean_abs_angle_past];
err = [sem_abs_angle_future sem_abs_angle_past];
b1 = bar(1,data(1)); hold on;
e1 = errorbar(1,data(1),err(1),'k','linestyle','none');
e1.CapSize = 4;
hold on
b2 = bar(2,data(2)); hold on;
e2 = errorbar(2,data(2),err(2),'k','linestyle','none');
e2.CapSize = 4;
b1.FaceColor = colors(1,:);
b1.EdgeColor = 'none';
b2.FaceColor = colors(2,:);
b2.EdgeColor = 'none';
ylim(fig_ylim)
xlim([0 3])
xticks([1 2])
xticklabels({'RF','RP'})
ylabel('|Displacement|')
box off
hold on

set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
saveas(gcf,fullfile(fig_path,'summary_angles_without_time'),'jpg')
saveas(gcf,fullfile(fig_path,'summary_angles_without_time'),'pdf')
%% Plot group types separately (home versus away, or laser off versus laser on)
figure('Position',[803 408 400 150])
tiledlayout(1,3)
angle_bins = linspace(0, 180, 19);
angle_bin_centers = mean([angle_bins(1:end-1)' angle_bins(2:end)'],2);
replay_past_path_angle_map = nan(length(angle_bin_centers),1);
replay_future_path_angle_map = nan(length(angle_bin_centers),1);

for j = 1:length(angle_bins)-1
    replay_past_path_angle_map(j) = sum(replay_group1.meanAngDisplacement_pastPath>=angle_bins(j) & replay_group1.meanAngDisplacement_pastPath<angle_bins(j+1));
    replay_future_path_angle_map(j) = sum(replay_group1.meanAngDisplacement_futPath>=angle_bins(j) & replay_group1.meanAngDisplacement_futPath<angle_bins(j+1));
end
nexttile(1)
plot(angle_bin_centers,replay_future_path_angle_map./sum(replay_future_path_angle_map).*100,'color',colors(1,:),'LineWidth',2)
xlim([0 180])
xticks([0 90 180])
xtickangle(0)
box off
ylabel('% Events')
xlabel('|Displacement|')

hold on;
nexttile(2)
plot(angle_bin_centers,replay_past_path_angle_map./sum(replay_past_path_angle_map).*100,'color',colors(2,:),'LineWidth',2)
xlim([0 180])
xticks([0 90 180])
xtickangle(0)
hold on;
xlabel('|Displacement|')
box off

% Group 2 Trials (Away, or Laser ON)
angle_bins = linspace(0, 180, 19);
angle_bin_centers = mean([angle_bins(1:end-1)' angle_bins(2:end)'],2);
replay_past_path_angle_map = nan(length(angle_bin_centers),1);
replay_future_path_angle_map = nan(length(angle_bin_centers),1);

for j = 1:length(angle_bins)-1
    replay_past_path_angle_map(j) = sum(replay_group2.meanAngDisplacement_pastPath>=angle_bins(j) & replay_group2.meanAngDisplacement_pastPath<angle_bins(j+1));
    replay_future_path_angle_map(j) = sum(replay_group2.meanAngDisplacement_futPath>=angle_bins(j) & replay_group2.meanAngDisplacement_futPath<angle_bins(j+1));
end

nexttile(1)
plot(angle_bin_centers,replay_future_path_angle_map./sum(replay_future_path_angle_map).*100,'color',colors(1,:),'LineStyle',':','LineWidth',2)
xlim([0 180])
xticks([0 90 180])
xtickangle(0)
box off
ylabel('% Events')
xlabel('|Displacement|')
hold on;
nexttile(2)
plot(angle_bin_centers,replay_past_path_angle_map./sum(replay_past_path_angle_map).*100,'color',colors(2,:),'LineStyle',':','LineWidth',2)
xlim([0 180])
xticks([0 90 180])
xtickangle(0)
xlabel('|Displacement|')
box off

group1_n = height(replay_group1)
group2_n = height(replay_group2)

mean_abs_angle_future_group1 = nanmean(replay_group1.meanAngDisplacement_futPath);
sem_abs_angle_future_group1 = nanstd(replay_group1.meanAngDisplacement_futPath)./sqrt(sum(~isnan(replay_group1.meanAngDisplacement_futPath)));
mean_abs_angle_past_group1 = nanmean(replay_group1.meanAngDisplacement_pastPath);
sem_abs_angle_past_group1 = nanstd(replay_group1.meanAngDisplacement_pastPath)./sqrt(sum(~isnan(replay_group1.meanAngDisplacement_pastPath)));
mean_abs_angle_future_group2 = nanmean(replay_group2.meanAngDisplacement_futPath);
sem_abs_angle_future_group2 = nanstd(replay_group2.meanAngDisplacement_futPath)./sqrt(sum(~isnan(replay_group2.meanAngDisplacement_futPath)));
mean_abs_angle_past_group2 = nanmean(replay_group2.meanAngDisplacement_pastPath);
sem_abs_angle_past_group2 = nanstd(replay_group2.meanAngDisplacement_pastPath)./sqrt(sum(~isnan(replay_group2.meanAngDisplacement_pastPath)));

data = [mean_abs_angle_future_group1, mean_abs_angle_future_group2, mean_abs_angle_past_group1, mean_abs_angle_past_group2];
err = [sem_abs_angle_future_group1, sem_abs_angle_future_group2, sem_abs_angle_past_group1, sem_abs_angle_past_group2];

ax3 = nexttile(3);
[p,h,z] = ranksum(replay_group1.meanAngDisplacement_futPath,replay_group2.meanAngDisplacement_futPath)
[p,h,z] = ranksum(replay_group1.meanAngDisplacement_pastPath,replay_group2.meanAngDisplacement_pastPath)

colors = [.4660 0.6740 0.1880; 0.4940 0.1840 0.5560;];
transparency_pcnt = 0.5;
colors_2 = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
    [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];
b1 = bar(1,data(1)); hold on;
e1 = errorbar(1,data(1),err(1),'k','linestyle','none');
e1.CapSize = 4;
hold on
b2 = bar(2,data(2)); hold on;
e2 = errorbar(2,data(2),err(2),'k','linestyle','none');
e2.CapSize = 4;
b3 = bar(4,data(3)); hold on;
e3 = errorbar(4,data(3),err(3),'k','linestyle','none');
e3.CapSize = 4;
b4 = bar(5,data(4)); hold on;
e4 = errorbar(5,data(4),err(4),'k','linestyle','none');
e4.CapSize = 4;
b1.FaceColor = colors(1,:);
b1.EdgeColor = 'none';
b2.FaceColor = colors_2(1,:);
b2.EdgeColor = 'none';
b3.FaceColor = colors(2,:);
b3.EdgeColor = 'none';
b4.FaceColor = colors_2(2,:);
b4.EdgeColor = 'none';


box off
xticks([1.5 4.5])
xticklabels({})
ylim(fig_ylim)
ylabel('|Displacement|')

set(gcf, 'Color', 'white','Renderer','painters', 'PaperPositionMode', 'auto');
fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures';
saveas(gcf,fullfile(fig_path,['hist_angles_without_time_' group1 '_' group2 ]),'jpg')
saveas(gcf,fullfile(fig_path,['hist_angles_without_time_' group1 '_' group2 ]),'pdf')


figure('Position',[1986 1051 100 100])
b1 = bar(1,data(1)); hold on;
e1 = errorbar(1,data(1),err(1),'k','linestyle','none');
e1.CapSize = 4;
hold on
b2 = bar(2,data(2)); hold on;
e2 = errorbar(2,data(2),err(2),'k','linestyle','none');
e2.CapSize = 4;
b3 = bar(4,data(3)); hold on;
e3 = errorbar(4,data(3),err(3),'k','linestyle','none');
e3.CapSize = 4;
b4 = bar(5,data(4)); hold on;
e4 = errorbar(5,data(4),err(4),'k','linestyle','none');
e4.CapSize = 4;
b1.FaceColor = colors(1,:);
b1.EdgeColor = 'none';
b2.FaceColor = colors_2(1,:);
b2.EdgeColor = 'none';
b3.FaceColor = colors(2,:);
b3.EdgeColor = 'none';
b4.FaceColor = colors_2(2,:);
b4.EdgeColor = 'none';


box off
xticks([1.5 4.5])
xticklabels({})
ylim(fig_ylim)
ylabel('|Displacement|')


set(gcf, 'Color', 'white','Renderer','painters', 'PaperPositionMode', 'auto');
fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures';
saveas(gcf,fullfile(fig_path,['summary_angles_without_time_' group1 '_' group2 ]),'jpg')
saveas(gcf,fullfile(fig_path,['summary_angles_without_time_' group1 '_' group2 ]),'pdf')

data_table_sub = table();
data_table_sub.angles = [replay_group1.meanAngDisplacement_futPath; replay_group2.meanAngDisplacement_futPath; replay_group1.meanAngDisplacement_pastPath; replay_group2.meanAngDisplacement_pastPath];
data_table_sub.group = [zeros(height(replay_group1),1); ones(height(replay_group2),1); zeros(height(replay_group1),1); ones(height(replay_group2),1)];
data_table_sub.fut_past = [ones(height(replay_group1),1); ones(height(replay_group2),1); 2*ones(height(replay_group1),1); 2*ones(height(replay_group2),1)];
groupOrder = ["1","2"];
namedGroup = categorical(data_table_sub.group,0:1,groupOrder);
directionOrder = ["Future","Past"];
namedDirection = categorical(data_table_sub.fut_past,[1,2],directionOrder);
data_sub = data_table_sub.angles;

% Anova effects:
% [p,tbl,stats] = anovan(data_sub',{namedGroup,namedDirection},'model','interaction','varnames',{'namedGroup','namedDirection'});
% [results,~,~,gnames] = multcompare(stats,"Dimension",[1 2])
%[results,~,~,gnames] = multcompare(stats,"Dimension",[1 2],'CType','lsd')

% NANs removed:
n1 = sum(~isnan(replay_group1.meanAngDisplacement_futPath));
n2 = sum(~isnan(replay_group1.meanAngDisplacement_pastPath));
n3 = sum(~isnan(replay_group2.meanAngDisplacement_futPath));
n4 = sum(~isnan(replay_group2.meanAngDisplacement_pastPath));
n_total = n1+n2+n3+n4;

n_total = 2*height(replay_group1)+2*height(replay_group2);
%% Shuffle (specifically to test significance of interaction)
data1_fut = replay_group1.meanAngDisplacement_futPath;
data2_fut = replay_group2.meanAngDisplacement_futPath;

data1_past = replay_group1.meanAngDisplacement_pastPath;
data2_past = replay_group2.meanAngDisplacement_pastPath;

num_shuffles = 10; %10000
real_interaction_diff = (nanmean(data2_past)-nanmean(data1_past)) - (nanmean(data2_fut)-nanmean(data1_fut));
real_future_diff = nanmean(data2_fut)-nanmean(data1_fut);
real_past_diff = nanmean(data2_past)-nanmean(data1_past);

shuffled_diffs = nan(num_shuffles,1);
shuffled_future_diff = nan(num_shuffles,1);
shuffled_past_diff = nan(num_shuffles,1);

all_replays = [replay_group1; replay_group2];
for shuffle = 1:num_shuffles
    data1_fut_shuffle_inds = randperm(height(all_replays),height(data1_fut));
    data2_fut_shuffle_inds = setdiff([1:height(all_replays)],data1_fut_shuffle_inds);

    data1_past_shuffle_inds = randperm(height(all_replays),height(data1_past));
    data2_past_shuffle_inds = setdiff([1:height(all_replays)],data1_past_shuffle_inds);

    shuffled_future_diff(shuffle) = nanmean(all_replays.meanAngDisplacement_futPath(data1_fut_shuffle_inds)) - nanmean(all_replays.meanAngDisplacement_futPath(data2_fut_shuffle_inds));
    shuffled_past_diff(shuffle) = nanmean(all_replays.meanAngDisplacement_pastPath(data1_past_shuffle_inds)) - nanmean(all_replays.meanAngDisplacement_pastPath(data2_past_shuffle_inds));
end
shuffled_interaction_differences = shuffled_past_diff - shuffled_future_diff;

shuffled_interaction_pval = (sum(shuffled_interaction_differences > abs(real_interaction_diff)) + sum(shuffled_interaction_differences < -1*abs(real_interaction_diff)) + 1 )/(num_shuffles+1);
shuffled_future_pval = (sum(shuffled_future_diff > abs(real_future_diff)) + sum(shuffled_future_diff < -1*abs(real_future_diff)) + 1 )/(num_shuffles+1);
shuffled_past_pval = (sum(shuffled_past_diff > abs(real_past_diff)) + sum(shuffled_past_diff < -1*abs(real_past_diff)) + 1 )/(num_shuffles+1);

[p,h,z] = ranksum(data1_fut,data2_fut)
[p,h,z] = ranksum(data1_past,data2_past)

n_group1_fut = sum(~isnan(data1_fut))
n_group2_fut = sum(~isnan(data2_fut))

n_group1_past = sum(~isnan(data1_past))
n_group2_past = sum(~isnan(data2_past))

n_group1_with_nans = height(replay_group1)
n_group2_with_nans = height(replay_group2)