% Setting axis limits for the plots
ylim_deg = [60 120]; % main axis limits for the y-axis
graph_x_low = 0; graph_x_high = 10; % x-axis limits for plotting
ylabels = {'|Displacement|'}; % Label for y-axis
fig_title = 'replay_angle_past_future_over_time'; % Figure title for saving the plot
% Defining colors for different data categories
color_future = [.4660 0.6740 0.1880]; 
color_past = [0.4940 0.1840 0.556];
color_shuffle_future = [0 0 0]; 
color_shuffle_past = [0.3 0.3 0.3]; 

% Plot the angle between replay and future path, and between replay and past path over time during the stopping period.


%% Binning data: creating time windows to aggregate data for each time bin
windowSize = 2; % Size of the window for binning in seconds
windowShift = 0.5; % Shift between windows in seconds
bin_start = (start_time:windowShift:(end_time-windowSize))'; % Start times for bins
bin_end = bin_start + windowSize; % End times for bins
bin_edges = [bin_start bin_end]; % Edges of the bins
bin_centers = mean(bin_edges, 2); % Centers of the bins

% Properties to be plotted
properties_to_plot = {'angDisplacement_futPath', 'angDisplacement_pastPath', 'meanAngDisplacement_futPath', 'meanAngDisplacement_pastPath'};
binned_properties = table(); % Table to store binned data

% Binning the data for each property
for property = 1:length(properties_to_plot)
    binned_data = table(); % Table to hold binned data for each property
    for time_bin = 1:length(bin_edges)
        % If data is a cell array, convert it to numeric array
        if iscell(replay.(properties_to_plot{property}))
        binned_data.data{time_bin} = abs(cell2mat(replay.(properties_to_plot{property})(replay.time_into_stopping_period>= bin_edges(time_bin,1) & replay.time_into_stopping_period<bin_edges(time_bin,2))));
        else
              binned_data.data{time_bin} = replay.(properties_to_plot{property})(replay.time_into_stopping_period>= bin_edges(time_bin,1) & replay.time_into_stopping_period<bin_edges(time_bin,2));
        end
         % Compute statistics for each bin
        binned_data.mean(time_bin) = nanmean(binned_data.data{time_bin});
        binned_data.median(time_bin) = nanmedian(binned_data.data{time_bin});
        binned_data.sem(time_bin) = nanstd(binned_data.data{time_bin})/sqrt(sum(~isnan(binned_data.data{time_bin})));
        binned_data.n(time_bin) = sum(~isnan(binned_data.data{time_bin}));
        binned_data.std(time_bin) = nanstd(binned_data.data{time_bin});
        binned_data.low(time_bin) = quantile(binned_data.data{time_bin},0.025);
        binned_data.high(time_bin) = quantile(binned_data.data{time_bin},0.975);
        binned_data.sum(time_bin) = nansum(binned_data.data{time_bin});
        binned_data.pcnt(time_bin) = nansum(binned_data.data{time_bin})./length(binned_data.data{time_bin});
    end
     % Store binned data in table
    binned_properties.(properties_to_plot{property}) = binned_data; 
end

% Calculate difference between future and past angles for each bin
binned_data = table();
for time_bin = 1:length(bin_centers)
    a = binned_properties.meanAngDisplacement_futPath(time_bin,:).data{:};
    b = binned_properties.meanAngDisplacement_pastPath(time_bin,:).data{:};
    c = a-b;
    binned_data.data{time_bin} = c;
    binned_data.mean(time_bin) = nanmean(c);
    binned_data.n(time_bin) = sum(~isnan(c));
    binned_data.sem(time_bin) = nanstd(c)/sqrt(sum(~isnan(c)));
    binned_data.std(time_bin) = nanstd(c);
    binned_data.low(time_bin) = quantile(c,0.025);
    binned_data.high(time_bin) = quantile(c,0.975);
end
binned_properties.future_minus_past_angle = binned_data; % Store difference
 
% Repeat the same for calculating "crossings" 
binned_data = table();
for time_bin = 1:length(bin_centers)
    a = binned_properties.angDisplacement_futPath(time_bin,:).data{:};
    b = binned_properties.angDisplacement_pastPath(time_bin,:).data{:};
    c = a-b;
    binned_data.data{time_bin} = c;
    binned_data.mean(time_bin) = nanmean(c);
    binned_data.n(time_bin) = sum(~isnan(c));
    binned_data.sem(time_bin) = nanstd(c)/sqrt(sum(~isnan(c)));
    binned_data.std(time_bin) = nanstd(c);
    binned_data.low(time_bin) = quantile(c,0.025);
    binned_data.high(time_bin) = quantile(c,0.975);
end
binned_properties.crossings_future_minus_past_angle = binned_data; % Store "crossings" data


%% Plots mean angular displacement from past and future paths over time, and the difference in angular displacement (future-past) over time.
fig = figure();
fig.Position = [600 600 300 250];
tiledlayout(2,2)

% Plot 1: Mean angular displacement of future and past paths
ax1 = nexttile(1);
h = shadedErrorBar(bin_centers,binned_properties.meanAngDisplacement_futPath.mean,...
    binned_properties.meanAngDisplacement_futPath.sem,'lineprops','g'); hold on;
h.patch.FaceColor = color_future;
h.mainLine.Color = color_future;
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';

h = shadedErrorBar(bin_centers,binned_properties.meanAngDisplacement_pastPath.mean,...
    binned_properties.meanAngDisplacement_pastPath.sem,'lineprops','m'); hold on;
h.patch.FaceColor = color_past;
h.mainLine.Color = color_past;
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
xticks(graph_x_low:2:graph_x_high)
xlim([graph_x_low graph_x_high])
ylabel(ylabels{1});
ylim(ylim_deg);
yticks([30 60 90 120])

% [rho1,p1] = nancorr(replay.time_into_stopping_period(replay.time_into_stopping_period<=10),replay.meanAngDisplacement_futPath(replay.time_into_stopping_period<=10));
% [rho2,p2] = nancorr(replay.time_into_stopping_period(replay.time_into_stopping_period<=10),replay.meanAngDisplacement_pastPath(replay.time_into_stopping_period<=10));

% Perform statistical comparisons (e.g., paired sign-rank test, unpaired ranksum test)
pval_paired = nan(size(bin_centers));
pval_unpaired = nan(size(bin_centers));

for time_bin = 1:length(bin_centers)
    a = binned_properties.meanAngDisplacement_futPath(time_bin,:).data{:};
    b = binned_properties.meanAngDisplacement_pastPath(time_bin,:).data{:};
    if ~isempty(a) && sum(~isnan(a)) > 0 && ~isempty(b) && sum(~isnan(b)) > 0
        c = binned_properties.future_minus_past_angle(time_bin,:).data{:};
        if sum(~isnan(c))>2
            pval_paired(time_bin) = signrank(c);  % Paired sign-rank test
            pval_unpaired(time_bin) = ranksum(a,b); % Unpaired ranksum test
        else
            pval_paired(time_bin) = nan;
            pval_unpaired(time_bin) = nan;
        end
    else
        pval_paired(time_bin) = nan;
        pval_unpaired(time_bin) = nan;
    end
end

xticks(graph_x_low:2:graph_x_high)
xtickangle(0)


% Plot p-values on the figure (significance indicated by dots)
plot(bin_centers(pval_unpaired<0.05),ylim_deg(2)*ones(sum(pval_unpaired<0.05),1),'.k')

% Plot 2: Plot difference in angular displacement (future-past) for each
% replay, over time
ax2 = nexttile(2);
h = shadedErrorBar(bin_centers,binned_properties.future_minus_past_angle.mean,...
    binned_properties.future_minus_past_angle.sem,'lineprops','k'); hold on;
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
hold on
plot([graph_x_low graph_x_high],[0 0],'--','color',[0.5 0.5 0.5],'LineWidth',1); hold on;
ylim(ylim_deg_difference_plot);
yticks([ylim_deg_difference_plot(1) 0 ylim_deg_difference_plot(2)])
xticks(graph_x_low:2:graph_x_high)
xlim([graph_x_low graph_x_high])
xtickangle(0)
ylabel('|Displacement|');
plot(bin_centers(pval_unpaired<0.05),ylim_deg_difference_plot(2)*ones(sum(pval_unpaired<0.05),1),'.k')

% Visualization only
ax3 = nexttile(3);
% ax1Chil = ax1.Children;
% copyobj(ax1Chil,ax3)
xticks(graph_x_low:2:graph_x_high)
xlim([graph_x_low graph_x_high])
xtickangle(0)
ylim(ylim_deg);
xlabel('Time since arrival (s)')
ylabel('|Displacement|');
yticks([30 60 90 120])

% Visualization only
ax4 = nexttile(4);
% ax1Chil = ax1.Children;
% copyobj(ax1Chil,ax4)
xticks(graph_x_low:2:graph_x_high)
xlim([graph_x_low graph_x_high])
xtickangle(0)
ylim(ylim_deg);
xlabel('Time since arrival (s)')
ylabel('|Displacement|');
yticks([30 60 90 120])

% Export and save the figure
set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
export_fig(fullfile(fig_path,[fig_title '_disp' num2str(replay_dispersionThr,2) '_rats_' num2str(rats)]),'-jpeg')
saveas(gcf,fullfile(fig_path,[fig_title '_disp' num2str(replay_dispersionThr,2) '_rats_' num2str(rats)]),'pdf')

total_n_nans_included = height(replay(replay.time_since_real_drink_onset >0 & replay.time_since_real_drink_onset<=10,:));
total_n_rf = sum(~isnan(replay.meanAngDisplacement_futPath(replay.time_since_real_drink_onset >0 & replay.time_since_real_drink_onset<=10)));
min_n_rf = min([min(binned_properties.meanAngDisplacement_futPath.n(bin_centers<=10))]);
max_n_rf = max([max(binned_properties.meanAngDisplacement_futPath.n(bin_centers<=10))]);
total_n_rp = sum(~isnan(replay.meanAngDisplacement_pastPath(replay.time_since_real_drink_onset >0 & replay.time_since_real_drink_onset<=10)));
min_n_rp = min([min(binned_properties.meanAngDisplacement_pastPath.n(bin_centers<=10))]);
max_n_rp = max([max(binned_properties.meanAngDisplacement_pastPath.n(bin_centers<=10))]);

%% Crossings: As above, but treats each 'crossing' as a unique datapoint (rather than averaging all crossings per replay).
fig = figure();
fig.Position = [600 600 300 250];
tiledlayout(2,2)
ax1 = nexttile(1);
h = shadedErrorBar(bin_centers,binned_properties.angDisplacement_futPath.mean,...
    binned_properties.angDisplacement_futPath.sem,'lineprops','g'); hold on;
h.patch.FaceColor = color_future;
h.mainLine.Color = color_future;
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';

h = shadedErrorBar(bin_centers,binned_properties.angDisplacement_pastPath.mean,...
    binned_properties.angDisplacement_pastPath.sem,'lineprops','m'); hold on;
h.patch.FaceColor = color_past;
h.mainLine.Color = color_past;
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
xticks(graph_x_low:2:graph_x_high)
xlim([graph_x_low graph_x_high])


% compare future and past similarity for each time bin
pval_paired = nan(size(bin_centers));
pval_unpaired = nan(size(bin_centers));

for time_bin = 1:length(bin_centers)
    a = binned_properties.angDisplacement_futPath(time_bin,:).data{:};
    b = binned_properties.angDisplacement_pastPath(time_bin,:).data{:};
    if ~isempty(a) && sum(~isnan(a)) > 0 && ~isempty(b) && sum(~isnan(b)) > 0
        c = binned_properties.crossings_future_minus_past_angle(time_bin,:).data{:};
        if sum(~isnan(c))>2
            pval_paired(time_bin) = signrank(c);
            pval_unpaired(time_bin) = ranksum(a,b);
        else
            pval_paired(time_bin) = nan;
            pval_unpaired(time_bin) = nan;
        end
    else
        pval_paired(time_bin) = nan;
        pval_unpaired(time_bin) = nan;
    end
end


xticks(graph_x_low:2:graph_x_high)
xtickangle(0)

ylim(ylim_deg);
yticks([30 60 90 120])
plot(bin_centers(pval_unpaired<0.05),ylim_deg(2)*ones(sum(pval_unpaired<0.05),1),'.k')
ylabel(ylabels{1});

ax2 = nexttile(2);

h = shadedErrorBar(bin_centers,binned_properties.crossings_future_minus_past_angle.mean,...
    binned_properties.crossings_future_minus_past_angle.sem,'lineprops','k'); hold on;
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
hold on
plot([graph_x_low graph_x_high],[0 0],'--','color',[0.5 0.5 0.5],'LineWidth',1); hold on;
ylim(ylim_deg_difference_plot);
yticks([ylim_deg_difference_plot(1) 0 ylim_deg_difference_plot(2)])
xticks(graph_x_low:2:graph_x_high)
xlim([graph_x_low graph_x_high])
xtickangle(0)
ylabel('|Displacement|');
plot(bin_centers(pval_unpaired<0.05),ylim_deg_difference_plot(2)*ones(sum(pval_unpaired<0.05),1),'.k')

ax3 = nexttile(3);
% ax1Chil = ax1.Children;
% copyobj(ax1Chil,ax3)
xticks(graph_x_low:2:graph_x_high)
xlim([graph_x_low graph_x_high])
xtickangle(0)
ylim(ylim_deg);
xlabel('Time since arrival (s)')
ylabel('|Displacement|');
yticks([30 60 90 120])

ax4 = nexttile(4);
% ax1Chil = ax1.Children;
% copyobj(ax1Chil,ax4)
xticks(graph_x_low:2:graph_x_high)
xlim([graph_x_low graph_x_high])
xtickangle(0)
ylim(ylim_deg);
xlabel('Time since arrival (s)')
ylabel('|Displacement|');
yticks([30 60 90 120])

set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
export_fig(fullfile(fig_path,[fig_title '_disp' num2str(replay_dispersionThr,2) '_rats_' num2str(rats) '_crossings']),'-jpeg')
saveas(gcf,fullfile(fig_path,[fig_title '_disp' num2str(replay_dispersionThr,2) '_rats_' num2str(rats) '_crossings']),'pdf')

min_n = min([min(binned_properties.crossings_future_minus_past_angle.n(bin_centers<=10)) min(binned_properties.crossings_future_minus_past_angle.n(bin_centers<=10))]);
max_n = max([max(binned_properties.crossings_future_minus_past_angle.n(bin_centers<=10)) max(binned_properties.crossings_future_minus_past_angle.n(bin_centers<=10))]);

%% Crossings and replays side by side 
fig = figure();
fig.Position = [600 600 300 250];
tiledlayout(2,2)
ax1 = nexttile(1);
h = shadedErrorBar(bin_centers,binned_properties.angDisplacement_futPath.mean,...
    binned_properties.angDisplacement_futPath.sem,'lineprops','g'); hold on;
h.patch.FaceColor = color_future;
h.mainLine.Color = color_future;
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';

h = shadedErrorBar(bin_centers,binned_properties.angDisplacement_pastPath.mean,...
    binned_properties.angDisplacement_pastPath.sem,'lineprops','m'); hold on;
h.patch.FaceColor = color_past;
h.mainLine.Color = color_past;
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
xticks(graph_x_low:2:graph_x_high)
xlim([graph_x_low graph_x_high])

% compare future and past similarity for each time bin
pval_paired = nan(size(bin_centers));
pval_unpaired = nan(size(bin_centers));

for time_bin = 1:length(bin_centers)
    a = binned_properties.angDisplacement_futPath(time_bin,:).data{:};
    b = binned_properties.angDisplacement_pastPath(time_bin,:).data{:};
    if ~isempty(a) && sum(~isnan(a)) > 0 && ~isempty(b) && sum(~isnan(b)) > 0
        c = binned_properties.crossings_future_minus_past_angle(time_bin,:).data{:};
        if sum(~isnan(c))>2
            pval_paired(time_bin) = signrank(c);
            pval_unpaired(time_bin) = ranksum(a,b);
        else
            pval_paired(time_bin) = nan;
            pval_unpaired(time_bin) = nan;
        end
    else
        pval_paired(time_bin) = nan;
        pval_unpaired(time_bin) = nan;
    end
end

xticks(graph_x_low:2:graph_x_high)
xtickangle(0)

ylim(ylim_deg);
yticks([30 60 90 120])
plot(bin_centers(pval_unpaired<0.05),ylim_deg(2)*ones(sum(pval_unpaired<0.05),1),'.k')
ylabel(ylabels{1}); 

ax2 = nexttile(2);
h = shadedErrorBar(bin_centers,binned_properties.meanAngDisplacement_futPath.mean,...
    binned_properties.meanAngDisplacement_futPath.sem,'lineprops','g'); hold on;
h.patch.FaceColor = color_future;
h.mainLine.Color = color_future;
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';

h = shadedErrorBar(bin_centers,binned_properties.meanAngDisplacement_pastPath.mean,...
    binned_properties.meanAngDisplacement_pastPath.sem,'lineprops','m'); hold on;
h.patch.FaceColor = color_past;
h.mainLine.Color = color_past;
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
xticks(graph_x_low:2:graph_x_high)
xlim([graph_x_low graph_x_high])

[rho1,p1] = nancorr(replay.time_into_stopping_period(replay.time_into_stopping_period<10),replay.meanAngDisplacement_futPath(replay.time_into_stopping_period<10));
[rho2,p2] = nancorr(replay.time_into_stopping_period(replay.time_into_stopping_period<10),replay.meanAngDisplacement_pastPath(replay.time_into_stopping_period<10));
ylabel(ylabels{1});


% compare future and past similarity for each time bin
pval_paired = nan(size(bin_centers));
pval_unpaired = nan(size(bin_centers));

for time_bin = 1:length(bin_centers)
    a = binned_properties.meanAngDisplacement_futPath(time_bin,:).data{:};
    b = binned_properties.meanAngDisplacement_pastPath(time_bin,:).data{:};
    if ~isempty(a) && sum(~isnan(a)) > 0 && ~isempty(b) && sum(~isnan(b)) > 0
        c = binned_properties.future_minus_past_angle(time_bin,:).data{:};
        if sum(~isnan(c))>2
            pval_paired(time_bin) = signrank(c);
            pval_unpaired(time_bin) = ranksum(a,b);
        else
            pval_paired(time_bin) = nan;
            pval_unpaired(time_bin) = nan;
        end
    else
        pval_paired(time_bin) = nan;
        pval_unpaired(time_bin) = nan;
    end
end

% plot([graph_x_low graph_x_high],[90 90],'--','color',[0.5 0.5 0.5],'LineWidth',2); hold on;
xticks(graph_x_low:2:graph_x_high)
xtickangle(0)

ylim(ylim_deg);
yticks([30 60 90 120])
plot(bin_centers(pval_unpaired<0.05),ylim_deg(2)*ones(sum(pval_unpaired<0.05),1),'.k')
ylabel(ylabels{1});

set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
saveas(gcf,fullfile(fig_path,['crossings_and_replays_rats_' num2str(rats) '_crossings']),'pdf')

a = abs(cell2mat(replay.angDisplacement_futPath(replay.time_into_stopping_period>= 0 & replay.time_into_stopping_period< 10)));
b = abs(cell2mat(replay.angDisplacement_pastPath(replay.time_into_stopping_period>= 0 & replay.time_into_stopping_period< 10)));

sum(~isnan(a))
sum(~isnan(b))
sum(replay.time_into_stopping_period>0 & replay.time_into_stopping_period<=10)