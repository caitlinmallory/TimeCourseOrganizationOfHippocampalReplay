% Plots Figure 3F 
ylim_deg = [70 110];
ylabels = {'|Displacement|'};
fig_title = 'replay_angle_past_versus_shuffle';

properties = {'meanAngDisplacement_futPath'; 'meanAngDisplacement_pastPath'; 'mean_shuffle_all'; 'mean_shuffle_future'; 'mean_shuffle_past'};
windowSize = 2;
windowShift = 0.5;
bin_start = (start_time:windowShift:(end_time-windowSize))';
bin_end = bin_start + windowSize;
bin_edges = [bin_start bin_end];
bin_centers = mean(bin_edges,2);

binned_properties = table();
for property = 1:length(properties)
    binned_data = table();
    for time_bin = 1:length(bin_edges)
        binned_data.data{time_bin} = replay.(properties{property})(replay.time_into_stopping_period>= bin_edges(time_bin,1) & replay.time_into_stopping_period<bin_edges(time_bin,2));
        binned_data.mean(time_bin) = nanmean(binned_data.data{time_bin});
        binned_data.median(time_bin) = nanmedian(binned_data.data{time_bin});

        binned_data.sem(time_bin) = nanstd(binned_data.data{time_bin})/sqrt(sum(~isnan(binned_data.data{time_bin})));
        binned_data.n(time_bin) = sum(~isnan(binned_data.data{time_bin}));
        binned_data.low(time_bin) = quantile(binned_data.data{time_bin},0.025);
        binned_data.high(time_bin) = quantile(binned_data.data{time_bin},0.975);
        binned_data.sum(time_bin) = nansum(binned_data.data{time_bin});
        binned_data.pcnt(time_bin) = nansum(binned_data.data{time_bin})./length(binned_data.data{time_bin});
    end
    binned_properties.(properties{property}) = binned_data;
end

binned_data = table();
for time_bin = 1:length(bin_edges)
    a = binned_properties.meanAngDisplacement_futPath(time_bin,:).data{:};
    b = binned_properties.meanAngDisplacement_pastPath(time_bin,:).data{:};
    c = a-b;
    binned_data.data{time_bin} = c;
    binned_data.mean(time_bin) = nanmean(c);
    binned_data.sem(time_bin) = nanstd(c)/sqrt(sum(~isnan(c)));
    binned_data.n(time_bin) = sum(~isnan(c));
    binned_data.low(time_bin) = quantile(c,0.025);
    binned_data.high(time_bin) = quantile(c,0.975);
end
binned_properties.future_minus_past_angle = binned_data;

sum(~isnan(replay.meanAngDisplacement_pastPath(replay.time_since_real_drink_onset>0 & replay.time_since_real_drink_onset<=10)) & ...
    ~isnan(replay.meanAngDisplacement_pastPath(replay.time_since_real_drink_onset>0 & replay.time_since_real_drink_onset<=10)) )

colors = [.4660 0.6740 0.1880; 0.4940 0.1840 0.5560;];

color_shuffle_future = [0 0 0];
color_shuffle_past = [0.3 0.3 0.3];
fig = figure();
fig.Position = [600 600 300 250];
tiledlayout(2,2)
ax1 = nexttile();
h = shadedErrorBar(bin_centers,binned_properties.meanAngDisplacement_pastPath.mean,...
    binned_properties.meanAngDisplacement_pastPath.sem,'lineprops','m'); hold on;
h.patch.FaceColor = colors(2,:);
h.mainLine.Color = colors(2,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
% Option 3: plot mean and sem of shuffled values at each time bin:
if plot_shuffles == 1 && strcmp(shuffle_type,'path_swap')
    h = shadedErrorBar(bin_centers,binned_properties.('mean_shuffle_all').mean,...
        binned_properties.('mean_shuffle_all').sem,'lineprops','k'); hold on;
    h.patch.FaceColor = [0.5 0.5 0.5];
    h.mainLine.Color = [0.5 0.5 0.5];
    h.mainLine.LineWidth = 1;
    h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
elseif plot_shuffles == 1 && strcmp(shuffle_type,'path_swap_future_past_preserved')
    h = shadedErrorBar(bin_centers,binned_properties.('mean_shuffle_future').mean,...
        binned_properties.('mean_shuffle_future').sem,'lineprops','k'); hold on;
    h.patch.FaceColor = color_shuffle_future;
    h.mainLine.Color = color_shuffle_future;
    h.mainLine.LineWidth = 1;
    h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
    h = shadedErrorBar(bin_centers,binned_properties.('mean_shuffle_past').mean,...
        binned_properties.('mean_shuffle_past').sem,'lineprops','k'); hold on;
    h.patch.FaceColor = color_shuffle_past;
    h.mainLine.Color = color_shuffle_past;
    h.mainLine.LineWidth = 1;
    h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
end
xticks(graph_x_low:2:graph_x_high)
xtickangle(0)
xticks()
xticks(graph_x_low:2:graph_x_high)
xlim([graph_x_low graph_x_high])
ylim(ylim_deg)
ylabel(ylabels{1});
% compare future and past similarity for each time bin
pval_paired = nan(size(bin_centers));
pval_unpaired = nan(size(bin_centers));
plot_n = nan(size(bin_centers));
for time_bin = 1:length(bin_centers)
    a = binned_properties.meanAngDisplacement_pastPath(time_bin,:).data{:};
    b = binned_properties.mean_shuffle_all(time_bin,:).data{:};
    plot_n(time_bin) = sum(~isnan(a));
    if ~isempty(a) && sum(~isnan(a)) > 0 && ~isempty(b) && sum(~isnan(b)) > 0
        pval_paired(time_bin) = signrank(a,b);
        pval_unpaired(time_bin) = ranksum(a,b);
    end
end

min(plot_n(bin_centers<=10))
max(plot_n(bin_centers<=10))
n_data_points = [min(plot_n(bin_centers<=10)) max(plot_n(bin_centers<=10))]
ot(bin_centers(pval_paired<0.05),ylim_deg(2)*ones(sum(pval_paired<0.05),1),'.k')

ax2 = nexttile();
ax1Chil = ax1.Children;
copyobj(ax1Chil,ax2);
xlim([0 10]);
xticks([0:2:10]);
xtickangle(0)
ylabel('|Displacement|')
ylim(ylim_deg)
ax3 = nexttile();
ax1Chil = ax1.Children;
copyobj(ax1Chil,ax3);
xlabel('Time since arrival (s)')
xlim([0 10]);
xticks([0:2:10]);
xtickangle(0)
ylabel('|Displacement|')
ylim(ylim_deg)
ax4 = nexttile();
ax1Chil = ax1.Children;
copyobj(ax1Chil,ax4);
xlabel('Time since arrival (s)')
xlim([0 10]);
xticks([0:2:10]);
xtickangle(0)
ylabel('|Displacement|')
ylim(ylim_deg)

set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
export_fig(fullfile(fig_path,[fig_title '_disp' num2str(replay_dispersionThr,2) '_rats_' num2str(rats)]),'-jpeg')
saveas(gcf,fullfile(fig_path,[fig_title '_disp' num2str(replay_dispersionThr,2) '_rats_' num2str(rats)]),'pdf')