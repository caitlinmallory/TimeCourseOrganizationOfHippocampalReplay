%% Plots angular displacements from past or future path over time, for 2 different groups. 
% Used to generate figure 3D
sig_test = 'medians';
start_time = 0;
end_time = 10;
windowSize = 2;
windowShift = 0.5;

home_only=0;
away_only=0;

ylim_deg = [40 120];
ylim_deg_ticks = ylim_deg(1):40:ylim_deg(2);
ylim_deg_diff_plot = [-60 40];

category1 = 'tight';
category2 = 'wide';
angle_bins = [0 60; 120 180];
inds1 = find(replay.angle_between_past_future_trajectory >=  angle_bins(1,1) & replay.angle_between_past_future_trajectory<=angle_bins(1,2));
inds2 = find(replay.angle_between_past_future_trajectory >=  angle_bins(2,1) & replay.angle_between_past_future_trajectory<=angle_bins(2,2));
properties = {'meanAngDisplacement_futPath'; 'meanAngDisplacement_pastPath'};

% category1 = 'home';
% category2 = 'away';
% inds1 = find(replay.home_event==1);
% inds2 = find(replay.home_event==0);

% category1 = 'off';
% category2 = 'on';
% inds1 = find(replay.laser_state_binary==0);
% inds2 = find(replay.laser_state_binary==1);
% if home_only==1
% inds1 = find(replay.laser_state_binary==0 & replay.home_event==1);
% inds2 = find(replay.laser_state_binary==1 & replay.home_event==1);
% end
% if away_only==1
% inds1 = find(replay.laser_state_binary==0 & replay.home_event==0);
% inds2 = find(replay.laser_state_binary==1 & replay.home_event==0);
% end

bin_start = (start_time:windowShift:(end_time-windowSize))';
bin_end = bin_start + windowSize;

bin_edges = [bin_start bin_end];
bin_centers = mean(bin_edges,2);
grouped_binned_properties = struct();

all_inds = [{inds1} {inds2}];
for i = 1:2
    inds = all_inds{i};
    replay_sub = replay(inds,:);
    binned_properties_sub = table();
    for property = 1:length(properties)
        binned_data = table();
        for time_bin = 1:length(bin_edges)
            binned_data.data{time_bin} = replay_sub.(properties{property})(replay_sub.time_into_stopping_period>= bin_edges(time_bin,1) & replay_sub.time_into_stopping_period<bin_edges(time_bin,2));
            binned_data.mean(time_bin) = nanmean(binned_data.data{time_bin});
            binned_data.median(time_bin) = nanmedian(binned_data.data{time_bin});

            binned_data.sem(time_bin) = nanstd(binned_data.data{time_bin})/sqrt(sum(~isnan(binned_data.data{time_bin})));
            binned_data.n(time_bin) = sum(~isnan(binned_data.data{time_bin}));
            binned_data.low(time_bin) = quantile(binned_data.data{time_bin},0.025);
            binned_data.high(time_bin) = quantile(binned_data.data{time_bin},0.975);
            binned_data.sum(time_bin) = nansum(binned_data.data{time_bin});
            binned_data.pcnt(time_bin) = nansum(binned_data.data{time_bin})./length(binned_data.data{time_bin});
        end
        binned_properties_sub.(properties{property}) = binned_data;
    end

    binned_data = table();
    for time_bin = 1:length(bin_centers)
        a = binned_properties_sub.meanAngDisplacement_futPath(time_bin,:).data{:};
        b = binned_properties_sub.meanAngDisplacement_pastPath(time_bin,:).data{:};
        c = b-a;
        binned_data.data{time_bin} = c;
        binned_data.mean(time_bin) = nanmean(c);
        binned_data.n(time_bin) = sum(~isnan(c));
        binned_data.sem(time_bin) = nanstd(c)/sqrt(sum(~isnan(c)));
        binned_data.low(time_bin) = quantile(c,0.025);
        binned_data.high(time_bin) = quantile(c,0.975);
    end
    binned_properties_sub.past_minus_future_angle = binned_data;

    binned_data = table();
    for time_bin = 1:length(bin_centers)
        a = binned_properties_sub.meanAngDisplacement_futPath(time_bin,:).data{:};
        b = binned_properties_sub.meanAngDisplacement_pastPath(time_bin,:).data{:};
        c = a-b;
        binned_data.data{time_bin} = c;
        binned_data.mean(time_bin) = nanmean(c);
        binned_data.sem(time_bin) = nanstd(c)/sqrt(sum(~isnan(c)));
        binned_data.n(time_bin) = sum(~isnan(c));
        binned_data.low(time_bin) = quantile(c,0.025);
        binned_data.high(time_bin) = quantile(c,0.975);
    end
    binned_properties_sub.future_minus_past_angle = binned_data;
    grouped_binned_properties(i).binned_properties = binned_properties_sub;
end

graph_x_low = 0;
graph_x_high = 10;

if align_to_drink_offset == 1
    bin_centers = -1*(bin_centers);
    graph_x_low = -10;
    graph_x_high = 0;
end

replay_sub = replay(inds1,:);
replay_sub(isnan(replay_sub.meanAngDisplacement_futPath) | isnan(replay_sub.meanAngDisplacement_pastPath),:) = [];
replay_sub = replay_sub(replay_sub.time_into_stopping_period >=0 & replay_sub.time_into_stopping_period <= 10,:);
[rho,p] = nancorr(replay_sub.time_into_stopping_period, replay_sub.meanAngDisplacement_futPath)
[rho,p] = nancorr(replay_sub.time_into_stopping_period, replay_sub.meanAngDisplacement_pastPath)
df = height(replay_sub)-2

replay_sub = replay(inds2,:);
replay_sub = replay_sub(replay_sub.time_into_stopping_period >=0 & replay_sub.time_into_stopping_period <= 10,:);
[rho,p] = nancorr(replay_sub.time_into_stopping_period, replay_sub.meanAngDisplacement_futPath)
[rho,p] = nancorr(replay_sub.time_into_stopping_period, replay_sub.meanAngDisplacement_pastPath)
df = height(replay_sub)-2;
%% Version 1:
colors = [.4660 0.6740 0.1880; 0.4940 0.1840 0.5560;];
transparency_pcnt = 1;
colors_2 = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
    [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];

figure('Position',[1921 560 350 125])
tiledlayout(1,3,'TileSpacing','tight')
ax1=nexttile(1);
h = shadedErrorBar(bin_centers,grouped_binned_properties(1).binned_properties.meanAngDisplacement_futPath.mean,...
    grouped_binned_properties(1).binned_properties.meanAngDisplacement_futPath.sem,'lineprops','k'); hold on;
h.patch.FaceColor = colors_2(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
hold on
h = shadedErrorBar(bin_centers,grouped_binned_properties(2).binned_properties.meanAngDisplacement_futPath.mean,...
    grouped_binned_properties(2).binned_properties.meanAngDisplacement_futPath.sem,'lineprops','--k'); hold on;
h.patch.FaceColor = colors_2(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
xlim([0 10])
xticks(0:2:10)
xtickangle(0)
ylim(ylim_deg)
yticks(ylim_deg_ticks);
ylabel('|Displacement|')

pvals = nan(length(bin_centers),1);
for i = 1:length(bin_centers)
    if strcmp(sig_test,'means')
        [h, pvals(i)] = ttest2(grouped_binned_properties(1).binned_properties.meanAngDisplacement_futPath.data{i},grouped_binned_properties(2).binned_properties.meanAngDisplacement_futPath.data{i});
    elseif strcmp(sig_test,'medians')
        [pvals(i)] = ranksum(grouped_binned_properties(1).binned_properties.meanAngDisplacement_futPath.data{i},grouped_binned_properties(2).binned_properties.meanAngDisplacement_futPath.data{i});
    end
end
hold on
ylimit = gca().YLim(2);
plot(bin_centers(pvals<0.05),ylimit*ones(sum(pvals<0.05)),'.k')

ax2=nexttile(2);
h = shadedErrorBar(bin_centers,grouped_binned_properties(1).binned_properties.meanAngDisplacement_pastPath.mean,...
    grouped_binned_properties(1).binned_properties.meanAngDisplacement_pastPath.sem,'lineprops','k'); hold on;
h.patch.FaceColor = colors_2(2,:);
h.mainLine.Color = colors(2,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
h = shadedErrorBar(bin_centers,grouped_binned_properties(2).binned_properties.meanAngDisplacement_pastPath.mean,...
    grouped_binned_properties(2).binned_properties.meanAngDisplacement_pastPath.sem,'lineprops','--k'); hold on;
h.patch.FaceColor = colors_2(2,:);
h.mainLine.Color = colors(2,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
xlabel('Time since arrival (s)')
xlim([0 10])
xticks(0:2:10)
xtickangle(0)
ylim(ylim_deg)
yticks(ylim_deg_ticks);

pvals = nan(length(bin_centers),1);
for i = 1:length(bin_centers)
    if strcmp(sig_test,'means')
        [h, pvals(i)] = ttest2(grouped_binned_properties(1).binned_properties.meanAngDisplacement_pastPath.data{i},grouped_binned_properties(2).binned_properties.meanAngDisplacement_pastPath.data{i});
    elseif strcmp(sig_test,'medians')
        [pvals(i)] = ranksum(grouped_binned_properties(1).binned_properties.meanAngDisplacement_pastPath.data{i},grouped_binned_properties(2).binned_properties.meanAngDisplacement_pastPath.data{i});
    end
end
hold on
ylimit = gca().YLim(2);
plot(bin_centers(pvals<0.05),ylimit*ones(sum(pvals<0.05)),'.k')

% Difference of differences, real versus shuffle:
ax3=nexttile(3);
colors = [0 0 0; 0 0 0];
transparency_pcnt = 1;
colors_2 = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
    [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];

h = shadedErrorBar(bin_centers,grouped_binned_properties(1).binned_properties.future_minus_past_angle.mean,...
    grouped_binned_properties(1).binned_properties.past_minus_future_angle.sem,'lineprops','g'); hold on;
h.patch.FaceColor = colors_2(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
hold on;
h = shadedErrorBar(bin_centers,grouped_binned_properties(2).binned_properties.future_minus_past_angle.mean,...
    grouped_binned_properties(2).binned_properties.past_minus_future_angle.sem,'lineprops','--k'); hold on;
h.patch.FaceColor = colors_2(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
yline(0)
ylim(ylim_deg_diff_plot)

pvals = nan(length(bin_centers),1);
for i = 1:length(bin_centers)
    if strcmp(sig_test,'means')
        [h, pvals(i)] = ttest2(grouped_binned_properties(1).binned_properties.future_minus_past_angle.data{i},grouped_binned_properties(2).binned_properties.future_minus_past_angle.data{i});
    elseif strcmp(sig_test,'medians')
        [pvals(i)] = ranksum(grouped_binned_properties(1).binned_properties.future_minus_past_angle.data{i},grouped_binned_properties(2).binned_properties.future_minus_past_angle.data{i});
    end
end
hold on
ylimit = gca().YLim(2);
plot(bin_centers(pvals<0.05),ylimit*ones(sum(pvals<0.05)),'.k')

xticks(0:2:10)
xtickangle(0)
xlim([0 10])
set(gcf, 'Color', 'white','PaperPositionMode','auto');
fig_title = [category1 '_v_' category2];

saveas(gcf,fullfile(fig_path,[fig_title '_disp' num2str(replay_dispersionThr,2) '_rats_' num2str(rats)]),'jpg')
saveas(gcf,fullfile(fig_path,[fig_title '_disp' num2str(replay_dispersionThr,2) '_rats_' num2str(rats)]),'pdf')

%% Version 2:
colors = [.4660 0.6740 0.1880; 0.4940 0.1840 0.5560;];
transparency_pcnt = 1;
colors_2 = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
    [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];

figure('Position',[1921 560 350 125])
tiledlayout(1,3,'TileSpacing','tight')
nexttile(1)
h = shadedErrorBar(bin_centers,grouped_binned_properties(1).binned_properties.meanAngDisplacement_futPath.mean,...
    grouped_binned_properties(1).binned_properties.meanAngDisplacement_futPath.sem,'lineprops','k'); hold on;
h.patch.FaceColor = colors_2(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
hold on
h = shadedErrorBar(bin_centers,grouped_binned_properties(1).binned_properties.meanAngDisplacement_pastPath.mean,...
    grouped_binned_properties(1).binned_properties.meanAngDisplacement_pastPath.sem,'lineprops','k'); hold on;
h.patch.FaceColor = colors_2(2,:);
h.mainLine.Color = colors(2,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
xlim([0 10])
xticks(0:2:10)
xtickangle(0)
ylim(ylim_deg)
yticks(ylim_deg_ticks);
ylabel('|Displacement|')

pvals = nan(length(bin_centers),1);
for i = 1:length(bin_centers)
    if strcmp(sig_test,'means')
        [h, pvals(i)] = ttest2(grouped_binned_properties(1).binned_properties.meanAngDisplacement_pastPath.data{i},grouped_binned_properties(1).binned_properties.meanAngDisplacement_futPath.data{i})
    elseif strcmp(sig_test,'medians')
        [pvals(i)] = ranksum(grouped_binned_properties(1).binned_properties.meanAngDisplacement_pastPath.data{i},grouped_binned_properties(1).binned_properties.meanAngDisplacement_futPath.data{i})
    end
end
hold on
ylimit = gca().YLim(2);
plot(bin_centers(pvals<0.05),ylimit*ones(sum(pvals<0.05),1),'.k')

nexttile(2)
h = shadedErrorBar(bin_centers,grouped_binned_properties(2).binned_properties.meanAngDisplacement_futPath.mean,...
    grouped_binned_properties(1).binned_properties.meanAngDisplacement_futPath.sem,'lineprops','--k'); hold on;
h.patch.FaceColor = colors_2(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
h = shadedErrorBar(bin_centers,grouped_binned_properties(2).binned_properties.meanAngDisplacement_pastPath.mean,...
    grouped_binned_properties(2).binned_properties.meanAngDisplacement_pastPath.sem,'lineprops','--k'); hold on;
h.patch.FaceColor = colors_2(2,:);
h.mainLine.Color = colors(2,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
xlabel('Time since arrival (s)')
xlim([0 10])
xticks(0:2:10)
xtickangle(0)
ylim(ylim_deg)
yticks(ylim_deg_ticks);
hold on
ylimit = gca().YLim(2);
pvals = nan(length(bin_centers),1);
for i = 1:length(bin_centers)
    if strcmp(sig_test,'means')
        [h, pvals(i)] = ttest2(grouped_binned_properties(2).binned_properties.meanAngDisplacement_pastPath.data{i},grouped_binned_properties(2).binned_properties.meanAngDisplacement_futPath.data{i})
    elseif strcmp(sig_test,'medians')
        [pvals(i)] = ranksum(grouped_binned_properties(2).binned_properties.meanAngDisplacement_pastPath.data{i},grouped_binned_properties(2).binned_properties.meanAngDisplacement_futPath.data{i})
    end
end
plot(bin_centers(pvals<0.05),ylimit*ones(sum(pvals<0.05),1),'.k')
% Difference of differences, real versus shuffle

nexttile(3)
colors = [0 0 0; 0 0 0];
transparency_pcnt = 0.5;
colors_2 = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
    [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];

h = shadedErrorBar(bin_centers,grouped_binned_properties(1).binned_properties.future_minus_past_angle.mean,...
    grouped_binned_properties(1).binned_properties.past_minus_future_angle.sem,'lineprops','g'); hold on;
h.patch.FaceColor = colors_2(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
hold on;
h = shadedErrorBar(bin_centers,grouped_binned_properties(2).binned_properties.future_minus_past_angle.mean,...
    grouped_binned_properties(2).binned_properties.past_minus_future_angle.sem,'lineprops','--k'); hold on;
h.patch.FaceColor = colors_2(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
yline(0)
ylim(ylim_deg_diff_plot)

pvals = nan(length(bin_centers),1);
for i = 1:length(bin_centers)
    if strcmp(sig_test,'means')
        [h, pvals(i)] = ttest2(grouped_binned_properties(1).binned_properties.future_minus_past_angle.data{i},grouped_binned_properties(2).binned_properties.future_minus_past_angle.data{i});
    elseif strcmp(sig_test,'medians')
        [pvals(i)] = ranksum(grouped_binned_properties(1).binned_properties.future_minus_past_angle.data{i},grouped_binned_properties(2).binned_properties.future_minus_past_angle.data{i});
    end
end
hold on
ylimit = gca().YLim(2);
plot(bin_centers(pvals<0.05),ylimit*ones(sum(pvals<0.05),1),'.k')

xticks(0:2:10)
xtickangle(0)
xlim([0 10])
set(gcf, 'Color', 'white','PaperPositionMode','auto');
fig_title = [category1 '_v_' category2];

saveas(gcf,fullfile(fig_path,[fig_title '_v2_disp' num2str(replay_dispersionThr,2) '_rats_' num2str(rats)]),'jpg')
saveas(gcf,fullfile(fig_path,[fig_title '_v2_disp' num2str(replay_dispersionThr,2) '_rats_' num2str(rats)]),'pdf')

%% Version 3 (Overlay)
colors = [.4660 0.6740 0.1880; 0.4940 0.1840 0.5560;];
transparency_pcnt = 1;
colors_2 = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
    [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];

figure('Position',[1921 560 350 125])
tiledlayout(1,3,'TileSpacing','tight')
ax1 = nexttile(1);
h = shadedErrorBar(bin_centers,grouped_binned_properties(1).binned_properties.meanAngDisplacement_futPath.mean,...
    grouped_binned_properties(1).binned_properties.meanAngDisplacement_futPath.sem,'lineprops','k'); hold on;
h.patch.FaceColor = colors_2(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
hold on
h = shadedErrorBar(bin_centers,grouped_binned_properties(1).binned_properties.meanAngDisplacement_pastPath.mean,...
    grouped_binned_properties(1).binned_properties.meanAngDisplacement_pastPath.sem,'lineprops','k'); hold on;
h.patch.FaceColor = colors_2(2,:);
h.mainLine.Color = colors(2,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
xlim([0 10])
xticks(0:2:10)
xtickangle(0)
ylim(ylim_deg)
yticks(ylim_deg_ticks);
ylabel('|Displacement|')

hold on
h = shadedErrorBar(bin_centers,grouped_binned_properties(2).binned_properties.meanAngDisplacement_futPath.mean,...
    grouped_binned_properties(1).binned_properties.meanAngDisplacement_futPath.sem,'lineprops','--k'); hold on;
h.patch.FaceColor = colors_2(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
h = shadedErrorBar(bin_centers,grouped_binned_properties(2).binned_properties.meanAngDisplacement_pastPath.mean,...
    grouped_binned_properties(2).binned_properties.meanAngDisplacement_pastPath.sem,'lineprops','--k'); hold on;
h.patch.FaceColor = colors_2(2,:);
h.mainLine.Color = colors(2,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
xlim([0 10])
xticks(0:2:10)
xtickangle(0)
ylim(ylim_deg)
yticks(ylim_deg_ticks);

ax1Chil = ax1.Children;
ax2 = nexttile(2);
copyobj(ax1Chil,ax2)
xlim([0 10])
xticks(0:2:10)
xtickangle(0)
ylim(ylim_deg)
yticks(ylim_deg_ticks);
xlabel('Time since arrival (s)')

% Difference of differences, real versus shuffle:
nexttile(3)
colors = [0 0 0; 0 0 0];
transparency_pcnt = 0.5;
colors_2 = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
    [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];

h = shadedErrorBar(bin_centers,grouped_binned_properties(1).binned_properties.future_minus_past_angle.mean,...
    grouped_binned_properties(1).binned_properties.past_minus_future_angle.sem,'lineprops','g'); hold on;
h.patch.FaceColor = colors_2(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
hold on;
h = shadedErrorBar(bin_centers,grouped_binned_properties(2).binned_properties.future_minus_past_angle.mean,...
    grouped_binned_properties(2).binned_properties.past_minus_future_angle.sem,'lineprops','--k'); hold on;
h.patch.FaceColor = colors_2(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
yline(0)
ylim(ylim_deg_diff_plot)

pvals = nan(length(bin_centers),1);
for i = 1:length(bin_centers)
    if strcmp(sig_test,'means')
        [h, pvals(i)] = ttest2(grouped_binned_properties(1).binned_properties.future_minus_past_angle.data{i},grouped_binned_properties(2).binned_properties.future_minus_past_angle.data{i});
    elseif strcmp(sig_test,'medians')
        [pvals(i)] = ranksum(grouped_binned_properties(1).binned_properties.future_minus_past_angle.data{i},grouped_binned_properties(2).binned_properties.future_minus_past_angle.data{i});
    end
end
hold on
ylimit = gca().YLim(2);
plot(bin_centers(pvals<0.05),ylimit*ones(sum(pvals<0.05),1),'.k')

xticks(0:2:10)
xtickangle(0)
xlim([0 10])
set(gcf, 'Color', 'white','PaperPositionMode','auto');
fig_title = [category1 '_v_' category2];

saveas(gcf,fullfile(fig_path,[fig_title '_v3_disp' num2str(replay_dispersionThr,2) '_rats_' num2str(rats)]),'jpg')
saveas(gcf,fullfile(fig_path,[fig_title '_v3_disp' num2str(replay_dispersionThr,2) '_rats_' num2str(rats)]),'pdf')

%%
group1_n_min = [min(grouped_binned_properties(1).binned_properties.meanAngDisplacement_futPath.n(bin_centers<=10))];
group1_n_max = [max(grouped_binned_properties(1).binned_properties.meanAngDisplacement_futPath.n(bin_centers<=10))];
group2_n_min = [min(grouped_binned_properties(2).binned_properties.meanAngDisplacement_futPath.n(bin_centers<=10))];
group2_n_max = [max(grouped_binned_properties(2).binned_properties.meanAngDisplacement_futPath.n(bin_centers<=10))];

group1_n_rf = [group1_n_min group1_n_max]
group2_n_rf = [group2_n_min group2_n_max]

group1_n_min = [min(grouped_binned_properties(1).binned_properties.meanAngDisplacement_pastPath.n(bin_centers<=10))];
group1_n_max = [max(grouped_binned_properties(1).binned_properties.meanAngDisplacement_pastPath.n(bin_centers<=10))];
group2_n_min = [min(grouped_binned_properties(2).binned_properties.meanAngDisplacement_pastPath.n(bin_centers<=10))];
group2_n_max = [max(grouped_binned_properties(2).binned_properties.meanAngDisplacement_pastPath.n(bin_centers<=10))];

group1_n_rp = [group1_n_min group1_n_max]
group2_n_rp = [group2_n_min group2_n_max]

group1_n_min = [min(grouped_binned_properties(1).binned_properties.future_minus_past_angle.n(bin_centers<=10))];
group1_n_max = [max(grouped_binned_properties(1).binned_properties.future_minus_past_angle.n(bin_centers<=10))];
group2_n_min = [min(grouped_binned_properties(2).binned_properties.future_minus_past_angle.n(bin_centers<=10))];
group2_n_max = [max(grouped_binned_properties(2).binned_properties.future_minus_past_angle.n(bin_centers<=10))];

group1_n_diff = [group1_n_min group1_n_max]
group2_n_diff = [group2_n_min group2_n_max]

replay_group1 = replay(inds1,:);
replay_group1 = replay_group1(replay_group1.time_since_real_drink_onset>0 & replay_group1.time_since_real_drink_onset<=10,:);
replay_group2 = replay(inds2,:);
replay_group2 = replay_group2(replay_group2.time_since_real_drink_onset>0 & replay_group2.time_since_real_drink_onset<=10,:);

replay_group1_n_total_rf = sum(~isnan(replay_group1.meanAngDisplacement_futPath))
replay_group1_n_total_rp = sum(~isnan(replay_group1.meanAngDisplacement_pastPath))
replay_group2_n_total_rf = sum(~isnan(replay_group2.meanAngDisplacement_futPath))
replay_group2_n_total_rp = sum(~isnan(replay_group2.meanAngDisplacement_pastPath))

replay_group1_n_total_nans_not_removed = height(replay_group1)
replay_group2_n_total_nans_not_removed = height(replay_group2)