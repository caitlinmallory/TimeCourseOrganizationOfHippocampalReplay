%% Plots the angle of all 'crossing points' from the past or future path, for 2 groups 
% Used in Fig S6
sig_test = 'medians';
ylim_deg = [40 120];
ylim_deg_ticks = ylim_deg(1):40:ylim_deg(2);
ylim_deg_diff_plot = [-60 40];

category1 = 'off';
category2 = 'on';
inds1 = find(replay.laser_state_binary==0);
inds2 = find(replay.laser_state_binary==1);

properties = {'angDisplacement_futPath'; 'angDisplacement_pastPath'};
windowSize = 2;
windowShift = 0.5;

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

        binned_data.data{time_bin} = abs(cell2mat(replay_sub.(properties{property})(replay_sub.time_into_stopping_period>= bin_edges(time_bin,1) & replay_sub.time_into_stopping_period<bin_edges(time_bin,2))));
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
    a = binned_properties_sub.angDisplacement_futPath(time_bin,:).data{:};
    b = binned_properties_sub.angDisplacement_pastPath(time_bin,:).data{:};
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

%% Version 1:
colors = [.4660 0.6740 0.1880; 0.4940 0.1840 0.5560;];
transparency_pcnt = 1;
colors_2 = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
    [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];

figure('Position',[1921 560 350 125])
tiledlayout(1,3,'TileSpacing','tight')
ax1=nexttile(1);
h = shadedErrorBar(bin_centers,grouped_binned_properties(1).binned_properties.angDisplacement_futPath.mean,...
    grouped_binned_properties(1).binned_properties.angDisplacement_futPath.sem,'lineprops','k'); hold on;
h.patch.FaceColor = colors_2(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none'; 
hold on
h = shadedErrorBar(bin_centers,grouped_binned_properties(2).binned_properties.angDisplacement_futPath.mean,...
    grouped_binned_properties(2).binned_properties.angDisplacement_futPath.sem,'lineprops','--k'); hold on;
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
    [h, pvals(i)] = ttest2(grouped_binned_properties(1).binned_properties.angDisplacement_futPath.data{i},grouped_binned_properties(2).binned_properties.angDisplacement_futPath.data{i});
    elseif strcmp(sig_test,'medians')
          [pvals(i)] = ranksum(grouped_binned_properties(1).binned_properties.angDisplacement_futPath.data{i},grouped_binned_properties(2).binned_properties.angDisplacement_futPath.data{i});
    end
end
hold on
ylimit = gca().YLim(2);
plot(bin_centers(pvals<0.05),ylimit*ones(sum(pvals<0.05)),'.k')

ax2=nexttile(2);
h = shadedErrorBar(bin_centers,grouped_binned_properties(1).binned_properties.angDisplacement_pastPath.mean,...
    grouped_binned_properties(1).binned_properties.angDisplacement_pastPath.sem,'lineprops','k'); hold on;
h.patch.FaceColor = colors_2(2,:);
h.mainLine.Color = colors(2,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none'; 
h = shadedErrorBar(bin_centers,grouped_binned_properties(2).binned_properties.angDisplacement_pastPath.mean,...
    grouped_binned_properties(2).binned_properties.angDisplacement_pastPath.sem,'lineprops','--k'); hold on;
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
    [h, pvals(i)] = ttest2(grouped_binned_properties(1).binned_properties.angDisplacement_pastPath.data{i},grouped_binned_properties(2).binned_properties.angDisplacement_pastPath.data{i});
    elseif strcmp(sig_test,'medians')
          [pvals(i)] = ranksum(grouped_binned_properties(1).binned_properties.angDisplacement_pastPath.data{i},grouped_binned_properties(2).binned_properties.angDisplacement_pastPath.data{i});
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
    grouped_binned_properties(1).binned_properties.future_minus_past_angle.sem,'lineprops','g'); hold on;
h.patch.FaceColor = colors_2(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none'; 
hold on;
h = shadedErrorBar(bin_centers,grouped_binned_properties(2).binned_properties.future_minus_past_angle.mean,...
    grouped_binned_properties(2).binned_properties.future_minus_past_angle.sem,'lineprops','--k'); hold on;
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


saveas(gcf,fullfile(fig_path,[fig_title '_disp' num2str(replay_dispersionThr,2) '_crossings_rats_' num2str(rats)]),'jpg')
saveas(gcf,fullfile(fig_path,[fig_title '_disp' num2str(replay_dispersionThr,2) '_crossings_rats_' num2str(rats)]),'pdf')

group1_n_min = [min(grouped_binned_properties(1).binned_properties.future_minus_past_angle.n)];
group1_n_max = [max(grouped_binned_properties(1).binned_properties.future_minus_past_angle.n)];
group2_n_min = [min(grouped_binned_properties(2).binned_properties.future_minus_past_angle.n)];
group2_n_max = [max(grouped_binned_properties(2).binned_properties.future_minus_past_angle.n)];

group1_n = [group1_n_min group1_n_max]
group2_n = [group2_n_min group2_n_max]

min(grouped_binned_properties(2).binned_properties.angDisplacement_futPath.n)

%% Version 2:
colors = [.4660 0.6740 0.1880; 0.4940 0.1840 0.5560;];
transparency_pcnt = 1;
colors_2 = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
    [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];

figure('Position',[1921 560 350 125])
tiledlayout(1,3,'TileSpacing','tight')
nexttile(1)
h = shadedErrorBar(bin_centers,grouped_binned_properties(1).binned_properties.angDisplacement_futPath.mean,...
    grouped_binned_properties(1).binned_properties.angDisplacement_futPath.sem,'lineprops','k'); hold on;
h.patch.FaceColor = colors_2(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none'; 
hold on
h = shadedErrorBar(bin_centers,grouped_binned_properties(1).binned_properties.angDisplacement_pastPath.mean,...
    grouped_binned_properties(1).binned_properties.angDisplacement_pastPath.sem,'lineprops','k'); hold on;
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
    [h, pvals(i)] = ttest2(grouped_binned_properties(1).binned_properties.angDisplacement_pastPath.data{i},grouped_binned_properties(1).binned_properties.angDisplacement_futPath.data{i})
    elseif strcmp(sig_test,'medians')
    [pvals(i)] = ranksum(grouped_binned_properties(1).binned_properties.angDisplacement_pastPath.data{i},grouped_binned_properties(1).binned_properties.angDisplacement_futPath.data{i})
    end
end
hold on
ylimit = gca().YLim(2);
plot(bin_centers(pvals<0.05),ylimit*ones(sum(pvals<0.05),1),'.k')

nexttile(2)
h = shadedErrorBar(bin_centers,grouped_binned_properties(2).binned_properties.angDisplacement_futPath.mean,...
    grouped_binned_properties(1).binned_properties.angDisplacement_futPath.sem,'lineprops','--k'); hold on;
h.patch.FaceColor = colors_2(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none'; 
h = shadedErrorBar(bin_centers,grouped_binned_properties(2).binned_properties.angDisplacement_pastPath.mean,...
    grouped_binned_properties(2).binned_properties.angDisplacement_pastPath.sem,'lineprops','--k'); hold on;
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
    [h, pvals(i)] = ttest2(grouped_binned_properties(2).binned_properties.angDisplacement_pastPath.data{i},grouped_binned_properties(2).binned_properties.angDisplacement_futPath.data{i})
    elseif strcmp(sig_test,'medians')
    [pvals(i)] = ranksum(grouped_binned_properties(2).binned_properties.angDisplacement_pastPath.data{i},grouped_binned_properties(2).binned_properties.angDisplacement_futPath.data{i})
    end
end
plot(bin_centers(pvals<0.05),ylimit*ones(sum(pvals<0.05),1),'.k')
% Difference of differences, real versus shuffle:


nexttile(3)
colors = [0 0 0; 0 0 0];
transparency_pcnt = 0.5;
colors_2 = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
    [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];

h = shadedErrorBar(bin_centers,grouped_binned_properties(1).binned_properties.future_minus_past_angle.mean,...
    grouped_binned_properties(1).binned_properties.future_minus_past_angle.sem,'lineprops','g'); hold on;
h.patch.FaceColor = colors_2(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none'; 
hold on;
h = shadedErrorBar(bin_centers,grouped_binned_properties(2).binned_properties.future_minus_past_angle.mean,...
    grouped_binned_properties(2).binned_properties.future_minus_past_angle.sem,'lineprops','--k'); hold on;
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

saveas(gcf,fullfile(fig_path,[fig_title '_v2_disp' num2str(replay_dispersionThr,2) '_crossings_rats_' num2str(rats)]),'jpg')
saveas(gcf,fullfile(fig_path,[fig_title '_v2_disp' num2str(replay_dispersionThr,2) '_crossings_rats_' num2str(rats)]),'pdf')

%% Version 3 (Overlay)

colors = [.4660 0.6740 0.1880; 0.4940 0.1840 0.5560;];
transparency_pcnt = 1;
colors_2 = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
    [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];

figure('Position',[1921 560 350 125])
tiledlayout(1,3,'TileSpacing','tight')
ax1 = nexttile(1);
h = shadedErrorBar(bin_centers,grouped_binned_properties(1).binned_properties.angDisplacement_futPath.mean,...
    grouped_binned_properties(1).binned_properties.angDisplacement_futPath.sem,'lineprops','k'); hold on;
h.patch.FaceColor = colors_2(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none'; 
hold on
h = shadedErrorBar(bin_centers,grouped_binned_properties(1).binned_properties.angDisplacement_pastPath.mean,...
    grouped_binned_properties(1).binned_properties.angDisplacement_pastPath.sem,'lineprops','k'); hold on;
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
h = shadedErrorBar(bin_centers,grouped_binned_properties(2).binned_properties.angDisplacement_futPath.mean,...
    grouped_binned_properties(1).binned_properties.angDisplacement_futPath.sem,'lineprops','--k'); hold on;
h.patch.FaceColor = colors_2(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none'; 
h = shadedErrorBar(bin_centers,grouped_binned_properties(2).binned_properties.angDisplacement_pastPath.mean,...
    grouped_binned_properties(2).binned_properties.angDisplacement_pastPath.sem,'lineprops','--k'); hold on;
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
    grouped_binned_properties(1).binned_properties.future_minus_past_angle.sem,'lineprops','g'); hold on;
h.patch.FaceColor = colors_2(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none'; 
hold on;
h = shadedErrorBar(bin_centers,grouped_binned_properties(2).binned_properties.future_minus_past_angle.mean,...
    grouped_binned_properties(2).binned_properties.future_minus_past_angle.sem,'lineprops','--k'); hold on;
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

saveas(gcf,fullfile(fig_path,[fig_title '_v3_disp' num2str(replay_dispersionThr,2) '_crossings_rats_' num2str(rats)]),'jpg')
saveas(gcf,fullfile(fig_path,[fig_title '_v3_disp' num2str(replay_dispersionThr,2) '_crossings_rats_' num2str(rats)]),'pdf')

