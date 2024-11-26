% ylim_rates = [0 0.35];
% ylim_rate_diffs = [-0.2 0.3];
correct_for_multiple_comparisons = 0;
alpha = 0.05;
if align_to_drink_onset==1
   num_comparisons = sum(bin_centers_plot>0 & bin_centers_plot<10);
elseif align_to_drink_offset==1
   num_comparisons = sum(bin_centers_plot>-10 & bin_centers_plot<0);
end

if correct_for_multiple_comparisons
    alpha = alpha/num_comparisons;
end

if home_trials_only==1
    title_ending = 'home_trials';
elseif away_trials_only==1
    title_ending = 'away_trials';
else
    title_ending = 'all_trials';
end

% Main Fig 4:
ylim_rates = [0 0.4];
ylim_rate_diffs = [-0.2 0.4];

% Main Fig 2:
% ylim_rates = [0 0.2];
% ylim_rate_diffs = [-0.2 0.2];

% Supplementary Figure (overlapping with SWRS):
% ylim_rates = [0 0.1];
% ylim_rate_diffs = [-0.05 0.1];


% ylim_rates = [0 0.3];
% ylim_rate_diffs = [-0.2 0.3];


transparency_pcnt = 1;
colors_2 = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
    [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];


data1_n = [min(sum(~isnan(combined_hists_rates.more_future(inds1,:)))) max(sum(~isnan(combined_hists_rates.more_future(inds1,:))))]
data2_n = [min(sum(~isnan(combined_hists_rates.more_future(inds2,:)))) max(sum(~isnan(combined_hists_rates.more_future(inds2,:))))]


figure('Position',[1921 560 350 125])
tiledlayout(1,3,'TileSpacing','tight')
% Plot home prospective versus away prospective
ax1 = nexttile(1);

data = combined_hists_rates.more_future(inds1,:);
x_mean = nanmean(data);
x_err = nanstd(data)./sqrt(sum(~isnan(data)));
if align_to_drink_offset==1
    if smooth_rate_plots==1
h = shadedErrorBar(bin_centers_plot((filter_length-1)/2:end),x_mean((filter_length-1)/2:end),...
    x_err((filter_length-1)/2:end),'lineprops','k'); hold on;
    else
         h = shadedErrorBar(bin_centers_plot,x_mean,...
    x_err,'lineprops','k'); hold on;
    end
else
 h = shadedErrorBar(bin_centers_plot,x_mean,...
    x_err,'lineprops','k'); hold on;
end
h.patch.FaceColor = colors(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.patch.FaceColor = colors_2(1,:);
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';
hold on
data = combined_hists_rates.more_future(inds2,:);
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
h.patch.FaceColor = colors(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.patch.FaceColor = colors_2(1,:);
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';
% xlabel('Time since arrival (s)')
ylabel('Events/s')
%legend({'Future, off', 'Past, off', 'Future, on', 'Past, on'})
ylim(ylim_rates)
yticks(ylim_rates(1):0.1:ylim_rates(2))
xticks(0:2:10)
xtickangle(0);
if strcmp(sig_test,'means')
    for i = 1:size(combined_hists_rates.more_future,2)
        a = combined_hists_rates.more_future(inds1,i);
        b = combined_hists_rates.more_future(inds2,i);
        [~,pvals(i)] = ttest2(a,b);
    end
    sig_times = pvals < alpha;
elseif strcmp(sig_test,'medians')
    for i = 1:size(combined_hists_rates.more_future,2)
        a = combined_hists_rates.more_future(inds1,i);
        b = combined_hists_rates.more_future(inds2,i);
        [pvals(i)] = ranksum(a,b);
    end
        sig_times = pvals < alpha;
end
sig_times = logical(sig_times);
plot(bin_centers_plot(sig_times),ylimit(2)*ones(sum(sig_times),1),'.k')

% Plot home retrospective versus away retrospective
ax2 = nexttile(2);
data = combined_hists_rates.more_past(inds1,:);
x_mean = nanmean(data);
x_err = nanstd(data)./sqrt(sum(~isnan(data)));
if align_to_drink_offset==1
    if smooth_rate_plots==1
h = shadedErrorBar(bin_centers_plot((filter_length-1)/2:end),x_mean((filter_length-1)/2:end),...
    x_err((filter_length-1)/2:end),'lineprops','-k'); hold on;
    else
         h = shadedErrorBar(bin_centers_plot,x_mean,...
    x_err,'lineprops','-m'); hold on;
    end
    xlim([-10 0]);
else
 h = shadedErrorBar(bin_centers_plot,x_mean,...
    x_err,'lineprops','-m'); hold on;
    xlim([0 10]);
end
h.patch.FaceColor = colors(2,:);
h.mainLine.Color = colors(2,:);
h.mainLine.LineWidth = 1;
h.patch.FaceColor = colors_2(2,:);
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';
hold on
data = combined_hists_rates.more_past(inds2,:);
x_mean = nanmean(data);
x_err = nanstd(data)./sqrt(sum(~isnan(data)));
if align_to_drink_offset==1
    if smooth_rate_plots==1
h = shadedErrorBar(bin_centers_plot((filter_length-1)/2:end),x_mean((filter_length-1)/2:end),...
    x_err((filter_length-1)/2:end),'lineprops','--k'); hold on;
    else
         h = shadedErrorBar(bin_centers_plot,x_mean,...
    x_err,'lineprops','--m'); hold on;
    end
    xlim([-10 0]);
else
 h = shadedErrorBar(bin_centers_plot,x_mean,...
    x_err,'lineprops','--k'); hold on;
    xlim([0 10]);
end
h.patch.FaceColor = colors(2,:);
h.mainLine.Color = colors(2,:);
h.mainLine.LineWidth = 1;
h.patch.FaceColor = colors_2(2,:);
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';
xlabel('Time since arrival (s)')
% ylabel('Events/s')
%legend({'Future, off', 'Past, off', 'Future, on', 'Past, on'})
ylim(ylim_rates)
yticks(ylim_rates(1):0.1:ylim_rates(2))
xtickangle(0);
xticks(0:2:10)
if strcmp(sig_test,'means')
    for i = 1:size(combined_hists_rates.more_future,2)
        a = combined_hists_rates.more_past(inds1,i);
        b = combined_hists_rates.more_past(inds2,i);
        [~,pvals(i)] = ttest2(a,b);
    end
    sig_times = pvals < alpha;
elseif strcmp(sig_test,'medians')
    for i = 1:size(combined_hists_rates.more_future,2)
        a = combined_hists_rates.more_past(inds1,i);
        b = combined_hists_rates.more_past(inds2,i);
        [pvals(i)] = ranksum(a,b);
    end
        sig_times = pvals < alpha;
end
sig_times = logical(sig_times);
plot(bin_centers_plot(sig_times),ylimit(2)*ones(sum(sig_times),1),'.k')

nexttile(3)
hold on
data1 = combined_hists_rates.more_future(inds1,:) - combined_hists_rates.more_past(inds1,:);
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
data2 = combined_hists_rates.more_future(inds2,:) - combined_hists_rates.more_past(inds2,:);
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
% xlabel('Time since arrival (s)')
% ylabel('Events/s')
yline(0)
 box off
ylim(ylim_rate_diffs)
yticks(ylim_rate_diffs(1):0.1:ylim_rate_diffs(2))
xtickangle(0);
xticks(0:2:10)

diff_of_diffs = abs(nanmean(data1)-nanmean(data2));
pvals = nan(size(combined_hists_rates.more_future,2),1);

if strcmp(sig_test,'shuffle') == 1
    shuffle_diff_of_diffs = nan(1000,size(diff_of_diffs,2));
    for i = 1:1000
        shuffled_inds1 = randperm(size(combined_hists_rates.more_future,1), size(inds1,1));
        shuffled_inds2 = setdiff(1:size(combined_hists_rates.more_future,1),shuffled_inds1);
        shuffle1_diff = combined_hists_rates.more_future(shuffled_inds1,:)-combined_hists_rates.more_past(shuffled_inds1,:);
        shuffle2_diff = combined_hists_rates.more_future(shuffled_inds2,:)-combined_hists_rates.more_past(shuffled_inds2,:);
        shuffle_diff_of_diffs(i,:) = abs(nanmean(shuffle1_diff)-nanmean(shuffle2_diff));
    end
    sig_times = diff_of_diffs>quantile(shuffle_diff_of_diffs,1-alpha/2);
elseif strcmp(sig_test,'means')
    for i = 1:size(combined_hists_rates.more_future,2)
        a = combined_hists_rates.more_future(inds1,i)-combined_hists_rates.more_past(inds1,i);
        b = combined_hists_rates.more_future(inds2,i)-combined_hists_rates.more_past(inds2,i);
        [~,pvals(i)] = ttest2(a,b);
    end
    sig_times = pvals < alpha;
elseif strcmp(sig_test,'medians')
    for i = 1:size(combined_hists_rates.more_future,2)
        a = combined_hists_rates.more_future(inds1,i)-combined_hists_rates.more_past(inds1,i);
        b = combined_hists_rates.more_future(inds2,i)-combined_hists_rates.more_past(inds2,i);
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
saveas(gcf,fullfile(fig_path,['replay_rates_' category1 '_versus_' category2 '_DegThr_' num2str(thrForCategorization_include) 'sdeThr_' num2str(sd_thr) '_disp' num2str(replay_dispersionThr) '_' title_ending]),'jpg')
saveas(gcf,fullfile(fig_path,['replay_rates_' category1 '_versus_' category2 '_DegThr_' num2str(thrForCategorization_include) 'sdeThr_' num2str(sd_thr) '_disp' num2str(replay_dispersionThr) '_' title_ending]),'pdf')


%% Version 2:

num_bins = size(combined_hists_rates.more_future,2);
figure('Position',[1921 560 350 125])
tiledlayout(1,3,'TileSpacing','tight')
% Plot home prospective versus away prospective
ax1 = nexttile(1);
data = combined_hists_rates.more_future(inds1,:);
x_mean = nanmean(data);
x_err = nanstd(data)./sqrt(sum(~isnan(data)));
if align_to_drink_offset==1
    if smooth_rate_plots==1
h = shadedErrorBar(bin_centers_plot((filter_length-1)/2:end),x_mean((filter_length-1)/2:end),...
    x_err((filter_length-1)/2:end),'lineprops','-k'); hold on;
    else
         h = shadedErrorBar(bin_centers_plot,x_mean,...
    x_err,'lineprops','-m'); hold on;
    end
    xlim([-10 0]);
else
 h = shadedErrorBar(bin_centers_plot,x_mean,...
    x_err,'lineprops','-k'); hold on;
    xlim([0 10]);
end
h.patch.FaceColor = colors(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.patch.FaceColor = colors_2(1,:);
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';
hold on
data = combined_hists_rates.more_past(inds1,:);
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
h.patch.FaceColor = colors(2,:);
h.mainLine.Color = colors(2,:);
h.mainLine.LineWidth = 1;
h.patch.FaceColor = colors_2(2,:);
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';
% xlabel('Time since arrival (s)')
ylabel('Events/s')
%legend({'Future, off', 'Past, off', 'Future, on', 'Past, on'})
ylim(ylim_rates)
yticks(ylim_rates(1):0.1:ylim_rates(2))
xtickangle(0);
xticks(0:2:10)

if strcmp(sig_test,'shuffle') == 1
%     shuffle_diff_of_diffs = nan(1000,size(diff_of_diffs,2));
%     for i = 1:1000
%         shuffled_inds1 = randperm(size(combined_hists_rates.more_future,1), size(inds1,1));
%         shuffled_inds2 = setdiff(1:size(combined_hists_rates.more_future,1),shuffled_inds1);
%         shuffle1_diff = combined_hists_rates.more_future(shuffled_inds1,:)-combined_hists_rates.more_past(shuffled_inds1,:);
%         shuffle2_diff = combined_hists_rates.more_future(shuffled_inds2,:)-combined_hists_rates.more_past(shuffled_inds2,:);
%         shuffle_diff_of_diffs(i,:) = abs(nanmean(shuffle1_diff)-nanmean(shuffle2_diff));
%     end
%     sig_times = diff_of_diffs>quantile(shuffle_diff_of_diffs,0.975);
elseif strcmp(sig_test,'means')
    for i = 1:size(combined_hists_rates.more_future,2)
        a = combined_hists_rates.more_future(inds1,i);
        b = combined_hists_rates.more_past(inds1,i);
        [~,pvals(i)] = ttest2(a,b);
    end
    sig_times = pvals < alpha;
elseif strcmp(sig_test,'medians')
    for i = 1:size(combined_hists_rates.more_future,2)
        a = combined_hists_rates.more_future(inds1,i);
        b = combined_hists_rates.more_past(inds1,i);
        [pvals(i)] = ranksum(a,b);
    end
    sig_times = pvals < alpha;
end
plot(bin_centers_plot(sig_times),ylimit(2)*ones(sum(sig_times),1),'.k')

% Plot home retrospective versus away retrospective
ax2 = nexttile(2);
data = combined_hists_rates.more_future(inds2,:);
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
h.patch.FaceColor = colors(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.patch.FaceColor = colors_2(1,:);
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';
hold on
data = combined_hists_rates.more_past(inds2,:);
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
h.patch.FaceColor = colors(2,:);
h.mainLine.Color = colors(2,:);
h.mainLine.LineWidth = 1;
h.patch.FaceColor = colors_2(2,:);
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';
xlabel('Time since arrival (s)')
% ylabel('Events/s')
%legend({'Future, off', 'Past, off', 'Future, on', 'Past, on'})
ylim(ylim_rates)
yticks(ylim_rates(1):0.1:ylim_rates(2))
xtickangle(0);
xticks(0:2:10)


if strcmp(sig_test,'shuffle') == 1
%     shuffle_diff_of_diffs = nan(1000,size(diff_of_diffs,2));
%     for i = 1:1000
%         shuffled_inds1 = randperm(size(combined_hists_rates.more_future,1), size(inds1,1));
%         shuffled_inds2 = setdiff(1:size(combined_hists_rates.more_future,1),shuffled_inds1);
%         shuffle1_diff = combined_hists_rates.more_future(shuffled_inds1,:)-combined_hists_rates.more_past(shuffled_inds1,:);
%         shuffle2_diff = combined_hists_rates.more_future(shuffled_inds2,:)-combined_hists_rates.more_past(shuffled_inds2,:);
%         shuffle_diff_of_diffs(i,:) = abs(nanmean(shuffle1_diff)-nanmean(shuffle2_diff));
%     end
%     sig_times = diff_of_diffs>quantile(shuffle_diff_of_diffs,0.975);
elseif strcmp(sig_test,'means')
    for i = 1:size(combined_hists_rates.more_future,2)
        a = combined_hists_rates.more_future(inds2,i);
        b = combined_hists_rates.more_past(inds2,i);
        [~,pvals(i)] = ttest2(a,b);
    end
    sig_times = pvals < alpha;
elseif strcmp(sig_test,'medians')
    for i = 1:size(combined_hists_rates.more_future,2)
        a = combined_hists_rates.more_future(inds2,i);
        b = combined_hists_rates.more_past(inds2,i);
        [pvals(i)] = ranksum(a,b);
    end
    sig_times = pvals < alpha;
end
plot(bin_centers_plot(sig_times),ylimit(2)*ones(sum(sig_times),1),'.k')


nexttile(3)
hold on
data1 = combined_hists_rates.more_future(inds1,:) - combined_hists_rates.more_past(inds1,:);
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
data2 = combined_hists_rates.more_future(inds2,:) - combined_hists_rates.more_past(inds2,:);
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

diff_of_diffs = abs(nanmean(data1)-nanmean(data2));
pvals = nan(size(combined_hists_rates.more_future,2),1);



if strcmp(sig_test,'shuffle') == 1
%     shuffle_diff_of_diffs = nan(1000,size(diff_of_diffs,2));
%     for i = 1:1000
%         shuffled_inds1 = randperm(size(combined_hists_rates.more_future,1), size(inds1,1));
%         shuffled_inds2 = setdiff(1:size(combined_hists_rates.more_future,1),shuffled_inds1);
%         shuffle1_diff = combined_hists_rates.more_future(shuffled_inds1,:)-combined_hists_rates.more_past(shuffled_inds1,:);
%         shuffle2_diff = combined_hists_rates.more_future(shuffled_inds2,:)-combined_hists_rates.more_past(shuffled_inds2,:);
%         shuffle_diff_of_diffs(i,:) = abs(nanmean(shuffle1_diff)-nanmean(shuffle2_diff));
%     end
%     sig_times = diff_of_diffs>quantile(shuffle_diff_of_diffs,0.975);
elseif strcmp(sig_test,'means')
    for i = 1:size(combined_hists_rates.more_future,2)
        a = combined_hists_rates.more_future(inds1,i)-combined_hists_rates.more_past(inds1,i);
        b = combined_hists_rates.more_future(inds2,i)-combined_hists_rates.more_past(inds2,i);
        [~,pvals(i)] = ttest2(a,b);
    end
    sig_times = pvals < alpha;
elseif strcmp(sig_test,'medians')
    for i = 1:size(combined_hists_rates.more_future,2)
        a = combined_hists_rates.more_future(inds1,i)-combined_hists_rates.more_past(inds1,i);
        b = combined_hists_rates.more_future(inds2,i)-combined_hists_rates.more_past(inds2,i);
        [pvals(i)] = ranksum(a,b);
    end
    sig_times = pvals < alpha;
end




hold on
ylim(ylim_rate_diffs)
yticks(ylim_rate_diffs(1):0.1:ylim_rate_diffs(2))
ylimit = gca().YLim;
plot(bin_centers_plot(sig_times),ylimit(2)*ones(sum(sig_times),1),'.k')

set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
saveas(gcf,fullfile(fig_path,['replay_rates_v2_' category1 '_versus_' category2 '_DegThr_' num2str(thrForCategorization_include) 'sdeThr_' num2str(sd_thr) '_disp' num2str(replay_dispersionThr) '_' title_ending]),'jpg')
saveas(gcf,fullfile(fig_path,['replay_rates_v2_' category1 '_versus_' category2 '_DegThr_' num2str(thrForCategorization_include) 'sdeThr_' num2str(sd_thr) '_disp' num2str(replay_dispersionThr) '_' title_ending]),'pdf')


%%
%% Version 3:

num_bins = size(combined_hists_rates.more_future,2);
figure('Position',[1921 560 350 125])
tiledlayout(1,3,'TileSpacing','tight')
% Plot home prospective versus away prospective
ax1 = nexttile(1);
data = combined_hists_rates.more_future(inds1,:);
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
h.patch.FaceColor = colors(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.patch.FaceColor = colors_2(1,:);
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';
hold on
data = combined_hists_rates.more_past(inds1,:);
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
h.patch.FaceColor = colors(2,:);
h.mainLine.Color = colors(2,:);
h.mainLine.LineWidth = 1;
h.patch.FaceColor = colors_2(2,:);
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';

% Plot home retrospective versus away retrospective
hold on
data = combined_hists_rates.more_future(inds2,:);
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
h.patch.FaceColor = colors(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.patch.FaceColor = colors_2(1,:);
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';
hold on
data = combined_hists_rates.more_past(inds2,:);
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
hold on
data1 = combined_hists_rates.more_future(inds1,:) - combined_hists_rates.more_past(inds1,:);
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
data2 = combined_hists_rates.more_future(inds2,:) - combined_hists_rates.more_past(inds2,:);
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
% xlabel('Time since arrival (s)')
% ylabel('Events/s')
yline(0)
 box off
ylim([-0.2 0.2])
yticks(ylim_rates(1):0.1:ylim_rates(2))
xtickangle(0);
xticks(0:2:10)

diff_of_diffs = abs(nanmean(data1)-nanmean(data2));
pvals = nan(size(combined_hists_rates.more_future,2),1);

if strcmp(sig_test,'shuffle') == 1
    shuffle_diff_of_diffs = nan(1000,size(diff_of_diffs,2));
    for i = 1:1000
        shuffled_inds1 = randperm(size(combined_hists_rates.more_future,1), size(inds1,1));
        shuffled_inds2 = setdiff(1:size(combined_hists_rates.more_future,1),shuffled_inds1);
        shuffle1_diff = combined_hists_rates.more_future(shuffled_inds1,:)-combined_hists_rates.more_past(shuffled_inds1,:);
        shuffle2_diff = combined_hists_rates.more_future(shuffled_inds2,:)-combined_hists_rates.more_past(shuffled_inds2,:);
        shuffle_diff_of_diffs(i,:) = abs(nanmean(shuffle1_diff)-nanmean(shuffle2_diff));
    end
    sig_times = diff_of_diffs>quantile(shuffle_diff_of_diffs,0.975);
elseif strcmp(sig_test,'means')
    for i = 1:size(combined_hists_rates.more_future,2)
        a = combined_hists_rates.more_future(inds1,i)-combined_hists_rates.more_past(inds1,i);
        b = combined_hists_rates.more_future(inds2,i)-combined_hists_rates.more_past(inds2,i);
        [~,pvals(i)] = ttest2(a,b);
    end
    sig_times = pvals < alpha;
elseif strcmp(sig_test,'medians')
    for i = 1:size(combined_hists_rates.more_future,2)
        a = combined_hists_rates.more_future(inds1,i)-combined_hists_rates.more_past(inds1,i);
        b = combined_hists_rates.more_future(inds2,i)-combined_hists_rates.more_past(inds2,i);
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

%%
data1_n = [min(sum(~isnan(combined_hists_rates.more_future(inds1,:)))) max(sum(~isnan(combined_hists_rates.more_future(inds1,:))))]
data2_n = [min(sum(~isnan(combined_hists_rates.more_future(inds2,:)))) max(sum(~isnan(combined_hists_rates.more_future(inds2,:))))]


figure('Position',[1921 560 350 125])
tiledlayout(1,3,'TileSpacing','tight')
% Plot home prospective versus away prospective
ax1 = nexttile(1);

data = combined_hists_rates.either_future_or_past(inds1,:);
x_mean = nanmean(data);
x_err = nanstd(data)./sqrt(sum(~isnan(data)));
if align_to_drink_offset==1
    if smooth_rate_plots==1
h = shadedErrorBar(bin_centers_plot((filter_length-1)/2:end),x_mean((filter_length-1)/2:end),...
    x_err((filter_length-1)/2:end),'lineprops','k'); hold on;
    else
         h = shadedErrorBar(bin_centers_plot,x_mean,...
    x_err,'lineprops','k'); hold on;
    end
else
 h = shadedErrorBar(bin_centers_plot,x_mean,...
    x_err,'lineprops','k'); hold on;
end
h.patch.FaceColor = [0 0 0];
h.mainLine.Color = [0 0 0];
h.mainLine.LineWidth = 1;
h.patch.FaceColor = [0 0 0];
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';
hold on

data = combined_hists_rates.more_future(inds2,:);
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
h.patch.FaceColor = colors(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.patch.FaceColor = colors_2(1,:);
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';
% xlabel('Time since arrival (s)')
ylabel('Events/s')
%legend({'Future, off', 'Past, off', 'Future, on', 'Past, on'})
ylim(ylim_rates)
yticks(ylim_rates(1):0.1:ylim_rates(2))
xticks(0:2:10)
xtickangle(0);

sig_times = logical(sig_times);
plot(bin_centers_plot(sig_times),ylimit(2)*ones(sum(sig_times),1),'.k')

% Plot home retrospective versus away retrospective

data = combined_hists_rates.more_past(inds2,:);
x_mean = nanmean(data);
x_err = nanstd(data)./sqrt(sum(~isnan(data)));
if align_to_drink_offset==1
    if smooth_rate_plots==1
h = shadedErrorBar(bin_centers_plot((filter_length-1)/2:end),x_mean((filter_length-1)/2:end),...
    x_err((filter_length-1)/2:end),'lineprops','--k'); hold on;
    else
         h = shadedErrorBar(bin_centers_plot,x_mean,...
    x_err,'lineprops','--m'); hold on;
    end
    xlim([-10 0]);
else
 h = shadedErrorBar(bin_centers_plot,x_mean,...
    x_err,'lineprops','--k'); hold on;
    xlim([0 10]);
end
h.patch.FaceColor = colors(2,:);
h.mainLine.Color = colors(2,:);
h.mainLine.LineWidth = 1;
h.patch.FaceColor = colors_2(2,:);
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';
xlabel('Time since arrival (s)')
% ylabel('Events/s')
%legend({'Future, off', 'Past, off', 'Future, on', 'Past, on'})
ylim(ylim_rates)
yticks(ylim_rates(1):0.1:ylim_rates(2))
xtickangle(0);
xticks(0:2:10)
if strcmp(sig_test,'means')
    for i = 1:size(combined_hists_rates.more_future,2)
        a = combined_hists_rates.more_past(inds1,i);
        b = combined_hists_rates.more_past(inds2,i);
        [~,pvals(i)] = ttest2(a,b);
    end
    sig_times = pvals < alpha;
elseif strcmp(sig_test,'medians')
    for i = 1:size(combined_hists_rates.more_future,2)
        a = combined_hists_rates.more_past(inds1,i);
        b = combined_hists_rates.more_past(inds2,i);
        [pvals(i)] = ranksum(a,b);
    end
        sig_times = pvals < alpha;
end
sig_times = logical(sig_times);
plot(bin_centers_plot(sig_times),ylimit(2)*ones(sum(sig_times),1),'.k')

set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
saveas(gcf,fullfile(fig_path,['replay_rates_' category1 '_versus_' category2 '_DegThr_' num2str(thrForCategorization_include) 'sdeThr_' num2str(sd_thr) '_disp' num2str(replay_dispersionThr) '_' title_ending]),'jpg')
saveas(gcf,fullfile(fig_path,['replay_rates_' category1 '_versus_' category2 '_DegThr_' num2str(thrForCategorization_include) 'sdeThr_' num2str(sd_thr) '_disp' num2str(replay_dispersionThr) '_' title_ending]),'pdf')

%%
figure('Position',[1921 560 350 125])
tiledlayout(1,3,'TileSpacing','tight')
% Plot home prospective versus away prospective
ax1 = nexttile(1);

data = combined_hists_rates.either_future_or_past(inds1,:);
x_mean = nanmean(data);
x_err = nanstd(data)./sqrt(sum(~isnan(data)));
if align_to_drink_offset==1
    if smooth_rate_plots==1
h = shadedErrorBar(bin_centers_plot((filter_length-1)/2:end),x_mean((filter_length-1)/2:end),...
    x_err((filter_length-1)/2:end),'lineprops','k'); hold on;
    else
         h = shadedErrorBar(bin_centers_plot,x_mean,...
    x_err,'lineprops','k'); hold on;
    end
else
 h = shadedErrorBar(bin_centers_plot,x_mean,...
    x_err,'lineprops','k'); hold on;
end
h.patch.FaceColor = [0 0 0];
h.mainLine.Color = [0 0 0];
h.mainLine.LineWidth = 1;
h.patch.FaceColor = [0 0 0];
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';
hold on


data = combined_hists_rates.either_future_or_past(inds2,:);
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
else
 h = shadedErrorBar(bin_centers_plot,x_mean,...
    x_err,'lineprops','--k'); hold on;
end
h.patch.FaceColor = [0 0 0];
h.mainLine.Color = [0 0 0];
h.mainLine.LineWidth = 1;
h.patch.FaceColor = [0 0 0];
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';
hold on

