% ylim_rates = [0 0.35];
% ylim_rate_diffs = [-0.2 0.3];

% Main Fig 4:
% ylim_rates = [0 0.4];
% ylim_rate_diffs = [-0.2 0.3];

% Main Fig 2:
% ylim_rates = [0 0.2];
% ylim_rate_diffs = [-0.2 0.2];

% Supplementary Figure (overlapping with SWRS):
% ylim_rates = [0 0.1];
% ylim_rate_diffs = [-0.05 0.1];

ylim_rates = [0 1.2];
data1_n = [min(sum(~isnan(combined_hists_rates.replay(inds1,:)))) max(sum(~isnan(combined_hists_rates.replay(inds1,:))))]
data2_n = [min(sum(~isnan(combined_hists_rates.replay(inds2,:)))) max(sum(~isnan(combined_hists_rates.replay(inds2,:))))]


figure('Position',[1921 560 350 125])
tiledlayout(1,3,'TileSpacing','tight')
% Plot home prospective versus away prospective
ax1 = nexttile(1);

data = combined_hists_rates.replay(inds1,:);
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
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';
hold on

data = combined_hists_rates.replay(inds2,:);
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

h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';
% xlabel('Time since arrival (s)')
ylabel('Events/s')
%legend({'Future, off', 'Past, off', 'Future, on', 'Past, on'})
ylim(ylim_rates)
yticks(ylim_rates(1):0.4:ylim_rates(2))
xtickangle(0);

data1 = combined_hists_rates.replay(inds1);
data2 = combined_hists_rates.replay(inds2);
diff_of_diffs = abs(nanmean(data1)-nanmean(data2));
sig_times = nan(size(combined_hists_rates.replay,2),1);

if strcmp(sig_test,'shuffle') == 1
    shuffle_diffs = nan(1000,size(diff_of_diffs,2));
    % not right- need to fix this;
%     for i = 1:1000
%         shuffled_inds1 = randperm(size(combined_hists_rates.replay,1), size(inds1,1));
%         shuffled_inds2 = setdiff(1:size(combined_hists_rates.replay,1),shuffled_inds1);
%         shuffle_diffs(i,:) = abs(combined_hists_rates.replay(shuffled_inds1)-combined_hists_rates.replay(shuffled_inds1));
%     end
%     sig_times = diff_of_diffs>quantile(shuffle_diffs,0.975);
elseif strcmp(sig_test,'means')
    for i = 1:size(combined_hists_rates.replay,2)
        a = combined_hists_rates.replay(inds1,i);
        b = combined_hists_rates.replay(inds2,i);
        sig_times(i) = ttest2(a,b);
    end
elseif strcmp(sig_test,'medians')
    for i = 1:size(combined_hists_rates.future,2)
        a = combined_hists_rates.replay(inds1,i);
        b = combined_hists_rates.replay(inds2,i);
        [~,sig_times(i)] = ranksum(a,b);
    end
end
sig_times = logical(sig_times);


hold on
ylimit = gca().YLim;
plot(bin_centers_plot(sig_times),ylimit(2)*ones(sum(sig_times),1),'.k')

set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
saveas(gcf,fullfile(fig_path,['overall_replay_rates_' category1 '_versus_' category2 '_DegThr_' num2str(thrForCategorization_include) 'sdeThr_' num2str(sd_thr) '_disp' num2str(replay_dispersionThr) '_' title_ending]),'jpg')
saveas(gcf,fullfile(fig_path,['overall_replay_rates_' category1 '_versus_' category2 '_DegThr_' num2str(thrForCategorization_include) 'sdeThr_' num2str(sd_thr) '_disp' num2str(replay_dispersionThr) '_' title_ending]),'pdf')

