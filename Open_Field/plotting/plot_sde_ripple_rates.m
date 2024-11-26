
figure('Position',[1921 560 350 125])
tiledlayout(1,3,'TileSpacing','tight')

% Plot SDE and Ripple rate on laser off or laser on trials
ax1 = nexttile(1);
inds_to_keep = find(t_sde.trial_duration >= duration_thr & t_sde.laser_state == 0);
x_mean  = nanmean(combined_hists_rates.sde(inds_to_keep,1:num_bins_to_plot));
x_err = nanstd(combined_hists_rates.sde(inds_to_keep,1:num_bins_to_plot)./sqrt((sum(~isnan(combined_hists_rates.sde(inds_to_keep,1:num_bins_to_plot))))));

h = shadedErrorBar(bin_centers,x_mean,...
    x_err,'lineprops','k'); hold on;
h.mainLine.Color = color(1,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
hold on

ylim([0 0.8])
xlabel('Time since arrival (s)')
ylabel('Events/s')


%%
ax2 = nexttile(2);
inds_to_keep = find(t_ripple.trial_duration >= duration_thr & t_ripple.laser_state == 0);
x_mean  = nanmean(combined_hists_rates.ripple(inds_to_keep,1:num_bins_to_plot));
x_err = nanstd(combined_hists_rates.ripple(inds_to_keep,1:num_bins_to_plot)./sqrt((sum(~isnan(combined_hists_rates.ripple(inds_to_keep,1:num_bins_to_plot))))));

h = shadedErrorBar(bin_centers,x_mean,...
    x_err,'lineprops','k'); hold on;
h.mainLine.Color = color(1,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
hold on

inds_to_keep = find(t_ripple.trial_duration >= duration_thr & t_ripple.laser_state == 1);
x_mean  = nanmean(combined_hists_rates.ripple(inds_to_keep,1:num_bins_to_plot));
x_err = nanstd(combined_hists_rates.ripple(inds_to_keep,1:num_bins_to_plot)./sqrt((sum(~isnan(combined_hists_rates.ripple(inds_to_keep,1:num_bins_to_plot))))));

h = shadedErrorBar(bin_centers,x_mean,...
    x_err,'lineprops','--k'); hold on;
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
hold on
ylim([0 0.8])
xticks(0:2:10)
xtickangle(0);
% title('All trials')
% xlabel('Time since arrival (s)')
% ylabel('Events/s')
yticks([0 0.4 0.8])

