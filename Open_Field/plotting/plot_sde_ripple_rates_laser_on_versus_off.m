
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

inds_to_keep = find(t_sde.trial_duration >= duration_thr & t_sde.laser_state == 1);
x_mean  = nanmean(combined_hists_rates.sde(inds_to_keep,1:num_bins_to_plot));
x_err = nanstd(combined_hists_rates.sde(inds_to_keep,1:num_bins_to_plot)./sqrt((sum(~isnan(combined_hists_rates.sde(inds_to_keep,1:num_bins_to_plot))))));

h = shadedErrorBar(bin_centers,x_mean,...
    x_err,'lineprops','--k'); hold on;
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
hold on
ylim([0 0.8])
xticks(0:2:10)
xtickangle(0);
yticks([0 0.4 0.8])
% title('All trials')
xlabel('Time since arrival (s)')
ylabel('Events/s')

sig_times_sde = nan(num_bins_to_plot,1);
for i = 1:num_bins_to_plot
    inds1 = find(t_sde.trial_duration >= duration_thr & t_sde.laser_state == 0);
    inds2 = find(t_sde.trial_duration >= duration_thr & t_sde.laser_state == 1);
    data1 = combined_hists_rates.sde(inds1,i);
    data2 = combined_hists_rates.sde(inds2,i);

    [~,sig_times_sde(i)] = ranksum(data1,data2);
end

sig_times_sde = logical(sig_times_sde);
hold on
ylimit = gca().YLim;
plot(bin_centers(sig_times_sde),ylimit(2)*ones(sum(sig_times_sde),1),'.k')

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
sig_times_ripple = nan(num_bins_to_plot,1);
for i = 1:num_bins_to_plot
    inds1 = find(t_ripple.trial_duration >= duration_thr & t_ripple.laser_state == 0);
    inds2 = find(t_ripple.trial_duration >= duration_thr & t_ripple.laser_state == 1);
    data1 = combined_hists_rates.ripple(inds1,i);
    data2 = combined_hists_rates.ripple(inds2,i);

    [~,sig_times_ripple(i)] = ranksum(data1,data2);
end

sig_times_ripple = logical(sig_times_ripple);
hold on
ylimit = gca().YLim;
plot(bin_centers(sig_times_ripple),ylimit(2)*ones(sum(sig_times_ripple),1),'.k')


ax3 = nexttile(3);
ax2Chil = ax2.Children;
copyobj(ax2Chil,ax3);
xticks(0:2:10)
xtickangle(0)
% ylabel('Events/s')
yticks([0 0.4 0.8])

set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
saveas(gcf,fullfile(fig_path,'sde_ripple_rates_on_v_off'),'jpg')
saveas(gcf,fullfile(fig_path,'sde_ripple_rates_on_v_off'),'pdf')


sde_n_off = [min(sum(~isnan(combined_hists_rates.sde(inds1,1:num_bins_to_plot)))) max(sum(~isnan(combined_hists_rates.sde(inds1,1:num_bins_to_plot))))]
sde_n_on = [min(sum(~isnan(combined_hists_rates.sde(inds2,1:num_bins_to_plot)))) max(sum(~isnan(combined_hists_rates.sde(inds2,1:num_bins_to_plot))))]

ripple_n_off = [min(sum(~isnan(combined_hists_rates.ripple(inds1,1:num_bins_to_plot)))) max(sum(~isnan(combined_hists_rates.ripple(inds1,1:num_bins_to_plot))))]
ripple_n_on = [min(sum(~isnan(combined_hists_rates.ripple(inds2,1:num_bins_to_plot)))) max(sum(~isnan(combined_hists_rates.ripple(inds2,1:num_bins_to_plot))))]