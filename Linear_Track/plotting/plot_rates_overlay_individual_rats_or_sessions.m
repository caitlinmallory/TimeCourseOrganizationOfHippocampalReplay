
fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures';
unit = 'rats';
unit = 'sessions';

figure('Position',[500 500 125 100]);
plot_individual_n = 0;
plot_mean_and_sem = 1;

if strcmp(unit,'rats')
    unique_n = unique(combined_table.rat_label);
    combined_table.unique_n = combined_table.rat_label;
elseif strcmp(unit,'sessions')
    unique_n = unique(combined_table.unique_session_id);
    combined_table.unique_n = combined_table.unique_session_id;
end

data_n_all = nan(length(unique_n),length(hist_bin_centers_plot));
for j = 1:length(unique_n)
    combined_table_sub = combined_table(combined_table.unique_n == unique_n(j) & combined_table.laser_state==0,:);
    sum(combined_table_sub.reverse_congruent_replays_spike)
    sum(combined_table_sub.forward_congruent_replays_spike)
    combined_hists_rates_off_sub = combined_hists_rates(combined_table.unique_n == unique_n(j) & combined_table.laser_state==0,:);

    data = combined_hists_rates_off_sub.congruent_forward_minus_reverse_spike;
    data_mean = nanmean(data);
    data_n_all(j,:) = data_mean;

    data_n = sum(~isnan(data));
    data_sem = nanstd(data)./sqrt(data_n);
    if plot_individual_n==1
    plot(hist_bin_centers_plot,data_mean); hold on;
    end
end

pvals = nan(length(hist_bin_centers_plot),1)
for j = 1:length(hist_bin_centers_plot)
   [h,pvals(j)] = ttest(data_n_all(:,j))
end

ylim([-0.05 0.08])
yticks([-0.05 0.05])
yticklabels({})
yline(0)
xlim([0 10])
xticks([0:2:10])
xtickangle(0)
hold on
if plot_mean_and_sem == 1
data_n_all_mean = nanmean(data_n_all);
data_n_all_sem = nanstd(data_n_all)./(sqrt(size(data_n_all,1)));
shadedErrorBar(hist_bin_centers_plot,data_n_all_mean,data_n_all_sem);
end
ylimit=gca().YLim; ylimit = ylimit(2);
plot(hist_bin_centers_plot(pvals<0.05),ylimit*ones(sum(pvals<0.05),1),'.k')


set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
saveas(gcf,fullfile(fig_path,['rate_difference_by_' unit]),'pdf')