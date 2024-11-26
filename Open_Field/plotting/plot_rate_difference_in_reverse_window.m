%% Look at the overall rate of retrospective or prospective replay on each trial type.
bins = 1:20;
rate_future_inds1 = nansum(combined_hists.future(inds1,bins),2)./(min(10,C.drink_period_time(inds1,:)));
rate_past_inds1 = nansum(combined_hists.past(inds1,bins),2)./(min(10,C.drink_period_time(inds1)));
rate_future_inds2 = nansum(combined_hists.future(inds2,bins),2)./(min(10,C.drink_period_time(inds2,:)));
rate_past_inds2 = nansum(combined_hists.past(inds2,bins),2)./(min(10,C.drink_period_time(inds2,:)));

[pval_future_inds1_inds2, h, zval_future_inds1_inds2] = ranksum(rate_future_inds1,rate_future_inds2)
[pval_past_inds1_inds2, h, zval_past_inds1_inds2] = ranksum(rate_past_inds1,rate_past_inds2)

bins = find(bin_centers> 3 & bin_centers<=10)
late_rate_future_inds1 = nansum(combined_hists.future(inds1,bins),2)./(min(10,C.drink_period_time(inds1,:)));
late_rate_past_inds1 = nansum(combined_hists.past(inds1,bins),2)./(min(10,C.drink_period_time(inds1)));
late_rate_future_inds2 = nansum(combined_hists.future(inds2,bins),2)./(min(10,C.drink_period_time(inds2,:)));
late_rate_past_inds2 = nansum(combined_hists.past(inds2,bins),2)./(min(10,C.drink_period_time(inds2,:)));


late_rate_difference_inds1 = late_rate_future_inds1-late_rate_past_inds1;
late_rate_difference_inds2 = late_rate_future_inds2-late_rate_past_inds2;
ranksum(late_rate_difference_inds1,late_rate_difference_inds2)

[p,h,z] = signrank(late_rate_difference_inds1)
[p,h,z] = signrank(late_rate_difference_inds2)
[p,h,z] = ranksum(late_rate_difference_inds1,late_rate_difference_inds2)


data = [nanmean(late_rate_difference_inds1),nanmean(late_rate_difference_inds2)];
err = [nanstd(late_rate_difference_inds1)/sqrt(nansum(~isnan(late_rate_difference_inds1))),...
    nanstd(late_rate_difference_inds2)/sqrt(nansum(~isnan(late_rate_difference_inds2)))];

figure('Position',[1986 1051 100 100])
b1 = bar(1,data(1)); hold on;
e1 = errorbar(1,data(1),err(1),'k','linestyle','none')
e1.CapSize = 4;
hold on
b2 = bar(2,data(2)); hold on;
e2 = errorbar(2,data(2),err(2),'k','linestyle','none')
e2.CapSize = 4;
b1.FaceColor = [0 0 0];
b1.EdgeColor = 'none';
b2.FaceColor = [0 0 0];
b2.EdgeColor = 'none';



box off
% xticks([1.5 4.5])
xticklabels({})
ylabel('Events/s')


set(gcf, 'Color', 'white','Renderer','painters', 'PaperPositionMode', 'auto');
fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures';
if home_trials_only==1
    title_ending = 'home_trials';
elseif away_trials_only==1
    title_ending = 'away_trials';
else
    title_ending = 'all_trials';
end
saveas(gcf,fullfile(fig_path,['late_window_rate_difference_' category1 '_' category2 'trials_' title_ending]),'pdf')

