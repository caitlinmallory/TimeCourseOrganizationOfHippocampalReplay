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
ylim_rates = [0 1.5];


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
h.patch.FaceColor = colors(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.patch.FaceColor = colors_2(1,:);
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';
hold on
ylabel('Events/s')
ylim(ylim_rates)
yticks(ylim_rates(1):0.1:ylim_rates(2))
xticks(0:2:10)
xtickangle(0);

hold on;
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
h.patch.FaceColor = colors(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.patch.FaceColor = colors_2(1,:);
h.edge(1).Color = 'none';
h.edge(2).Color = 'none';
hold on


