%% Generates Figure S8F, comparing the rates of forward and reverse replays over time in epochs with MEC active or inactive.
% NOTE: you must restrict analysis to rats 1,2,3,6, as only these animals have both laser ON and OFF epochs, and Jaws expression.

events_to_plot = [{'reverse_congruent_replays_spike'},{'forward_congruent_replays_spike'},{'sde'},{'ripples'}];
   % events_to_plot = [{'reverse_congruent_replays_ripple'},{'forward_congruent_replays_ripple'},{'sde'},{'ripples'}];

fig1 = figure();
fig1.Position = [600 600 300 250];
tiledlayout(2,2);

tiledlayout(2,2);
% Plot forward rates and reverse rates on top of each other
nexttile(2)
event_types = [{'reverse_congruent_replays_spike'},{'forward_congruent_replays_spike'}];

mean_off_rev = mean_combined_hists_rates_off.(event_types{1});
mean_on_rev = mean_combined_hists_rates_on.(event_types{1});
sem_off_rev = sem_combined_hists_rates_off.(event_types{1});
sem_on_rev = sem_combined_hists_rates_on.(event_types{1});

mean_off_for = mean_combined_hists_rates_off.(event_types{2});
mean_on_for = mean_combined_hists_rates_on.(event_types{2});
sem_off_for = sem_combined_hists_rates_off.(event_types{2});
sem_on_for = sem_combined_hists_rates_on.(event_types{2});

colors = [.4660 0.6740 0.1880; 0.4940 0.1840 0.5560];
transparency_pcnt = 0.99;
colors_transparent = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
    [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];

% Plot forward, off
H = shadedErrorBar(hist_bin_centers_plot,mean_off_for,sem_off_for,'lineprops','g'); hold on
H.patch.FaceColor = colors_transparent(1,:);
H.mainLine.Color = colors(1,:);
H.mainLine.LineWidth = 1;
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';

% Plot forward, on
H = shadedErrorBar(hist_bin_centers_plot,mean_on_for,sem_on_for,'lineprops','--g'); hold on
H.patch.FaceColor = colors_transparent(1,:);
H.mainLine.Color = colors(1,:);
H.mainLine.LineWidth = 1;
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';

% Plot reverse, off
H = shadedErrorBar(hist_bin_centers_plot,mean_off_rev,sem_off_rev,'lineprops','m'); hold on
H.patch.FaceColor = colors_transparent(2,:);
H.mainLine.Color = colors(2,:);
H.mainLine.LineWidth = 1;
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';

% Plot reverse, on
H = shadedErrorBar(hist_bin_centers_plot,mean_on_rev,sem_on_rev,'lineprops','--m'); hold on
H.patch.FaceColor = colors_transparent(2,:);
H.mainLine.Color = colors(2,:);
H.mainLine.LineWidth = 1;
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';
ylim(fig_ylim)
xlim([0 10])
xticks(0:2:10);
hold on
xticklabels({})
xtickangle(0);
ylabel('Events/s')

% Plot the difference in rate (forward-reverse), laser off and on
nexttile(4)
colors = [0 0 0; 0 0 0];
transparency_pcnt = 0.99;
colors_transparent = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
    [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];

event_types = {'congruent_forward_minus_reverse_spike'};
mean_off = mean_combined_hists_rates_off.(event_types{1});
mean_on = mean_combined_hists_rates_on.(event_types{1});
sem_off = sem_combined_hists_rates_off.(event_types{1});
sem_on = sem_combined_hists_rates_on.(event_types{1});
% Plot difference in rate, off
H = shadedErrorBar(hist_bin_centers_plot,mean_off,sem_off,'lineprops','k'); hold on
H.patch.FaceColor = colors_transparent(1,:);
H.mainLine.Color = colors(1,:);
H.mainLine.LineWidth = 1;
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';

% Plot difference in rate, on
H = shadedErrorBar(hist_bin_centers_plot,mean_on,sem_on,'lineprops','--k'); hold on
H.patch.FaceColor = colors_transparent(2,:);
H.mainLine.Color = colors(2,:);
H.mainLine.LineWidth = 1;
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';
yline(0);
ylim(fig_ylim_diff)
xticks(0:2:10);
xlim([0 10])
xtickangle(0);
xlabel('Time since arrival (s)')
ylabel('Events/s')

% Plot SDE rate, laser off and on
nexttile(1)
colors = [0 0 0; 0 0 0];
transparency_pcnt = 0.99;
colors_transparent = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
    [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];

event_types = {'sde'};
mean_off = mean_combined_hists_rates_off.(event_types{1});
mean_on = mean_combined_hists_rates_on.(event_types{1});
sem_off = sem_combined_hists_rates_off.(event_types{1});
sem_on = sem_combined_hists_rates_on.(event_types{1});

% Plot sde, off
H = shadedErrorBar(hist_bin_centers_plot,mean_off,sem_off,'lineprops','k'); hold on
H.patch.FaceColor = colors_transparent(1,:);
H.mainLine.Color = colors(1,:);
H.mainLine.LineWidth = 1;
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';

% Plot sde, on
H = shadedErrorBar(hist_bin_centers_plot,mean_on,sem_on,'lineprops','--k'); hold on
H.patch.FaceColor = colors_transparent(2,:);
H.mainLine.Color = colors(2,:);
H.mainLine.LineWidth = 1;
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';
xticks(0:2:10)
xticklabels({})
ylim([0 0.5])
xlim([0 10])
ylabel('Events/s')

off = combined_hists_rates.(event_types{1})(combined_table.laser_state==0,:);
on = combined_hists_rates.(event_types{1})(combined_table.laser_state==1,:); 
pval=nan(size(off,2),1);
for bin = 1:size(off,2);
    pval(bin) = ranksum(off(:,bin),on(:,bin));
end
plot(hist_bin_centers_plot(pval<0.05),ones(sum(pval<0.05),1),'.k')

% Plot ripple rate, laser off and on
nexttile(3)
colors = [0 0 0; 0 0 0];
transparency_pcnt = 0.99;
colors_transparent = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
    [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];

event_types = {'ripples'};
mean_off = mean_combined_hists_rates_off.(event_types{1});
mean_on = mean_combined_hists_rates_on.(event_types{1});
sem_off = sem_combined_hists_rates_off.(event_types{1});
sem_on = sem_combined_hists_rates_on.(event_types{1});

% Plot ripples, off
H = shadedErrorBar(hist_bin_centers_plot,mean_off,sem_off,'lineprops','k'); hold on
H.patch.FaceColor = colors_transparent(1,:);
H.mainLine.Color = colors(1,:);
H.mainLine.LineWidth = 1;
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';

% Plot ripples, on
H = shadedErrorBar(hist_bin_centers_plot,mean_on,sem_on,'lineprops','--k'); hold on
H.patch.FaceColor = colors_transparent(2,:);
H.mainLine.Color = colors(2,:);
H.mainLine.LineWidth = 1;
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';
xticks(0:2:10);
xlabel('Time since arrival (s)')
xtickangle(0);
ylim([0 0.5])
xlim([0 10])
ylabel('Events/s')

off = combined_hists_rates.(event_types{1})(combined_table.laser_state==0,:);
on = combined_hists_rates.(event_types{1})(combined_table.laser_state==1,:); 
pval=nan(size(off,2),1);
for bin = 1:size(off,2)
    pval(bin) = ranksum(off(:,bin),on(:,bin));
end

ylimit = gca().YLim;
plot(hist_bin_centers_plot(pval<0.05),ones(sum(pval<0.05),1),'.k')

set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode', 'auto');
fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures';
saveas(gcf,fullfile(fig_path,'Fig_S2'),'pdf')
saveas(gcf,fullfile(fig_path,'Fig_S2'),'jpg')


