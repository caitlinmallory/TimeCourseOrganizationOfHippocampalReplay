
fig1 = figure();
fig1.Position = [600 600 300 250];
tiledlayout(2,2);

% Plot forward rates and reverse rates on top of each other
nexttile(2)
event_types = [{'reverse_congruent_replays_spike'},{'forward_congruent_replays_spike'}];

mean_rev_group1 = mean_combined_hists_rates_group1.(event_types{1});
mean_rev_group2 = mean_combined_hists_rates_group2.(event_types{1});
sem_rev_group1 = sem_combined_hists_rates_group1.(event_types{1});
sem_rev_group2 = sem_combined_hists_rates_group2.(event_types{1});

mean_for_group1 = mean_combined_hists_rates_group1.(event_types{2});
mean_for_group2 = mean_combined_hists_rates_group2.(event_types{2});
sem_for_group1 = sem_combined_hists_rates_group1.(event_types{2});
sem_for_group2 = sem_combined_hists_rates_group2.(event_types{2});


colors = [.4660 0.6740 0.1880; 0.4940 0.1840 0.5560];
transparency_pcnt = 0.99;
colors_transparent = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
    [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];


% Plot forward, group1
H = shadedErrorBar(hist_bin_centers_plot,mean_for_group1,sem_for_group1,'lineprops','g'); hold on
H.patch.FaceColor = colors_transparent(1,:);
H.mainLine.Color = colors(1,:);
H.mainLine.LineWidth = 1;
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';

% Plot forward, group2
H = shadedErrorBar(hist_bin_centers_plot,mean_for_group2,sem_for_group2,'lineprops','--g'); hold on
H.patch.FaceColor = colors_transparent(1,:);
H.mainLine.Color = colors(1,:);
H.mainLine.LineWidth = 1;
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';

% Plot reverse, group1
H = shadedErrorBar(hist_bin_centers_plot,mean_rev_group1,sem_rev_group1,'lineprops','m'); hold on
H.patch.FaceColor = colors_transparent(2,:);
H.mainLine.Color = colors(2,:);
H.mainLine.LineWidth = 1;
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';

% Plot reverse, group2
H = shadedErrorBar(hist_bin_centers_plot,mean_rev_group2,sem_rev_group2,'lineprops','--m'); hold on
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

% Plot the difference in rate (forward-reverse), group1 and group2
nexttile(4)
colors = [0 0 0; 0 0 0];
transparency_pcnt = 0.99;
colors_transparent = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
    [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];

event_types = {'congruent_forward_minus_reverse_spike'};
% event_types = {'congruent_forward_minus_reverse_ripple'};
mean_group1 = mean_combined_hists_rates_group1.(event_types{1});
mean_group2 = mean_combined_hists_rates_group2.(event_types{1});
sem_group1 = sem_combined_hists_rates_group1.(event_types{1});
sem_group2 = sem_combined_hists_rates_group2.(event_types{1});
% Plot difference in rate, off
H = shadedErrorBar(hist_bin_centers_plot,mean_group1,sem_group1,'lineprops','k'); hold on
H.patch.FaceColor = colors_transparent(1,:);
H.mainLine.Color = colors(1,:);
H.mainLine.LineWidth = 1;
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';

% Plot difference in rate, on
H = shadedErrorBar(hist_bin_centers_plot,mean_group2,sem_group2,'lineprops','--k'); hold on
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

group1 = combined_hists_rates.(event_types{1})(combined_table.laser_state==0&combined_table.session_half==1,:);
group2 = combined_hists_rates.(event_types{1})(combined_table.laser_state==0&combined_table.session_half==2,:); 
pval=nan(size(group1,2),1);
for bin = 1:size(group1,2);
    pval(bin) = ranksum(group1(:,bin),group2(:,bin));
end

ylimit = gca().YLim;
plot(hist_bin_centers_plot(pval<0.05),fig_ylim_diff(2).*ones(sum(pval<0.05),1),'.k')

% Plot SDE rate, group1 and on
nexttile(1)
colors = [0 0 0; 0 0 0];
transparency_pcnt = 0.99;
colors_transparent = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
    [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];

event_types = {'sde'};
mean_group1 = mean_combined_hists_rates_group1.(event_types{1});
mean_group2 = mean_combined_hists_rates_group2.(event_types{1});
sem_group1 = sem_combined_hists_rates_group1.(event_types{1});
sem_group2 = sem_combined_hists_rates_group2.(event_types{1});

% Plot sde, off
H = shadedErrorBar(hist_bin_centers_plot,mean_group1,sem_group1,'lineprops','k'); hold on
H.patch.FaceColor = colors_transparent(1,:);
H.mainLine.Color = colors(1,:);
H.mainLine.LineWidth = 1;
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';

% Plot sde, on
H = shadedErrorBar(hist_bin_centers_plot,mean_group2,sem_group2,'lineprops','--k'); hold on
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

group1 = combined_hists_rates.(event_types{1})(combined_table.laser_state==0&combined_table.session_half==1,:);
group2 = combined_hists_rates.(event_types{1})(combined_table.laser_state==0&combined_table.session_half==2,:); 
pval=nan(size(group1,2),1);
for bin = 1:size(group1,2);
    pval(bin) = ranksum(group1(:,bin),group2(:,bin));
end


ylimit = gca().YLim;
plot(hist_bin_centers_plot(pval<0.05),0.5*ones(sum(pval<0.05),1),'.k')


% Plot ripple rate, group1 and on
nexttile(3)
colors = [0 0 0; 0 0 0];
transparency_pcnt = 0.99;
colors_transparent = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
    [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];

event_types = {'ripples'};
mean_group1 = mean_combined_hists_rates_group1.(event_types{1});
mean_group2 = mean_combined_hists_rates_group2.(event_types{1});
sem_group1 = sem_combined_hists_rates_group1.(event_types{1});
sem_group2 = sem_combined_hists_rates_group2.(event_types{1});


% Plot ripples, off
H = shadedErrorBar(hist_bin_centers_plot,mean_group1,sem_group1,'lineprops','k'); hold on
H.patch.FaceColor = colors_transparent(1,:);
H.mainLine.Color = colors(1,:);
H.mainLine.LineWidth = 1;
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';

% Plot ripples, on
H = shadedErrorBar(hist_bin_centers_plot,mean_group2,sem_group2,'lineprops','--k'); hold on
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

group1 = combined_hists_rates.(event_types{1})(combined_table.laser_state==0 & combined_table.session_half==1,:);
group2 = combined_hists_rates.(event_types{1})(combined_table.laser_state==0 & combined_table.session_half==2,:); 
pval=nan(size(group1,2),1);
for bin = 1:size(group1,2);
    pval(bin) = ranksum(group1(:,bin),group2(:,bin));
end

ylimit = gca().YLim;
plot(hist_bin_centers_plot(pval<0.05),0.5.*ones(sum(pval<0.05),1),'.k')


set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode', 'auto');
fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures';
saveas(gcf,fullfile(fig_path,'Fig_Supp_FirstSecondHalf'),'pdf')
saveas(gcf,fullfile(fig_path,'Fig_Supp_FirstSecondHalf'),'jpg')


