%% Version 1:

colors = [.4660 0.6740 0.1880; 0.4940 0.1840 0.5560;];
transparency_pcnt = 1;
colors_2 = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
    [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];

figure('Position',[1921 560 700 250])
tiledlayout(1,3,'TileSpacing','tight')


ax1 = nexttile(1);
hist(replay.meanAngDisplacement_HD_and_past_path);
xlim([0 180])
xticks([0 60 120 180])
xlabel('Angle between HD and past path')
ylabel('Num replays')

ax2=nexttile(2);
h = shadedErrorBar(bin_centers,grouped_binned_properties(1).binned_properties.meanAngDisplacement_HD.mean,...
    grouped_binned_properties(1).binned_properties.meanAngDisplacement_HD.sem,'lineprops','k'); hold on;
h.patch.FaceColor = [0 0 0];
h.mainLine.Color = [0 0 0];
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none'; 
hold on
h = shadedErrorBar(bin_centers,grouped_binned_properties(1).binned_properties.meanAngDisplacement_pastPath.mean,...
    grouped_binned_properties(1).binned_properties.meanAngDisplacement_pastPath.sem,'lineprops','k'); hold on;
h.patch.FaceColor = colors_2(2,:);
h.mainLine.Color = colors(2,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none'; 
ylim(ylim_deg)
yticks(ylim_deg_ticks);
ylabel('|Displacement|')
xlabel('Time (s)')
xlim([0 10])
xticks(0:2:10)
xtickangle(0)
legend({'Replay:HD';'Replay:past path'})


ax3=nexttile(3);

h = shadedErrorBar(bin_centers,grouped_binned_properties(2).binned_properties.meanAngDisplacement_HD.mean,...
    grouped_binned_properties(2).binned_properties.meanAngDisplacement_futPath.sem,'lineprops','--k'); hold on;
h.patch.FaceColor = [0 0 0];
h.mainLine.Color = [0 0 0];
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none'; 
xlim([0 10])
xticks(0:2:10)
xtickangle(0)
ylim(ylim_deg)
yticks(ylim_deg_ticks);


h = shadedErrorBar(bin_centers,grouped_binned_properties(2).binned_properties.meanAngDisplacement_pastPath.mean,...
    grouped_binned_properties(2).binned_properties.meanAngDisplacement_pastPath.sem,'lineprops','--k'); hold on;
h.patch.FaceColor = colors_2(2,:);
h.mainLine.Color = colors(2,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none'; 
xlabel('Time (s)')
xlim([0 10])
xticks(0:2:10)
xtickangle(0)
ylim(ylim_deg)
yticks(ylim_deg_ticks);

fig_title = 'replay_hd_angle_and_replay_past_angle';
set(gcf, 'Color', 'white','PaperPositionMode','auto');
saveas(gcf,fullfile(fig_path,[fig_title '_disp' num2str(replay_dispersionThr,2) '_rats_' num2str(rats)]),'jpg')
saveas(gcf,fullfile(fig_path,[fig_title '_disp' num2str(replay_dispersionThr,2) '_rats_' num2str(rats)]),'pdf')


[rho1,p1] = nancorr(replay.time_into_stopping_period(replay.time_into_stopping_period<=10),replay.meanAngDisplacement_futPath(replay.time_into_stopping_period<=10))
[rho1,p1] = nancorr(replay.time_into_stopping_period(replay.time_into_stopping_period<=10),replay.meanAngDisplacement_pastPath(replay.time_into_stopping_period<=10))
[rho1,p1] = nancorr(replay.time_into_stopping_period(replay.time_into_stopping_period<=10),replay.meanAngDisplacement_HD(replay.time_into_stopping_period<=10))


%% Compute the HD deviation for replays across a stopping period:
[C,ia, ic] = unique(replay(:,{'unique_session_id','drink_period_number','home_event','rat_label','drink_period_time','laser_state_binary'}));
unique_trials = height(C);
hd_deviations = nan(unique_trials,1);
for i = 1:unique_trials
    sub_replays = replay(replay.unique_session_id==C.unique_session_id(i)&replay.drink_period_number==C.drink_period_number(i),:);
    hd_deviations(i) = nanstd(sub_replays.ratHD);
end

a = rad2deg(hd_deviations);
mean(a)
median(a)
quantile(a,0.25)
quantile(a,0.75)