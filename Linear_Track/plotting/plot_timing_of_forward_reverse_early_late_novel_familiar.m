% LINEAR TRACK
colors = [0.4660    0.6740    0.1880; 0.4940    0.1840    0.5560];

time_limits = [0 10];
trial_num_limits = [1 10; 30 40];
f_early = t_replay.time_since_reward_zone_entry(t_replay.laser_state==0 & t_replay.forward_congruent_with_rat_location==1 & t_replay.pass_number >= trial_num_limits(1,1) & t_replay.pass_number <= trial_num_limits(1,2) & t_replay.time_since_reward_zone_entry > time_limits(1) & t_replay.time_since_reward_zone_entry < time_limits(2));
r_early = t_replay.time_since_reward_zone_entry(t_replay.laser_state==0 & t_replay.reverse_congruent_with_rat_location==1 & t_replay.pass_number >= trial_num_limits(1,1) & t_replay.pass_number <= trial_num_limits(1,2) & t_replay.time_since_reward_zone_entry > time_limits(1) & t_replay.time_since_reward_zone_entry < time_limits(2));
f_late = t_replay.time_since_reward_zone_entry(t_replay.laser_state==0 & t_replay.forward_congruent_with_rat_location==1 & t_replay.pass_number >= trial_num_limits(2,1) & t_replay.pass_number <= trial_num_limits(2,2) & t_replay.time_since_reward_zone_entry > time_limits(1) & t_replay.time_since_reward_zone_entry < time_limits(2));
r_late = t_replay.time_since_reward_zone_entry(t_replay.laser_state==0 & t_replay.reverse_congruent_with_rat_location==1 & t_replay.pass_number >= trial_num_limits(2,1) & t_replay.pass_number <= trial_num_limits(2,2) & t_replay.time_since_reward_zone_entry > time_limits(1) & t_replay.time_since_reward_zone_entry < time_limits(2));

data = [f_early; r_early; f_late; r_late];
labels = [ones(size(f_early)); 2*ones(size(r_early)); 3*ones(size(f_late)); 4*ones(size(r_late))];

tbl = table();
tbl.time_since_reward_zone_entry = data;
tbl.direction = [ones(size(f_early)); 2*ones(size(r_early)); ones(size(f_late)); 2*ones(size(r_late))];
tbl.early_late = [ones(size(f_early)); ones(size(r_early)); 2*ones(size(f_late)); 2*ones(size(r_late))];

figure('Position', [1335 707 100 100])
b = boxchart(tbl.early_late,tbl.time_since_reward_zone_entry,'GroupByColor',tbl.direction,'MarkerSize',2);
colororder(colors)

yticks([0:2:10])
ylabel('Time since arrival (s)')
xticks([1 2])
xticklabels({})
xlim([0 3])
[p_early,~,z_early]  = ranksum(f_early,r_early)
[p_late,~,z_late] = ranksum(f_late,r_late)
set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
saveas(gcf,fullfile(final_fig_path,'timing_of_forward_reverse_replays_early_or_late_trials'),'jpeg')
saveas(gcf,fullfile(final_fig_path,'timing_of_forward_reverse_replays_early_or_late_trials'),'pdf')

f_novel= t_replay.time_since_reward_zone_entry(t_replay.laser_state==0 & t_replay.forward_congruent_with_rat_location==1 & cellfun(@(x) ismember(10.1,x), t_replay.session_flags)  &  t_replay.time_since_reward_zone_entry > time_limits(1) & t_replay.time_since_reward_zone_entry < time_limits(2));
r_novel = t_replay.time_since_reward_zone_entry(t_replay.laser_state==0 & t_replay.reverse_congruent_with_rat_location==1 & cellfun(@(x) ismember(10.1,x), t_replay.session_flags) &  t_replay.time_since_reward_zone_entry > time_limits(1) & t_replay.time_since_reward_zone_entry < time_limits(2));

f_familiar = t_replay.time_since_reward_zone_entry(t_replay.laser_state==0 & t_replay.forward_congruent_with_rat_location==1 & ~cellfun(@(x) ismember(10.1,x), t_replay.session_flags) & t_replay.time_since_reward_zone_entry > time_limits(1) & t_replay.time_since_reward_zone_entry < time_limits(2));
r_familiar = t_replay.time_since_reward_zone_entry(t_replay.laser_state==0 & t_replay.reverse_congruent_with_rat_location==1 & ~cellfun(@(x) ismember(10.1,x), t_replay.session_flags) & t_replay.time_since_reward_zone_entry > time_limits(1) & t_replay.time_since_reward_zone_entry < time_limits(2));

data = [f_novel; r_novel; f_familiar; r_familiar];
labels = [ones(size(f_novel)); 2*ones(size(r_novel)); 3*ones(size(f_familiar)); 4*ones(size(r_familiar))];

tbl = table();
tbl.time_since_reward_zone_entry = data;
tbl.direction = [ones(size(f_novel)); 2*ones(size(r_novel)); ones(size(f_familiar)); 2*ones(size(r_familiar))];
tbl.novel_familiar = [ones(size(f_novel)); ones(size(r_novel)); 2*ones(size(f_familiar)); 2*ones(size(r_familiar))];

figure('Position', [1335 707 100 100])
b = boxchart(tbl.novel_familiar,tbl.time_since_reward_zone_entry,'GroupByColor',tbl.direction,'MarkerSize',2);
colororder(colors)
xticklabels({})
xticks([])
yticks([0:2:10])
ylabel('Time since arrival (s)')
xticks([1 2])
% xticklabels({'1','2-5'})
xlim([0 3])

[p_novel, h, z_novel]  = ranksum(f_novel,r_novel)
[p_familiar, h, z_familiar] = ranksum(f_familiar,r_familiar)

set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
saveas(gcf,fullfile(final_fig_path,'timing_of_forward_reverse_replays_novel_or_familiar_trials'),'jpeg')
saveas(gcf,fullfile(final_fig_path,'timing_of_forward_reverse_replays_novel_or_familiar_trials'),'pdf')










%Forward congruent, reverse incongruent; Reverse congruent, forward
%incongruent
% figure()
% boxplot(data,labels)
% xticklabels({'Forward','Reverse','Forward','Reverse'})
% ylabel('Time since arrival (s)')
% box off
% saveas(gcf,'timing_of_forward_reverse_local_nonlocal_replays','jpeg')

