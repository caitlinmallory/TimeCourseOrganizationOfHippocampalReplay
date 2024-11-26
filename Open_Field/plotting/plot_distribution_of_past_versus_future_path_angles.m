add_pseudo_replays_for_rate_analysis
%%
downsample_stopping_periods = 0;
downsample_replays = 0;

replay_with_pseudo_replays_for_trial(isnan(replay_with_pseudo_replays_for_trial.drink_period_number),:) = [];
[C,ia, ic] = unique(replay_with_pseudo_replays_for_trial(:,{'unique_session_id','drink_period_number','home_event','rat_label','angle_between_past_future_trajectory','laser_state_binary'}));
height(C)

group1 = 'home';
group2 = 'away';
stopping_period_inds1 = find(C.home_event==1);
stopping_period_inds2 = find(C.home_event==0);
replay_inds1 = find(replay.home_event==1);
replay_inds2 = find(replay.home_event==0);

%%
figure('Position',[1137 776 150 125])

unique_angles = unique(C.angle_between_past_future_trajectory);

angle_bin_edges = linspace(0,180,11);
angle_bin_centers = (angle_bin_edges(1:end-1)' + angle_bin_edges(2:end)')./2;
h=hist(unique_angles,angle_bin_centers);
b = bar(angle_bin_centers,h);
b.BarWidth = 1;
b.FaceColor = [0.5 0.5 0.5];
b.EdgeColor = 'none';
hold on;
xlim([0 180])
xticks([0 60 120 180])
xtickangle(0)
xlabel('PF')
ylabel('Num. trials')
unique_angles1 = unique(C.angle_between_past_future_trajectory(stopping_period_inds1));
unique_angles2= unique(C.angle_between_past_future_trajectory(stopping_period_inds2));
yyaxis right
hold on
[f1,x1]=ecdf(unique_angles1);
plot(x1,f1,'b','LineWidth',1)
if ~isempty(unique_angles2)
[f2,x2]=ecdf(unique_angles2);
plot(x2,f2,'k','LineWidth',1)
hold off
box off
end
set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures';
saveas(gcf,fullfile(fig_path,['distribution_past_future_angles_with_overlay' group1 '_' group2 '_' num2str(rats)]),'jpeg')
saveas(gcf,fullfile(fig_path,['distribution_past_future_angles_with_overlay' group1 '_' group2 '_' num2str(rats)]),'pdf')

%%
figure('Position',[1137 776 200 150])

unique_angles = unique(C.angle_between_past_future_trajectory);

angle_bin_edges = linspace(0,180,11);
angle_bin_centers = (angle_bin_edges(1:end-1)' + angle_bin_edges(2:end)')./2;
h=hist(unique_angles,angle_bin_centers);
b = bar(angle_bin_centers,h);
b.BarWidth = 1;
b.FaceColor = [0.5 0.5 0.5];
b.EdgeColor = 'none';
hold on;
xlim([0 180])
xticks([0 60 120 180])
xtickangle(0)
xlabel('PF')
ylabel('Num. trials')
unique_angles1 = unique(C.angle_between_past_future_trajectory(stopping_period_inds1));
unique_angles2= unique(C.angle_between_past_future_trajectory(stopping_period_inds2));
yyaxis right
hold on
[f1,x1]=ecdf(unique_angles1);
plot(x1,f1,'b','LineWidth',2)
if ~isempty(unique_angles2)
[f2,x2]=ecdf(unique_angles2);
plot(x2,f2,'k','LineWidth',2)
hold off
end
box off
set(gcf,'PaperPositionMode','auto')
fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures';
saveas(gcf,fullfile(fig_path,['distribution_past_future_angles_without_overlay' group1 '_' group2 '_' num2str(rats)]),'jpeg')

disp('replay differences group1 v group2')
a = replay.angle_between_past_future_trajectory(replay_inds1);
b = replay.angle_between_past_future_trajectory(replay_inds2);
ecdf(a); hold on;
ecdf(b);
ranksum(a,b)

disp('stopping period differences group1 v group2')
a = C.angle_between_past_future_trajectory(stopping_period_inds1);
b = C.angle_between_past_future_trajectory(stopping_period_inds2);
figure()
ecdf(a); hold on;
ecdf(b);
[p, ~ , z] = ranksum(a,b);


%% Downsample to match angle bins (not performed)
% 
% if downsample_to_match_angle_bins == 1
%     if downsample_stopping_periods == 1
%         % Method 1: Downsample Stopping Periods
%         [downsampleindices1, downsampleindices2] = downsample_match_angles(C.angle_between_past_future_trajectory(stopping_period_inds1),C.angle_between_past_future_trajectory(stopping_period_inds2));
%         new_stopping_period_inds1 = stopping_period_inds1(downsampleindices1);
%         new_stopping_period_inds2 = stopping_period_inds2(downsampleindices2);
% 
%         new_stopping_period_inds = [new_stopping_period_inds1; new_stopping_period_inds2];
%         C.unique_session_id([new_stopping_period_inds1;new_stopping_period_inds2]);
%         C.rat_label([new_stopping_period_inds1;new_stopping_period_inds2]);
%         C.drink_period_number([new_stopping_period_inds1; new_stopping_period_inds2]);
% 
%         replay.meets_ds_criterion=zeros(height(replay),1);
%         for i = 1:length([new_stopping_period_inds1; new_stopping_period_inds2])
%             replay.meets_ds_criterion(replay.unique_session_id == C.unique_session_id(new_stopping_period_inds(i)) & ...
%                 replay.rat_label == C.rat_label(new_stopping_period_inds(i)) & ...
%                 replay.drink_period_number == C.drink_period_number([new_stopping_period_inds(i)]))=1;
%         end
% 
%         figure()
%         subplot(2,2,1)
%         a = C.angle_between_past_future_trajectory(stopping_period_inds1);
%         b = C.angle_between_past_future_trajectory(stopping_period_inds2);
%         ecdf(a); hold on;
%         ecdf(b);
%         ranksum(a,b)
%         [h,p,z] = kstest2(a,b);
%         subplot(2,2,2)
%         a = C.angle_between_past_future_trajectory(new_stopping_period_inds1);
%         b = C.angle_between_past_future_trajectory(new_stopping_period_inds2);
%         ecdf(a); hold on;
%         ecdf(b);
%         ranksum(a,b)
%         [h,p,z] = kstest2(a,b);
%         subplot(2,2,3)
%         a = replay.angle_between_past_future_trajectory(replay_inds1);
%         b = replay.angle_between_past_future_trajectory(replay_inds2);
%         ecdf(a); hold on;
%         ecdf(b);
%         ranksum(a,b)
%         [h,p,z] = kstest2(a,b);
%         subplot(2,2,4)
%         a = replay.angle_between_past_future_trajectory(intersect(find(replay.meets_ds_criterion==1), replay_inds1));
%         b = replay.angle_between_past_future_trajectory(intersect(find(replay.meets_ds_criterion==1), replay_inds2));
%         ecdf(a); hold on;
%         ecdf(b);
%         ranksum(a,b)
%         [h,p,z] = kstest2(a,b);
% 
%         replay = replay(replay.meets_ds_criterion==1,:);
% 
%     elseif downsample_replays==1
%         % Method 2: Downsample Replays
%         [downsampleindices1, downsampleindices2] = downsample_match_angles(replay.angle_between_past_future_trajectory(replay_inds1),replay.angle_between_past_future_trajectory(replay_inds2));
%         new_replay_inds1 = replay_inds1(downsampleindices1);
%         new_replay_inds2 = replay_inds2(downsampleindices2);
%         new_replay_inds = [new_replay_inds1; new_replay_inds2];
% 
%         figure()
%         subplot(1,2,1)
%         a = replay.angle_between_past_future_trajectory(replay_inds1);
%         b = replay.angle_between_past_future_trajectory(replay_inds2);
%         ecdf(a); hold on;
%         ecdf(b);
%         ranksum(a,b);
%         subplot(1,2,2)
%         a = replay.angle_between_past_future_trajectory(new_replay_inds1);
%         b = replay.angle_between_past_future_trajectory(new_replay_inds2);
%         ecdf(a); hold on;
%         ecdf(b);
%         ranksum(a,b);
% 
%         replay = replay(new_replay_inds,:);
%     end
% end