fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures';

% green_colormap = customcolormap(linspace(0,1,100), {'#014419','#1c7735','#5aae60','#a6db9d','#d7f1d6'});
% purple_colormap = customcolormap(linspace(0,1,100), {'#410149','#762a84','#9b6fac','#c1a5cd','#e7d4e8'});
% colorlimits = [0 1];

normalize_plots = 1;
color_limits_normalized = [0 12];
color_limits = [0 200]; % main
% color_limits = [0 40]; % sup
plot_mean_angular_hist = 0;
time_bin_1 =[0 3];

plot_time_bin_1_angular_hist = 1;
plot_time_bin_2_angular_hist = 1;
% windowSize = 0.5;
% windowShift = 0.5;
% start_time = -10;
% end_time = 10;
bin_start = (start_time:windowShift:(end_time-windowSize))';
bin_end = bin_start + windowSize;

bin_edges = [bin_start bin_end];
time_bin_centers = mean(bin_edges,2);

path_path_angle_bins = linspace(0, 180, 19);
replay_path_angle_bins = linspace(0, 180, 19);
angle_bin_centers = mean([replay_path_angle_bins(1:end-1)' replay_path_angle_bins(2:end)'],2);
replay_past_path_angle_map = nan(length(angle_bin_centers),length(angle_bin_centers));
replay_future_path_angle_map = nan(length(angle_bin_centers),length(angle_bin_centers));

for i = 1:length(angle_bin_centers)
    replay_past_angles = replay.meanAngDisplacement_pastPath(replay.time_since_real_drink_onset>time_bin_1(1) & replay.time_since_real_drink_onset<time_bin_1(2) ... 
        & replay.angle_between_past_future_trajectory >=path_path_angle_bins(i) & replay.angle_between_past_future_trajectory < path_path_angle_bins(i+1));
    replay_future_angles = replay.meanAngDisplacement_futPath(replay.time_since_real_drink_onset>time_bin_1(1) & replay.time_since_real_drink_onset<time_bin_1(2) ...
        & replay.angle_between_past_future_trajectory >=path_path_angle_bins(i) & replay.angle_between_past_future_trajectory < path_path_angle_bins(i+1));

    for j = 1:length(angle_bin_centers)
        replay_past_path_angle_map(i,j) = sum(replay_past_angles>=replay_path_angle_bins(j) & replay_past_angles<replay_path_angle_bins(j+1));
        replay_future_path_angle_map(i,j) = sum(replay_future_angles>=replay_path_angle_bins(j) & replay_future_angles<replay_path_angle_bins(j+1));  
    end
end

total_number_events_past = sum(replay_past_path_angle_map,2);
total_number_angles_past = sum(replay_past_path_angle_map,1);
pcnt_angles_past = total_number_angles_past./sum(total_number_events_past);

total_number_events_future = sum(replay_future_path_angle_map,2);
total_number_angles_future = sum(replay_future_path_angle_map,1);
pcnt_angles_future = total_number_angles_future./sum(total_number_events_future);


%%

figure('Position', [552 561 150 300])
tiledlayout(2,1,'TileSpacing','tight')
nexttile()
if normalize_plots==1
    imagesc(100.*((replay_past_path_angle_map./sum(replay_past_path_angle_map,2)))');
    caxis(color_limits_normalized)
else
imagesc(replay_past_path_angle_map')
caxis(color_limits)
end
colormap(magma)
set(gca,'YDir','normal');
yticks(0.5:6:(length(angle_bin_centers)+1))
yticklabels(num2str(replay_path_angle_bins(1:6:end)'));
xticks(0.5:6:(length(angle_bin_centers)+1))
xticklabels(num2str(path_path_angle_bins(1:6:end)'));
xtickangle(0)
xlabel('Path:path displacement')
ylabel('Replay:past displacement')
% title('Replay verus future path')

if normalize_plots==1
    imagesc(100.*((replay_past_path_angle_map./sum(replay_past_path_angle_map,2)))');
    caxis(color_limits_normalized)
else
imagesc(replay_past_path_angle_map')
caxis(color_limits)
end
colormap(magma)
set(gca,'YDir','normal');
axis square

yticks(0.5:6:(length(angle_bin_centers)+1))
yticklabels(num2str(replay_path_angle_bins(1:6:end)'));
xticks(0.5:6:(length(angle_bin_centers)+1))
xticklabels(num2str(path_path_angle_bins(1:6:end)'));
xtickangle(0)
xlabel('Path:path displacement')
ylabel('Replay:past displacement')
% title('Replay verus future path')

nexttile()
if normalize_plots==1
    imagesc(100.*((replay_future_path_angle_map./sum(replay_future_path_angle_map,2)))');
    caxis(color_limits_normalized)
else
imagesc(replay_future_path_angle_map')
caxis(color_limits)
end
colormap(magma)
set(gca,'YDir','normal');
yticks(0.5:6:(length(angle_bin_centers)+1))
yticklabels(num2str(replay_path_angle_bins(1:6:end)'));
xticks(0.5:6:(length(angle_bin_centers)+1))
xticklabels(num2str(path_path_angle_bins(1:6:end)'));
xtickangle(0)
xlabel('Path:path displacement')
ylabel('Replay:future displacement')
% title('Replay verus future path')
axis square

set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');

saveas(gcf,fullfile(fig_path,'path_path_angle_versus_replay_path_angle'),'pdf')

