%%
normalize_plots = 1;
figure('Position', [552 561 200 250])
tiledlayout(1,2,'TileSpacing','tight')
nexttile(1)
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


%%
normalize_plots = 1;
nexttile()
if normalize_plots==1
    imagesc(100.*((replay_past_path_angle_map./sum(replay_past_path_angle_map,2)))');
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
ylabel('Replay:past displacement')
% title('Replay verus future path')



set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');

saveas(gcf,fullfile(fig_path,'path_path_angle_versus_replay_path_angle'),'pdf')