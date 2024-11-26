function [corr_value_left,p_value_left,corr_value_right,p_value_right] = rank_order_correlation(clusters,event_times,x_centers)
peak_fr_thr = 1;

place_fields_left = nan(length(clusters),length(clusters(1).runs(1).directions(1).rateMap_smoothed));
place_fields_right = nan(length(clusters),length(clusters(1).runs(1).directions(2).rateMap_smoothed));
peak_fr = nan(length(clusters),1);

for i = 1:length(clusters)
place_fields_left(i,:) = [clusters(i).runs(1).directions(1).rateMap_smoothed];
place_fields_right(i,:) = [clusters(i).runs(1).directions(2).rateMap_smoothed];
peak_fr(i,:) = [clusters(i).runs(1).directional_peak_firing_rate];
end

clusters_to_keep =  find([clusters.Excitatory]==1 & peak_fr' >=peak_fr_thr);
clusters_to_discard = find([clusters.Excitatory]==0 | peak_fr'<peak_fr_thr);


% The sort_idx are the 'templates'.
[sort_idx_left, sort_idx_right] = plot_all_directional_fields(place_fields_left', place_fields_right', [], [], x_centers, 0);


% Reverse the direction of the template for leftward runs!
sort_idx_left = flipud(sort_idx_left);

sort_idx_left(ismember(sort_idx_left,clusters_to_discard)) = [];
sort_idx_right(ismember(sort_idx_right,clusters_to_discard)) = [];


[corr_value_left, p_value_left] = spearman_median(clusters,event_times,sort_idx_left,0);
[corr_value_right, p_value_right] = spearman_median(clusters,event_times,sort_idx_right,0);






















% 
% keyboard
% % sort_idx_left(ismember(sort_idx_left,clusters_to_discard)) = [];
% % sort_idx_right(ismember(sort_idx_right,clusters_to_discard)) = [];
% 
% % for i = 1:length(event_times)
% %     event_times = [candidateEvents.spike_events(i,1) candidateEvents.spike_events(i,2)];
% 
%     %determine which cells spiked during this candidate event
%     participating_cells = [];
%     replay_spike_times = [];
%     for cluster = 1:length(clusters)
%         cell_spikes = compute_dataTemporalConcatenation(clusters(cluster).spkTime,event_times);
%         if ~isempty(cell_spikes)
%             %cell_spikes = cell_spikes(1);
%             cell_spikes = median(cell_spikes);
%             participating_cells = [participating_cells; cluster];
%             replay_spike_times = [replay_spike_times; cell_spikes];
%         end
%     end
% 
%     replay_spike_times(ismember(participating_cells,clusters_to_discard)) = [];
%     participating_cells(ismember(participating_cells,clusters_to_discard)) = [];
%     
%     [a,sort_idx] = sort(replay_spike_times);
%     observed_order = participating_cells(sort_idx);
% 
%     template_left = sort_idx_left(ismember(sort_idx_left,participating_cells));
%     template_right = sort_idx_right(ismember(sort_idx_right,participating_cells));
%     
%     template_left = [template_left  [1:length(template_left)]'];
%     template_right = [template_right  [1:length(template_right)]'];
% 
% 
%     order_with_respect_to_left_template = nan(length(template_left),1);
%     order_with_respect_to_right_template = nan(length(template_right),1);
% 
% 
%     for cell_num = 1:length(observed_order)
%          order_with_respect_to_left_template(cell_num) = template_left(template_left(:,1)==observed_order(cell_num),2);
%          order_with_respect_to_right_template(cell_num) = template_right(template_right(:,1)==observed_order(cell_num),2);
%     end
% 
% 
%     [rho_left,p_left] = corr(template_left(:,2),order_with_respect_to_left_template);
%     [rho_right,p_right] = corr(template_right(:,2),order_with_respect_to_right_template); 
% 
%     num_shuffles = 5000;
%     shuffled_rho_left = nan(num_shuffles,1);
%     shuffled_rho_right = nan(num_shuffles,1);
%     for num_shuffle = 1:num_shuffles
%         shuffled_cell_order = randperm(size(template_left,1),size(template_left,1));
%         shuffled_rho_left = corr(template_left(:,2),shuffled_cell_order');
%         shuffled_rho_right = corr(template_right(:,2),shuffled_cell_order');
%     end
% 
%     spearmans = struct();
%     spearmans.rho_left = rho_left;
%     spearmans.rho_right = rho_right;
%     spearmans.p_left = p_left;
%     spearmans.p_right = p_right;
%     spearmans.left_rho_shuffles = shuffled_rho_left;
%     spearmans.right_rho_shuffles = shuffled_rho_right;
% 
%     spearmans.left_map_reverse = 0;
%     spearmans.left_map_forward = 0;
%     spearmans.right_map_reverse = 0;
%     spearmans.right_map_forward = 0;
%     replay = 0;
% 
%     if abs(rho_left) >= abs(rho_right)
%         map = 1;
%     elseif abs(rho_right) > abs(rho_left)
%         map = 2;
%     end
% 
%     if map == 1 && rho_left < quantile(shuffled_rho_left,0.025) 
%         spearmans.left_map_reverse = 1;
%         spearmans.replay = 1;
%     end
%     if map == 1 && rho_left > quantile(shuffled_rho_left,0.975)
%         spearmans.left_map_forward = 1;
%         spearmans.replay = 1;
%     end
%     if map == 2 && rho_right < quantile(shuffled_rho_right,0.025) 
%         spearmans.right_map_reverse = 1;
%         spearmans.replay = 1;
%     end
%     if map == 2 && rho_right > quantile(shuffled_rho_right,0.975)
%         spearmans.right_map_forward = 1;
%         spearmans.replay = 1;
%     end
% 
% 
% 
% 
% %     
% %     subplot(1,2,1)
% %     plot(template_left(:,2),order_with_respect_to_left_template,'ok')
% %     subplot(1,2,2)
% %     plot(template_right(:,2),order_with_respect_to_right_template,'ok')
% % 
% %     if spearmans.left_map_forward == 1
% %         title('left_forward');
% %     end
% %     if spearmans.right_map_forward == 1
% %         title('right_forward');
% %     end
% %     if spearmans.left_map_reverse == 1
% %         title('left_reverse');
% %     end
% %     if spearmans.right_map_reverse == 1
% %         title('right_reverse');
% %     end
% 
% 
