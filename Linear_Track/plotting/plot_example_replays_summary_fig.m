
% Probably this should be done on one session.
% Select which replays you want to look at:
candidate_events_to_plot = 'spike_filtered';
lower_limit_time_since_reward_zone_entry = 0;
upper_limit_time_since_reward_zone_entry = 10;
use_duration_og = 1;
override_coverage_thr = 0.2;
override_weighted_r_thr = 0.6;
override_posterior_diff_thr = 0.33;

if strcmp(candidate_events_to_plot,'ripple_filtered')
    ripple_power_criterion = 3;
    sde_amplitude_criterion = -inf;
elseif strcmp(candidate_events_to_plot,'spike_filtered')
    ripple_power_criterion = -inf;
    sde_amplitude_criterion = 3;
else
    ripple_power_criterion = -inf;
    sde_amplitude_criterion = -inf;
end

t = load_replays_from_individual_session(candidate_events_to_plot,use_duration_og,override_coverage_thr,override_weighted_r_thr,override_posterior_diff_thr,sde_amplitude_criterion,ripple_power_criterion,lower_limit_time_since_reward_zone_entry,upper_limit_time_since_reward_zone_entry);
t.weighted_r = abs(t.weighted_r);
t_replay = t(t.replay==1,:);

%%
shiftSizeDecoding = 0.005;
windowSizeDecoding = 0.02;
spikeSampRate = 30000;
sessionNum=1;

plot_text = 1;
plot_ripple_power = 1;
plot_spikeDensity = 1;

num_examples = 4;
sorting_column = 'weighted_r';
inds_f0_map1 = find(t_replay.best_map == 1 & t_replay.laser_state == 0 & t_replay.direction == 1 & t_replay.incongruent_with_rat_location==1);
[~,index] = sortrows(t_replay(inds_f0_map1,:),sorting_column,'descend');
inds_f0_map1 = inds_f0_map1(index);

inds_r0_map1 = find(t_replay.best_map == 1 & t_replay.laser_state == 0 & t_replay.direction == 2 & t_replay.incongruent_with_rat_location==1);
[~,index] = sortrows(t_replay(inds_r0_map1,:),sorting_column,'descend');
inds_r0_map1 = inds_r0_map1(index);

inds_f0_map2 = find(t_replay.best_map == 2 & t_replay.laser_state == 0 & t_replay.direction == 1 & t_replay.incongruent_with_rat_location==1);
[~,index] = sortrows(t_replay(inds_f0_map2,:),sorting_column,'descend');
inds_f0_map2 = inds_f0_map2(index);

inds_r0_map2 = find(t_replay.best_map == 2 & t_replay.laser_state == 0 & t_replay.direction == 2 & t_replay.incongruent_with_rat_location==1);
[~,index] = sortrows(t_replay(inds_r0_map2,:),sorting_column,'descend');
inds_r0_map2 = inds_r0_map2(index);

inds_f1_map1 = find(t_replay.best_map == 1 & t_replay.laser_state == 1 & t_replay.direction == 1 & t_replay.incongruent_with_rat_location==1);
[~,index] = sortrows(t_replay(inds_f1_map1,:),sorting_column,'descend');
inds_f1_map1 = inds_f1_map1(index);

inds_r1_map1 = find(t_replay.best_map == 1 & t_replay.laser_state == 1 & t_replay.direction == 2 & t_replay.incongruent_with_rat_location==1);
[~,index] = sortrows(t_replay(inds_r1_map1,:),sorting_column,'descend');
inds_r1_map1 = inds_r1_map1(index);

inds_f1_map2 = find(t_replay.best_map == 2 & t_replay.laser_state == 1 & t_replay.direction == 1 & t_replay.incongruent_with_rat_location==1);
[~,index] = sortrows(t_replay(inds_f1_map2,:),sorting_column,'descend');
inds_f1_map2 = inds_f1_map2(index);

inds_r1_map2 = find(t_replay.best_map == 2 & t_replay.laser_state == 1 & t_replay.direction == 2 & t_replay.incongruent_with_rat_location==1);
[~,index] = sortrows(t_replay(inds_r1_map2,:),sorting_column,'descend');
inds_r1_map2 = inds_r1_map2(index);

% inds_f0_plot = inds_f0(randperm(numel(inds_f0),num_examples));
% inds_r0_plot = inds_r0(randperm(numel(inds_r0),num_examples));
inds_f0_map1_plot = inds_f0_map1(1:min(num_examples,length(inds_f0_map1)));
inds_r0_map1_plot = inds_r0_map1(1:min(num_examples,length(inds_r0_map1)));
inds_f0_map2_plot = inds_f0_map2(1:min(num_examples,length(inds_f0_map2)));
inds_r0_map2_plot = inds_r0_map2(1:min(num_examples,length(inds_r0_map2)));
inds_f1_map1_plot = inds_f1_map1(1:min(num_examples,length(inds_f1_map1)));
inds_r1_map1_plot = inds_r1_map1(1:min(num_examples,length(inds_r1_map1)));
inds_f1_map2_plot = inds_f1_map2(1:min(num_examples,length(inds_f1_map2)));
inds_r1_map2_plot = inds_r1_map2(1:min(num_examples,length(inds_r1_map2)));

figure('Position',[2094 569 1001 645])
tiledlayout(4,8)
ind_replay = [inds_f0_map1_plot; inds_f1_map1_plot; inds_r0_map1_plot; inds_r1_map1_plot];
for i = 1:length(ind_replay)

    session_path = t_replay.session_str{ind_replay(i)};
    cd(session_path);
    load clusters
    load Experiment_Information.mat
    load Analysis_Information.mat
    load zscored_sd_pyr
    load zscored_ripple_power

    min_plot_range = 0.5*numSpatialBins - 20;
    max_plot_range = 0.5*numSpatialBins + 20;
    % scaling_factor_bandpass = (max_plot_range-min_plot_range)/(max(ripple_replay(:,3))-min(ripple_replay(:,3)));
    % scaling_factor_raw = (max_plot_range-min_plot_range)/(max(ripple_replay(:,2))-min(ripple_replay(:,2)));
    scaling_factor_power = (max_plot_range-min_plot_range)/(15);

    %% Load and sort place fields for plotting raster plots
    [ind_cluster,rateMap_smoothed,resort_index] = load_rateMapsForDecoding_1D_cm(clusters,1,Experiment_Information.Run_Times,meanFiringRateThr,spatialInfoThr,noiseOverlapThr,isolationThr,peakSnrThr,numSpatialBins,1);

    left_map_cells = ind_cluster{1};
    right_map_cells = ind_cluster{2};

    left_maps = rateMap_smoothed(left_map_cells,1:numSpatialBins(2));
    [~, left_maps_peak] = max(left_maps,[],2);
    [left_map_sorted_peaks,left_maps_sort_idx] = sort(left_maps_peak);

    right_maps = rateMap_smoothed(right_map_cells,numSpatialBins(2)+1:end);
    [~, right_maps_peak] = max(right_maps,[],2);
    [right_map_sorted_peaks,right_maps_sort_idx] = sort(right_maps_peak);

    sorted_left_map_cells = left_map_cells(left_maps_sort_idx);
    sorted_right_map_cells = right_map_cells(right_maps_sort_idx);

    timePoints_sub = t_replay.timePoints_og(ind_replay(i),:);

    %limit spike density to the replay:
    spikeDensity_replay = compute_dataTemporalConcatenation(zscored_sd,timePoints_sub);
    maxSpikeDensity_replay = max(spikeDensity_replay(:,3));

    %limit spike density to the replay:
    ripple_replay = compute_dataTemporalConcatenation(zscored_ripple_power,timePoints_sub);
    maxRipplePower_replay = max(ripple_replay (:,7));

    posterior_left = t_replay.posterior_full_left{ind_replay(i,:)};
    posterior_right = t_replay.posterior_full_right{ind_replay(i,:)};

    numTimeBins = size(posterior_left,1);
    numSpatialBins = size(posterior_left,2);

    ripple_pseudo_time = linspace(1,numTimeBins,size(ripple_replay,1));
    %     ripple_bandpass_pseudo_plot = ripple_replay(:,3)*scaling_factor_bandpass + numSpatialBins/2;
    %     ripple_raw_pseudo_plot = ripple_replay(:,2)*scaling_factor_raw + numSpatialBins/2;
    ripple_power_pseudo_plot = ripple_replay(:,7)*scaling_factor_power + numSpatialBins/2;

    spikeDensity_pseudo_time = linspace(1,numTimeBins,size(spikeDensity_replay,1));
    spikeDensity_pseudo_plot = spikeDensity_replay(:,3)*scaling_factor_power + numSpatialBins/2;
    % Pull out replay metrics to add to plot:

    weighted_r_sub = t_replay.weighted_r(ind_replay(i));
    max_jump_distance_sub = t_replay.max_jump_distance(ind_replay(i));
    ave_jump_distance_sub = t_replay.mean_jump_distance(ind_replay(i));
    range_sub = t_replay.range(ind_replay(i));
    range_normalized_sub = t_replay.range_normalized(ind_replay(i));
    coverage_wu_sub = t_replay.coverage_wu(ind_replay(i));
    ave_sharpness_at_peak_sub = t_replay.mean_sharpness_at_peak(ind_replay(i));
    best_map_sub = t_replay.best_map(ind_replay(i));
    duration_sub = t_replay.duration(ind_replay(i));
    ratPos_sub = t_replay.ratPos(ind_replay(i),1);
    startsLocal_sub = t_replay.startsLocal(ind_replay(i));
    total_posterior_right_sub = t_replay.total_posterior_right(ind_replay(i));
    total_posterior_left_sub = t_replay.total_posterior_left(ind_replay(i));
    time_since_drink_onset_sub = t_replay.time_since_drink_onset(ind_replay(i));
    cum_coverage_sub = t_replay.cum_coverage(ind_replay(i));
    replay_ripple_power_sub = t_replay.replay_ripple_power(ind_replay(i));
    replay_spikeDensity_power_sub = t_replay.replay_spikeDensity_power(ind_replay(i));
    laser_state = t_replay.laser_state(ind_replay(i));

    time = {floor(((timePoints_sub(1)/spikeSampRate)-(Experiment_Information.Segments(sessionNum).Times(1)/spikeSampRate))/60),round(rem(((timePoints_sub(1)/spikeSampRate)-(Experiment_Information.Segments(sessionNum).Times(1)/spikeSampRate))/60,1)*60)};
    if time{2}>=10
        time_string = strcat(num2str(time{1}),':',num2str(time{2}));
    else
        time_string = strcat(num2str(time{1}),':0',num2str(time{2}));
    end
    candidate_event_num_string = num2str(ind_replay(i));
    info_string1 = ['PD ' num2str(time_since_drink_onset_sub,2)];
    info_string2 = ['Ripple ' num2str(replay_ripple_power_sub,2)];
    info_string3 = ['Spike ' num2str(replay_spikeDensity_power_sub,2)];
    info_string4 = ['R ' num2str(weighted_r_sub,2)];
    info_string5 = ['Cov ' num2str(range_sub,2)];
    %                 info_string2 = ['Spike ' num2str(replay_spikeDensity_power_sub,2)];

    if ratPos_sub(1) < 1
        ratPos_sub(1) = 1;
    end
    if ratPos_sub(1) > size(posterior_left,2)
        ratPos_sub(1) = size(posterior_left,2);
    end

    clear posterior_color
    posterior_color(:,:,1) = (1-posterior_left' + ones(size(posterior_left')))/2;
    posterior_color(:,:,2) = (1-posterior_left' + 1-posterior_right')/2;
    posterior_color(:,:,3) = (ones(size(posterior_left')) + 1-posterior_right')/2;

    %change saturation
    HSV = rgb2hsv(posterior_color);
    HSV(:,:,2) = HSV(:,:,2)*10;
    posterior_color = hsv2rgb(HSV);

    subplot(2,num_examples,i)
    imagesc(posterior_color)
    set(gca,'YDir','normal');
    textColor = 'k';
    hold on
    %     if plot_bandpass_filtered_ripple == 1
    %         plot(ripple_pseudo_time,ripple_bandpass_pseudo_plot,'color',[0.5 0.5 0.5])
    %     end
    %     if plot_raw_ripple ==1
    %         plot(ripple_pseudo_time,ripple_raw_pseudo_plot,'color',[0.5 0.5 0.5])
    %     end
    if plot_ripple_power==1
        plot(ripple_pseudo_time,ripple_power_pseudo_plot,'color','k')
    end
    if plot_spikeDensity==1
        plot(spikeDensity_pseudo_time,spikeDensity_pseudo_plot,'m')
    end
    hold on;

    line([0 size(posterior_right,1)+0.5],[ratPos_sub(1),ratPos_sub(1)],'Color',textColor,'LineWidth',2)

    xticks([0.5 size(posterior_right,1)+0.5])
    xticklabels({num2str(0) num2str(windowSizeDecoding + (size(posterior_right,1)-1)*shiftSizeDecoding,2)})
    yticks([0.5 size(posterior_right,2) + 0.5])
    yticklabels({num2str(0) num2str(size(posterior_right,2)*2)})
    set(gca,'FontSize',12)
    text(size(posterior_right,1)/2,numSpatialBins-4,time_string,'HorizontalAlignment', 'center', 'VerticalAlignment', 'top','Color',textColor,'fontsize',12)
    if plot_text==1
        %text(1,numSpatialBins-4,candidate_event_num_string,'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Color',textColor,'fontsize',10)
        text(1,numSpatialBins-16, info_string1,'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Color',textColor,'fontsize',10)
        text(1,numSpatialBins-28, info_string2,'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Color',textColor,'fontsize',10)
        text(1,numSpatialBins-40, info_string3,'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Color',textColor,'fontsize',10)
        text(1,numSpatialBins-52, info_string4, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Color',textColor,'fontsize',10)
        text(1,numSpatialBins-64, info_string5, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Color',textColor,'fontsize',10)
        %                         text(1,numSpatialBins(2)-76, info_string6, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Color',textColor,'fontsize',10)
    end
    drawnow
    subplot(2,num_examples,i+num_examples)

    if best_map_sub == 2
        for cell_id = 1:length(sorted_right_map_cells)

            cell_spikes = compute_dataTemporalConcatenation(clusters(sorted_right_map_cells(cell_id)).spkTime,timePoints_sub);
            cell_spikes = cell_spikes./30000-timePoints_sub(1)/30000;

            for spike = 1:length(cell_spikes)
                y_pos = [cell_id-1 cell_id];
                x_pos = [cell_spikes(spike) cell_spikes(spike)];
                plot(x_pos,y_pos,'-k','LineWidth',2); hold on

            end
            ylim([0 length(sorted_right_map_cells)])
            xlim([0 (timePoints_sub(2)/30000 - timePoints_sub(1)/30000)]);
            xticks([0 (timePoints_sub(2)/30000 - timePoints_sub(1)/30000)]);
            xticklabels({num2str(0) num2str(windowSizeDecoding + (size(posterior_right,1)-1)*shiftSizeDecoding,2)})
            set(gca,'FontSize',12)
        end

    else
        for cell_id = 1:length(sorted_left_map_cells)

            cell_spikes = compute_dataTemporalConcatenation(clusters(sorted_left_map_cells(cell_id)).spkTime,timePoints_sub);
            cell_spikes = cell_spikes./30000-timePoints_sub(1)/30000;

            for spike = 1:length(cell_spikes)
                y_pos = [cell_id-1 cell_id];
                x_pos = [cell_spikes(spike) cell_spikes(spike)];
                plot(x_pos,y_pos,'-k','LineWidth',2); hold on
            end

            ylim([0 length(sorted_left_map_cells)])
            xlim([0 (timePoints_sub(2)/30000 - timePoints_sub(1)/30000)]);
            xticks([0 (timePoints_sub(2)/30000 - timePoints_sub(1)/30000)]);
            xticklabels({num2str(0) num2str(windowSizeDecoding + (size(posterior_right,1)-1)*shiftSizeDecoding,2)})
            set(gca,'FontSize',12)
        end
    end
    drawnow
end

saveas(gcf,['best_replays_direction_' num2str(direction) '_map' num2str(map)],'png')



