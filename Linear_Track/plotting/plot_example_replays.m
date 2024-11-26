
% Probably this should be done on one session.
% Select which replays you want to look at:
candidate_events_to_plot = 'spike_filtered';
lower_limit_time_since_reward_zone_entry = 0;
upper_limit_time_since_reward_zone_entry = inf;
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
t_replay.score = t_replay.weighted_r.*t_replay.replay_ripple_power.*t_replay.range_normalized.*t_replay.fraction_of_excitatory_cells_participating;


%%
shiftSizeDecoding = 0.005;
windowSizeDecoding = 0.02;
spikeSampRate = 30000;
sessionNum=1;

plot_text = 0;
plot_ripple_power_zscored = 0;
plot_ripple_power = 0;
plot_bandpass_filtered_ripple = 0;
plot_spikeDensity = 0;

num_examples = 8;
% sorting_column1 = 'weighted_r';
sorting_column1 = 'replay_ripple_power';
sorting_column2 = 'weighted_r';
sorting_column3 = 'replay_ripple_power';
% sorting_column3 = 'fraction_of_excitatory_cells_participating';
% sorting_column4 = 'coverage_wu_normalized';


inds_f0_congruent_map1 = find(t_replay.best_map == 1 & t_replay.laser_state == 0 & t_replay.direction == 1 & t_replay.congruent_with_rat_location==1);
[~,index] = sortrows(t_replay(inds_f0_congruent_map1,:),{sorting_column1,sorting_column2},'descend');
inds_f0_congruent_map1 = inds_f0_congruent_map1(index);

inds_r0_congruent_map1 = find(t_replay.best_map == 1 & t_replay.laser_state == 0 & t_replay.direction == 2 & t_replay.congruent_with_rat_location==1);
[~,index] = sortrows(t_replay(inds_r0_congruent_map1,:),{sorting_column1,sorting_column2},'descend');
inds_r0_congruent_map1 = inds_r0_congruent_map1(index);

inds_f0_incongruent_map1 = find(t_replay.best_map == 1 & t_replay.laser_state == 0 & t_replay.direction == 1 & t_replay.incongruent_with_rat_location==1);
[~,index] = sortrows(t_replay(inds_f0_incongruent_map1,:),{sorting_column1,sorting_column2},'descend');
inds_f0_incongruent_map1 = inds_f0_incongruent_map1(index);

inds_r0_incongruent_map1 = find(t_replay.best_map == 1 & t_replay.laser_state == 0 & t_replay.direction == 2 & t_replay.incongruent_with_rat_location==1);
[~,index] = sortrows(t_replay(inds_r0_incongruent_map1,:),{sorting_column1,sorting_column2},'descend');
inds_r0_incongruent_map1 = inds_r0_incongruent_map1(index);


inds_f0_congruent_map2 = find(t_replay.best_map == 2 & t_replay.laser_state == 0 & t_replay.direction == 1 & t_replay.congruent_with_rat_location==1);
[~,index] = sortrows(t_replay(inds_f0_congruent_map2,:),{sorting_column1,sorting_column2},'descend');
inds_f0_congruent_map2 = inds_f0_congruent_map2(index);

inds_r0_congruent_map2 = find(t_replay.best_map == 2 & t_replay.laser_state == 0 & t_replay.direction == 2 & t_replay.congruent_with_rat_location==1);
[~,index] = sortrows(t_replay(inds_r0_congruent_map2,:),{sorting_column1,sorting_column2},'descend');
inds_r0_congruent_map2 = inds_r0_congruent_map2(index);

inds_f0_incongruent_map2 = find(t_replay.best_map == 2 & t_replay.laser_state == 0 & t_replay.direction == 1 & t_replay.incongruent_with_rat_location==1);
[~,index] = sortrows(t_replay(inds_f0_incongruent_map2,:),{sorting_column1,sorting_column2},'descend');
inds_f0_incongruent_map2 = inds_f0_incongruent_map2(index);

inds_r0_incongruent_map2 = find(t_replay.best_map == 2 & t_replay.laser_state == 0 & t_replay.direction == 2 & t_replay.incongruent_with_rat_location==1);
[~,index] = sortrows(t_replay(inds_r0_incongruent_map2,:),{sorting_column1,sorting_column2},'descend');
inds_r0_incongruent_map2 = inds_r0_incongruent_map2(index);


inds_f1_congruent_map1 = find(t_replay.best_map == 1 & t_replay.laser_state == 1 & t_replay.direction == 1 & t_replay.congruent_with_rat_location==1);
[~,index] = sortrows(t_replay(inds_f1_congruent_map1,:),{sorting_column1,sorting_column2},'descend');
inds_f1_congruent_map1 = inds_f1_congruent_map1(index);

inds_r1_congruent_map1 = find(t_replay.best_map == 1 & t_replay.laser_state == 1 & t_replay.direction == 2 & t_replay.congruent_with_rat_location==1);
[~,index] = sortrows(t_replay(inds_r1_congruent_map1,:),{sorting_column1,sorting_column2},'descend');
inds_r1_congruent_map1 = inds_r1_congruent_map1(index);

inds_f1_incongruent_map1 = find(t_replay.best_map == 1 & t_replay.laser_state == 1 & t_replay.direction == 1 & t_replay.incongruent_with_rat_location==1);
[~,index] = sortrows(t_replay(inds_f1_incongruent_map1,:),{sorting_column1,sorting_column2},'descend');
inds_f1_incongruent_map1 = inds_f1_incongruent_map1(index);

inds_r1_incongruent_map1 = find(t_replay.best_map == 1 & t_replay.laser_state == 1 & t_replay.direction == 2 & t_replay.incongruent_with_rat_location==1);
[~,index] = sortrows(t_replay(inds_r1_incongruent_map1,:),{sorting_column1,sorting_column2},'descend');
inds_r1_incongruent_map1 = inds_r1_incongruent_map1(index);


inds_f1_congruent_map2 = find(t_replay.best_map == 2 & t_replay.laser_state == 1 & t_replay.direction == 1 & t_replay.congruent_with_rat_location==1);
[~,index] = sortrows(t_replay(inds_f1_congruent_map2,:),{sorting_column1,sorting_column2},'descend');
inds_f1_congruent_map2 = inds_f1_congruent_map2(index);

inds_r1_congruent_map2 = find(t_replay.best_map == 2 & t_replay.laser_state == 1 & t_replay.direction == 2 & t_replay.congruent_with_rat_location==1);
[~,index] = sortrows(t_replay(inds_r1_congruent_map2,:),{sorting_column1,sorting_column2},'descend');
inds_r1_congruent_map2 = inds_r1_congruent_map2(index);

inds_f1_incongruent_map2 = find(t_replay.best_map == 2 & t_replay.laser_state == 1 & t_replay.direction == 1 & t_replay.incongruent_with_rat_location==1);
[~,index] = sortrows(t_replay(inds_f1_incongruent_map2,:),{sorting_column1,sorting_column2},'descend');
inds_f1_incongruent_map2 = inds_f1_incongruent_map2(index);

inds_r1_incongruent_map2 = find(t_replay.best_map == 2 & t_replay.laser_state == 1 & t_replay.direction == 2 & t_replay.incongruent_with_rat_location==1);
[~,index] = sortrows(t_replay(inds_r1_incongruent_map2,:),{sorting_column1,sorting_column2},'descend');
inds_r1_incongruent_map2 = inds_r1_incongruent_map2(index);

% if isempty(inds_r0_incongruent_map2)
%     return
% end

inds_f0_congruent_map1_plot = inds_f0_congruent_map1(1:min(num_examples,length(inds_f0_congruent_map1)));
inds_r0_congruent_map1_plot = inds_r0_congruent_map1(1:min(num_examples,length(inds_r0_congruent_map1)));
inds_f0_congruent_map2_plot = inds_f0_congruent_map2(1:min(num_examples,length(inds_f0_congruent_map2)));
inds_r0_congruent_map2_plot = inds_r0_congruent_map2(1:min(num_examples,length(inds_r0_congruent_map2)));
inds_f1_congruent_map1_plot = inds_f1_congruent_map1(1:min(num_examples,length(inds_f1_congruent_map1)));
inds_r1_congruent_map1_plot = inds_r1_congruent_map1(1:min(num_examples,length(inds_r1_congruent_map1)));
inds_f1_congruent_map2_plot = inds_f1_congruent_map2(1:min(num_examples,length(inds_f1_congruent_map2)));
inds_r1_congruent_map2_plot = inds_r1_congruent_map2(1:min(num_examples,length(inds_r1_congruent_map2)));




inds_f0_incongruent_map1_plot = inds_f0_incongruent_map1(1:min(num_examples,length(inds_f0_incongruent_map1)));
inds_r0_incongruent_map1_plot = inds_r0_incongruent_map1(1:min(num_examples,length(inds_r0_incongruent_map1)));
inds_f0_incongruent_map2_plot = inds_f0_incongruent_map2(1:min(num_examples,length(inds_f0_incongruent_map2)));
inds_r0_incongruent_map2_plot = inds_r0_incongruent_map2(1:min(num_examples,length(inds_r0_incongruent_map2)));
inds_f1_incongruent_map1_plot = inds_f1_incongruent_map1(1:min(num_examples,length(inds_f1_incongruent_map1)));
inds_r1_incongruent_map1_plot = inds_r1_incongruent_map1(1:min(num_examples,length(inds_r1_incongruent_map1)));
inds_f1_incongruent_map2_plot = inds_f1_incongruent_map2(1:min(num_examples,length(inds_f1_incongruent_map2)));
inds_r1_incongruent_map2_plot = inds_r1_incongruent_map2(1:min(num_examples,length(inds_r1_incongruent_map2)));

colors = [[115 82 68]./260; [194 150 130]./260];

load clusters
load Experiment_Information.mat
load Analysis_Information.mat
load zscored_sd_pyr
load zscored_ripple_power
load spikeDensity_pyr

% redo the zscoring of spike density
% Smooth spike density
smoothing_sigma = 0.001; % desired standard deviation of the Gaussian used for smoothing
stepSize = mean(diff(spikeDensity(:,1)))/spikeSampRate;
w = setUp_gaussFilt_sigma(smoothing_sigma,stepSize);
spikeDensity_smoothed = conv(spikeDensity(:,2),w,'same');
spikeDensity_smoothed = [spikeDensity(:,1),spikeDensity_smoothed];



for laser_state = 0:1 %0:1
    for direction = 1:2 %1:2
        for map = 1:2 %1:2
            for congruent = 1:2 %1:2
                figure('Position',[2084 744 800 400])
                tiledlayout(3,num_examples,"TileSpacing",'tight','Padding','tight')
                if direction==1 && map==1 && laser_state==0 && congruent==1
                    ind_replay = inds_f0_congruent_map1_plot;
                elseif direction==1 && map==1 && laser_state==0 && congruent==2
                    ind_replay = inds_f0_incongruent_map1_plot;
                elseif direction==2 && map==1 && laser_state==0 && congruent==1
                    ind_replay = inds_r0_congruent_map1_plot;
                elseif direction==2 && map==1 && laser_state==0 && congruent==2
                    ind_replay = inds_r0_incongruent_map1_plot;
                elseif direction==1 && map==2 && laser_state==0 && congruent == 1
                    ind_replay = inds_f0_congruent_map2_plot;
                elseif direction==1 && map==2 && laser_state==0 && congruent == 2
                    ind_replay = inds_f0_incongruent_map2_plot;
                elseif direction==2 && map==2 && laser_state==0 && congruent == 1
                    ind_replay = inds_r0_congruent_map2_plot;
                elseif direction==2 && map==2 && laser_state==0 && congruent == 2
                    ind_replay = inds_r0_incongruent_map2_plot;
                elseif direction==1 && map==1 && laser_state==1 && congruent == 1
                    ind_replay = inds_f1_congruent_map1_plot;
                elseif direction==1 && map==1 && laser_state==1 && congruent == 2
                    ind_replay = inds_f1_incongruent_map1_plot;
                elseif direction==2 && map==1 && laser_state==1 && congruent == 1
                    ind_replay = inds_r1_congruent_map1_plot;
                elseif direction==2 && map==1 && laser_state==1 && congruent == 2
                    ind_replay = inds_r1_incongruent_map1_plot;
                elseif direction==1 && map==2 && laser_state==1 && congruent == 1
                    ind_replay = inds_f1_congruent_map2_plot;
                elseif direction==1 && map==2 && laser_state==1 && congruent == 2
                    ind_replay = inds_f1_incongruent_map2_plot;
                elseif direction==2 && map==2 && laser_state==1 && congruent == 1
                    ind_replay = inds_r1_congruent_map2_plot;
                elseif direction==2 && map==2 && laser_state==1 && congruent == 2
                    ind_replay = inds_r1_incongruent_map2_plot;
                end
                num_examples_to_plot = min(num_examples,length(ind_replay));
                for i = 1:num_examples_to_plot
                    session_path = t_replay.session_str{ind_replay(i)};
                    %                 cd(session_path);
                    %                 load clusters
                    %                 load Experiment_Information.mat
                    %                 load Analysis_Information.mat
                    %                 load zscored_sd_pyr
                    %                 load zscored_ripple_power
                    %                 load spikeDensity_pyr

                    min_plot_range = 0.5*numSpatialBins - 20;
                    max_plot_range = 0.5*numSpatialBins + 20;
                    scaling_factor_bandpass = (max_plot_range-min_plot_range)/(500);
                    scaling_factor_power = (max_plot_range-min_plot_range)/(500);
                    % scaling_factor_raw = (max_plot_range-min_plot_range)/(max(ripple_replay(:,2))-min(ripple_replay(:,2)));
                    scaling_factor_power_zscored = (max_plot_range-min_plot_range)/(15);

                    %% Load and sort place fields for plotting raster plots
                    numSpatialBins = [1 numSpatialBins];
                    [ind_cluster,rateMap_smoothed,resort_index] = load_rateMapsForDecoding_1D_cm(clusters,1,Experiment_Information.Run_Times,meanFiringRateThr,spatialInfoThr,noiseOverlapThr,isolationThr,peakSnrThr,[numSpatialBins],1);

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
                    spikeDensity_replay_raw = compute_dataTemporalConcatenation(spikeDensity,timePoints_sub);
                    spikeDensity_replay = compute_dataTemporalConcatenation(zscored_sd,timePoints_sub);
                    maxSpikeDensity_replay = max(spikeDensity_replay(:,3));

                    spikeDensity_replay_2 = compute_dataTemporalConcatenation(spikeDensity_smoothed,timePoints_sub);




                    %limit ripple power to the replay:
                    ripple_replay = compute_dataTemporalConcatenation(zscored_ripple_power,timePoints_sub);
                    maxRipplePower_replay = max(ripple_replay (:,7));
                    posterior_left = t_replay.posterior_full_left{ind_replay(i,:)};
                    posterior_right = t_replay.posterior_full_right{ind_replay(i,:)};

                    numTimeBins = size(posterior_left,1);
                    numSpatialBins = size(posterior_left,2);

                    ripple_pseudo_time = linspace(1,numTimeBins,size(ripple_replay,1));
                    ripple_bandpass_pseudo_plot = ripple_replay(:,3)*scaling_factor_bandpass + numSpatialBins/2;
                    ripple_power_pseudo_plot = ripple_replay(:,5)*scaling_factor_power + numSpatialBins/2;
                    %     ripple_raw_pseudo_plot = ripple_replay(:,2)*scaling_factor_raw + numSpatialBins/2;
                    ripple_power_zscored_pseudo_plot = ripple_replay(:,7)*scaling_factor_power_zscored + numSpatialBins/2;

                    spikeDensity_pseudo_time = linspace(1,numTimeBins,size(spikeDensity_replay,1));
                    spikeDensity_pseudo_plot = spikeDensity_replay(:,3)*scaling_factor_power_zscored + numSpatialBins/2;
                    % Pull out replay metrics to add to plot:

                    ratPos_sub = t_replay.ratPos(ind_replay(i),1);
                    best_map_sub = t_replay.best_map(ind_replay(i));
                    time_since_drink_onset_sub = t_replay.time_since_drink_onset(ind_replay(i));


                    time_string = seconds(time_since_drink_onset_sub);
                    time_string.Format = 'mm:ss';
                    time_string = string(time_string);
                    candidate_event_num_string = num2str(ind_replay(i));


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
                    HSV(:,:,2) = HSV(:,:,2)*15;
                    posterior_color = hsv2rgb(HSV);

                    %subplot(2,num_examples,i)
                    nexttile(i)
                    imagesc(posterior_color)
                    set(gca,'YDir','normal');

                    hold on
                    if plot_bandpass_filtered_ripple == 1
                        plot(ripple_pseudo_time,ripple_bandpass_pseudo_plot,'color',[0.5 0.5 0.5])
                    end
                    %     if plot_raw_ripple ==1
                    %         plot(ripple_pseudo_time,ripple_raw_pseudo_plot,'color',[0.5 0.5 0.5])
                    %     end
                    if plot_ripple_power_zscored==1
                        plot(ripple_pseudo_time,ripple_power_zscored_pseudo_plot,'color',[194 150 130]./260,'LineWidth',2)
                    end
                    if plot_ripple_power==1
                        plot(ripple_pseudo_time,ripple_power_pseudo_plot,'color',[194 150 130]./260,'LineWidth',2)
                    end
                    if plot_spikeDensity==1
                        plot(spikeDensity_pseudo_time,spikeDensity_pseudo_plot,'color',[115 82 68]./260,'LineWidth',2)
                    end
                    hold on;
                    textColor = [0 0 0]
                    line([0 size(posterior_right,1)+0.5],[ratPos_sub(1),ratPos_sub(1)],'Color','k','LineWidth',2)

                    xticks([0.5 size(posterior_right,1)+0.5])
                    xticklabels({num2str(0) num2str(windowSizeDecoding + (size(posterior_right,1)-1)*shiftSizeDecoding,2)})
                    yticks([0.5 size(posterior_right,2) + 0.5])
                    yticklabels({num2str(0) num2str(size(posterior_right,2)*2)})
                    xtickangle(0);
                    set(gca,'FontSize',8)
                    text(size(posterior_right,1)/2,numSpatialBins-4,time_string,'HorizontalAlignment', 'center', 'VerticalAlignment', 'top','Color',textColor,'fontsize',10)
                    if plot_text==1
                        properties_to_label = {'weighted_r','replay_ripple_power','fraction_of_excitatory_cells_participating','mvl_bias_corrected',...
                            'coverage_wu_normalized','end_distance_from_rat','mean_HPD95_posterior_cropped_normalized'};
                        labels = cell(length(properties_to_label),1);
                        for property = 1:length(properties_to_label)
                            labels{property} = num2str(t_replay.(properties_to_label{property})(ind_replay(i)),2);
                        end
                        for text_line = 1:length(labels)
                            pos = 8 + 12*(text_line-1);
                            text(1,numSpatialBins-pos,labels{text_line},'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Color',textColor,'fontsize',8)
                        end

                    end
                    drawnow
                    % subplot(2,num_examples,i+num_examples)
                    nexttile(i+num_examples);

                    replay_spikes = [];
                    if best_map_sub == 2

                        for cell_id = 1:length(sorted_right_map_cells)
                            excitatory_cell = clusters(sorted_right_map_cells(cell_id)).Excitatory;
                            cell_spikes_phases = [clusters(sorted_right_map_cells(cell_id)).spkTime clusters(sorted_right_map_cells(cell_id)).spkRipplePhase];
                            cell_spikes_phases = compute_dataTemporalConcatenation(cell_spikes_phases,timePoints_sub);
                            cell_spikes = cell_spikes_phases(:,1);

                            replay_spikes = [replay_spikes; [cell_spikes_phases repmat(cell_id,[size(cell_spikes_phases,1),1]) repmat(excitatory_cell,[size(cell_spikes_phases,1),1])]];



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
                            xtickangle(0)
                            set(gca,'FontSize',8)
                        end

                    else
                        for cell_id = 1:length(sorted_left_map_cells)
                            excitatory_cell = clusters(sorted_left_map_cells(cell_id)).Excitatory;
                            cell_spikes_phases = [clusters(sorted_left_map_cells(cell_id)).spkTime clusters(sorted_left_map_cells(cell_id)).spkRipplePhase];
                            cell_spikes_phases = compute_dataTemporalConcatenation(cell_spikes_phases,timePoints_sub);
                            cell_spikes = cell_spikes_phases(:,1);

                            replay_spikes = [replay_spikes; [cell_spikes_phases repmat(cell_id,[size(cell_spikes_phases,1),1]) repmat(excitatory_cell,[size(cell_spikes_phases,1),1])]];




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
                            xtickangle(0)
                            set(gca,'FontSize',8)
                        end
                    end
                    drawnow


                    nexttile(i+2*num_examples)
                    plot(ripple_replay(:,1),ripple_replay(:,3)); hold on;
                    replay_spikes_ripple_phase = compute_dataInterpolation([ripple_replay(:,1),ripple_replay(:,3)],replay_spikes(:,1),[]);
                    hold on;
                    my_colors = parula(max(replay_spikes(:,3)));
                    replay_spikes_ripple_phase_colors = my_colors(replay_spikes(:,3),:);
                    replay_spikes_ripple_phase_colors(replay_spikes(:,4)==0,:) = repmat([0 0 0],[sum(replay_spikes(:,4)==0),1]);

                    scatter(replay_spikes_ripple_phase(:,1),replay_spikes_ripple_phase(:,2),20,replay_spikes_ripple_phase_colors,'filled')

                end
                set(gcf, 'Color', 'white','Renderer','painters');
                set(gcf, 'PaperPositionMode', 'auto');
                export_fig(['best_replays_direction_' num2str(direction) '_map' num2str(map), '_congruent' num2str(congruent) 'laser',num2str(laser_state)],'-jpeg')
                saveas(gcf,['best_replays_direction_' num2str(direction) '_map' num2str(map) '_congruent' num2str(congruent),'laser',num2str(laser_state)],'pdf')
                %             saveas(gcf,['best_replays_direction_' num2str(direction) '_map' num2str(map),'laser',num2str(laser_state)],'png')
            end
        end
    end
end