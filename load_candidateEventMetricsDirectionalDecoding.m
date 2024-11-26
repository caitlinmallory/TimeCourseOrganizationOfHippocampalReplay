function [eventMetrics] = load_candidateEventMetricsDirectionalDecoding(params,posterior_full,posterior,posterior_spread_full,posterior_spread_cropped,event_timepoints,event_inds,filtering_map,indNaN,Position_Data,Position_Data_scaled,Experiment_Information)
%filtering map is option. If left empty, compute which map has the greater
%posterior, and assign it as 'best_map'. If filtering_map is included, it
%will take the place of 'best_map'.

%indNaN is also optional. It's a list of bins that should be ignored (due
%to high posterior spread or high animal velocity). These were calculated
%in load_replayEvents_cm (John's bottom up filtering method). If left
%empty, we'll consider all time points during the replay when calculating
%things like weighted_r. If indNan is included, those bins will be removed.

% params.num_bins_to_mask_track_ends = 0; % number of bins to exclude on either side of the track when calculating total posterior in each map.
% params.num_bins_to_mask_local_posterior = 15; % number of bins to exlude around the rat's current location when calculated total posterior in each map.
% params.num_time_bins_to_check_for_local_start = 4;
% distance_thr_for_local_start;  % max number of bins between rat's current location and the start of the replay sequence to be considered local.
% params.spikeSampRate = 30000;
% params.posBinWidth = 2;
% params.decodingWindowShift = 0.005;
% params.davidson_distance_from_line_to_search = 6; % num of bins above OR below the line of best fit to sum posterior over.

% remove time bins from the posterior specified indNaN;

plot_fig = 0;

eventMetrics = struct;

posterior(indNaN,:) = [];
event_timepoints(indNaN,:) = [];

ratPos = [Position_Data(event_inds(1),2),Position_Data(event_inds(1),3)];
ratPos_scaled = [Position_Data_scaled(event_inds(1),2),Position_Data_scaled(event_inds(1),3)];

reward_locations = [0 Experiment_Information.maze_size];
[~,reward_zone] = min([abs(reward_locations(1) - ratPos(1)),abs(reward_locations(2) - ratPos(1))]);

num_spatial_bins = size(posterior_full,2)/2;
left_bins = 1:num_spatial_bins;
right_bins = num_spatial_bins+1:size(posterior_full,2);
eventMetrics.posterior_full_left = posterior_full(:,left_bins);
eventMetrics.posterior_full_right = posterior_full(:,right_bins);

if isempty(posterior)

    dispersion = nan;
    replay_score = [nan; nan];
    posterior_left = nan;
    posterior_right = nan;
    weighted_r = nan;
    replay_sequence_com = nan;
    replay_sequence_peak = nan;
    sharpness = nan;
    mean_sharpness_at_peak = nan;
    median_sharpness_at_peak = nan;
    max_jump_distance_peak = nan;
    mean_jump_distance_peak = nan;
    median_jump_distance_peak = nan;
    total_posterior_left = nan;
    total_posterior_right = nan;
    best_map = nan;
    coverage_wu = nan;
    coverage_wu_normalized = nan;
    coverage_wu_range = nan;
    coverage_wu_range_normalized = nan;
    posterior_range = nan;
    posterior_range_normalized = nan;
    best_map_distance_over_time = nan;
    cum_coverage = nan;
    cum_coverage_normalized = nan;
    %     replay_score = [nan nan nan nan];

    percent_posterior_in_left_map_full = nan;
    percent_posterior_in_right_map_full = nan;
    percent_non_local_posterior_in_left_map_full = nan;
    percent_non_local_posterior_in_right_map_full = nan;
    percent_posterior_ends_masked_in_left_map_full = nan;
    percent_posterior_ends_masked_in_right_map_full = nan;


    percent_posterior_in_left_map = nan;
    percent_posterior_in_right_map = nan;
    percent_non_local_posterior_in_left_map = nan;
    percent_non_local_posterior_in_right_map = nan;
    percent_posterior_ends_masked_in_left_map = nan;
    percent_posterior_ends_masked_in_right_map = nan;

    startsLocal = nan;
    start_distance_from_rat = nan;
    start_distance_from_rat_normalized = nan;
    endsLocal = nan;
    end_distance_from_rat = nan;
    end_distance_from_rat_normalized = nan;
    direction = nan;
    forward_congruent_with_rat_location = nan;
    reverse_congruent_with_rat_location = nan;
    congruent_with_rat_location = nan;
    forward_incongruent_with_rat_location = nan;
    reverse_incongruent_with_rat_location = nan;
    incongruent_with_rat_location = nan;
    posterior_spread_full = nan;
    posterior_spread_cropped = nan;

    mean_std_posterior_cropped = nan;
    median_std_posterior_cropped = nan;
    mean_std_posterior_cropped_normalized = nan;
    median_std_posterior_cropped_normalized = nan;
    mean_CI95_posterior_cropped = nan;
    median_CI95_posterior_cropped = nan;
    mean_CI95_posterior_cropped_normalized = nan;
    median_CI95_posterior_cropped_normalized = nan;
    mean_HPD95_posterior_cropped = nan;
    median_HPD95_posterior_cropped = nan;
    mean_HPD95_posterior_cropped_normalized = nan;
    median_HPD95_posterior_cropped_normalized = nan;


    mean_std_posterior_full = nan;
    median_std_posterior_full = nan;
    mean_std_posterior_full_normalized = nan;
    median_std_posterior_full_normalized = nan;
    mean_CI95_posterior_full = nan;
    median_CI95_posterior_full = nan;
    mean_CI95_posterior_full_normalized = nan;
    median_CI95_posterior_full_normalized = nan;
    mean_HPD95_posterior_full = nan;
    median_HPD95_posterior_full = nan;
    mean_HPD95_posterior_full_normalized = nan;
    median_HPD95_posterior_full_normalized = nan;

else

    % compute the duration of the replay event
    replay_duration = (max(event_timepoints) - min(event_timepoints))/params.spikeSampRate;

    posterior_left = posterior(:,left_bins);
    posterior_right = posterior(:,right_bins);
    posterior_full_left = posterior_full(:,left_bins);
    posterior_full_right = posterior_full(:,right_bins);

    % compute the amount of posterior in each map, using either the full
    % posterior, or best segment:

    % Best segment:
    total_posterior = sum(sum(posterior_left(:,(params.num_bins_to_mask_track_ends+1):(num_spatial_bins-params.num_bins_to_mask_track_ends))))...
        + sum(sum(posterior_right(:,(params.num_bins_to_mask_track_ends+1):(num_spatial_bins-params.num_bins_to_mask_track_ends))));
    % total_posterior will sum to 1*num_time_bins if you don't mask the
    % ends of the track
    total_posterior_left = sum(sum(posterior_left(:,(params.num_bins_to_mask_track_ends+1):(num_spatial_bins-params.num_bins_to_mask_track_ends))))/total_posterior;
    total_posterior_right = sum(sum(posterior_right(:,(params.num_bins_to_mask_track_ends+1):(num_spatial_bins-params.num_bins_to_mask_track_ends))))/total_posterior;
    if total_posterior_left >= total_posterior_right
        best_map = 1;
    else
        best_map = 2;
    end
    percent_posterior_in_left_map = total_posterior_left/(total_posterior_left+total_posterior_right);
    percent_posterior_in_right_map = total_posterior_right/(total_posterior_left+total_posterior_right);

    % Calculate the amount of posterior in the left or right maps, but mask
    % out the ends
    total_posterior = sum(sum(posterior_left(:,(params.num_bins_to_mask_track_ends+1):(num_spatial_bins-params.num_bins_to_mask_track_ends))))...
        + sum(sum(posterior_right(:,(params.num_bins_to_mask_track_ends+1):(num_spatial_bins-params.num_bins_to_mask_track_ends))));
    % total_posterior will sum to 1*num_time_bins if you don't mask the
    % ends of the track
    total_posterior_ends_masked_left = sum(sum(posterior_left(:,(params.num_bins_to_mask_track_ends+1):(num_spatial_bins-params.num_bins_to_mask_track_ends))))/total_posterior;
    total_posterior_ends_masked_right = sum(sum(posterior_right(:,(params.num_bins_to_mask_track_ends+1):(num_spatial_bins-params.num_bins_to_mask_track_ends))))/total_posterior;
    total_posterior_ends_masked = total_posterior_ends_masked_left + total_posterior_ends_masked_right;
    percent_posterior_ends_masked_in_left_map = total_posterior_ends_masked_left/total_posterior_ends_masked;
    percent_posterior_ends_masked_in_right_map = total_posterior_ends_masked_right/total_posterior_ends_masked;


    % Calculate the amount of posterior in the left or right maps, but mask
    % out the rat's current location.
    bins_to_mask_nonlocal_posterior = find(abs((1:num_spatial_bins)-ratPos_scaled(1)) < params.num_bins_to_mask_local_posterior);
    posterior_left_nonlocal = posterior_left;
    posterior_left_nonlocal(:,bins_to_mask_nonlocal_posterior) = [];
    posterior_right_nonlocal = posterior_right;
    posterior_right_nonlocal(:,bins_to_mask_nonlocal_posterior) = [];
    total_nonlocal_posterior_left = sum(sum(posterior_left_nonlocal));
    total_nonlocal_posterior_right = sum(sum(posterior_right_nonlocal));
    total_nonlocal_posterior = total_nonlocal_posterior_left + total_nonlocal_posterior_right;
    percent_non_local_posterior_in_left_map = total_nonlocal_posterior_left/total_nonlocal_posterior;
    percent_non_local_posterior_in_right_map = total_nonlocal_posterior_right/total_nonlocal_posterior;

    %%
       total_posterior_full = sum(sum(posterior_full_left(:,(params.num_bins_to_mask_track_ends+1):(num_spatial_bins-params.num_bins_to_mask_track_ends))))...
        + sum(sum(posterior_full_right(:,(params.num_bins_to_mask_track_ends+1):(num_spatial_bins-params.num_bins_to_mask_track_ends))));
    % total_posterior will sum to 1*num_time_bins if you don't mask the
    % ends of the track
    total_posterior_full_left = sum(sum(posterior_full_left(:,(params.num_bins_to_mask_track_ends+1):(num_spatial_bins-params.num_bins_to_mask_track_ends))))/total_posterior_full;
    total_posterior_full_right = sum(sum(posterior_full_right(:,(params.num_bins_to_mask_track_ends+1):(num_spatial_bins-params.num_bins_to_mask_track_ends))))/total_posterior_full;
    if total_posterior_full_left >= total_posterior_full_right
        best_map = 1;
    else
        best_map = 2;
    end
    percent_posterior_in_left_map_full = total_posterior_full_left/(total_posterior_full_left+total_posterior_full_right);
    percent_posterior_in_right_map_full = total_posterior_full_right/(total_posterior_full_left+total_posterior_full_right);

    % Calculate the amount of posterior in the left or right maps, but mask
    % out the ends
    total_posterior_full = sum(sum(posterior_full_left(:,(params.num_bins_to_mask_track_ends+1):(num_spatial_bins-params.num_bins_to_mask_track_ends))))...
        + sum(sum(posterior_full_right(:,(params.num_bins_to_mask_track_ends+1):(num_spatial_bins-params.num_bins_to_mask_track_ends))));
    % total_posterior will sum to 1*num_time_bins if you don't mask the
    % ends of the track
    total_posterior_ends_masked_left_full = sum(sum(posterior_full_left(:,(params.num_bins_to_mask_track_ends+1):(num_spatial_bins-params.num_bins_to_mask_track_ends))))/total_posterior_full;
    total_posterior_ends_masked_right_full = sum(sum(posterior_full_right(:,(params.num_bins_to_mask_track_ends+1):(num_spatial_bins-params.num_bins_to_mask_track_ends))))/total_posterior_full;
    total_posterior_ends_masked_full = total_posterior_ends_masked_left_full + total_posterior_ends_masked_right_full;
    percent_posterior_ends_masked_in_left_map_full = total_posterior_ends_masked_left_full/total_posterior_ends_masked_full;
    percent_posterior_ends_masked_in_right_map_full = total_posterior_ends_masked_right_full/total_posterior_ends_masked_full;


    % Calculate the amount of posterior in the left or right maps, but mask
    % out the rat's current location.
    bins_to_mask_nonlocal_posterior_full = find(abs((1:num_spatial_bins)-ratPos_scaled(1)) < params.num_bins_to_mask_local_posterior);
    posterior_left_nonlocal_full = posterior_full_left;
    posterior_left_nonlocal_full(:,bins_to_mask_nonlocal_posterior_full) = [];
    posterior_right_nonlocal_full = posterior_full_right;
    posterior_right_nonlocal_full(:,bins_to_mask_nonlocal_posterior_full) = [];
    total_nonlocal_posterior_full_left = sum(sum(posterior_left_nonlocal_full));
    total_nonlocal_posterior_full_right = sum(sum(posterior_right_nonlocal_full));
    total_nonlocal_posterior_full = total_nonlocal_posterior_full_left + total_nonlocal_posterior_full_right;
    percent_non_local_posterior_in_left_map_full = total_nonlocal_posterior_full_left/total_nonlocal_posterior_full;
    percent_non_local_posterior_in_right_map_full = total_nonlocal_posterior_full_right/total_nonlocal_posterior_full;%%



    % If filtering map was included, ignore the above and compute all metrics
    % on the desired map.
    if ~isempty(filtering_map)
        best_map = filtering_map;
    end

    if best_map == 1
        best_posterior = posterior_left;
        best_posterior_full = posterior_full_left;
    else
        best_posterior = posterior_right;
        best_posterior_full = posterior_full_right;
    end

    % For each time bin, find the location corresponding to the peak posterior,
    % and the center of mass of the posterior. Also find the posterior
    % associated with each of these locations.
    [posterior_at_peak, replay_sequence_peak] = max(best_posterior,[],2);
    replay_sequence_com = nan(size(best_posterior,1),1);
    posterior_at_com = nan(size(best_posterior,1),1);
    for t = 1:size(best_posterior,1)
        replay_sequence_com(t) = round(compute_centerOfMass_1D(best_posterior(t,:)));
        posterior_at_com(t) = best_posterior(t,replay_sequence_com(t));
    end

    %% Posterior spread needs WORK! NOT WORKING WELL.
    % Add the posterior spread:
    posterior_spread_full = mean(posterior_spread_full(:,best_map))/num_spatial_bins;
    posterior_spread_cropped = mean(posterior_spread_cropped(:,best_map))/num_spatial_bins;

    % sum the posterior in each bin, and find what's 80% of that value.
    % then see how many

    best_posterior_cropped_normalized = best_posterior./sum(best_posterior,2);
    best_posterior_full_normalized = best_posterior_full./sum(best_posterior_full,2);

    [mean_posterior_cropped,std_posterior_cropped,CI95_posterior_cropped,HPD95_posterior_cropped] = compute_posterior_statistics(best_posterior_cropped_normalized);
    [mean_posterior_full,std_posterior_full,CI95_posterior_full,HPD95_posterior_full] = compute_posterior_statistics(best_posterior_full_normalized);

    if plot_fig==1
        subplot(1,3,1)
        imagesc(best_posterior_cropped_normalized'); set(gca,'YDir','normal'); hold on
        plot([1:size(best_posterior_cropped_normalized,1)], mean_posterior_cropped,'.k','MarkerSize',30); hold on
        for i = 1:size(best_posterior,1)
            line([i, i], [mean_posterior_cropped(i) - std_posterior_cropped(i), mean_posterior_cropped(i) + std_posterior_cropped(i)],'color','k')
        end
        title('mean and STD')
        subplot(1,3,2)
        imagesc(best_posterior_cropped_normalized'); set(gca,'YDir','normal'); hold on
        plot([1:size(best_posterior_cropped_normalized,1)], mean_posterior_cropped,'.k','MarkerSize',30); hold on
        for i = 1:size(best_posterior,1)
            line([i, i], [CI95_posterior_cropped(i,1) CI95_posterior_cropped(i,2)],'color','k')
        end
        title('mean and 95% CI')
        subplot(1,3,3)
        imagesc(best_posterior_cropped_normalized'); set(gca,'YDir','normal'); hold on
        plot([1:size(best_posterior_cropped_normalized,1)], mean_posterior_cropped,'.k','MarkerSize',30); hold on
        for i = 1:size(best_posterior,1)
            line([i, i], [HPD95_posterior_cropped(i,1) HPD95_posterior_cropped(i,2)],'color','k')
        end
        title('mean and 95HPD')
    end

    mean_std_posterior_cropped = mean(std_posterior_cropped);
    median_std_posterior_cropped = median(std_posterior_cropped);
    mean_std_posterior_cropped_normalized = mean_std_posterior_cropped/num_spatial_bins;
    median_std_posterior_cropped_normalized = median_std_posterior_cropped/num_spatial_bins;
    mean_CI95_posterior_cropped = mean(abs(CI95_posterior_cropped(:,2)-CI95_posterior_cropped(:,1)));
    median_CI95_posterior_cropped = median(abs(CI95_posterior_cropped(:,2)-CI95_posterior_cropped(:,1)));
    mean_CI95_posterior_cropped_normalized = mean_CI95_posterior_cropped/num_spatial_bins;
    median_CI95_posterior_cropped_normalized = median_CI95_posterior_cropped/num_spatial_bins;
    mean_HPD95_posterior_cropped = mean(abs(HPD95_posterior_cropped(:,2)-HPD95_posterior_cropped(:,1)));
    median_HPD95_posterior_cropped = median(abs(HPD95_posterior_cropped(:,2)-HPD95_posterior_cropped(:,1)));
    mean_HPD95_posterior_cropped_normalized = mean_HPD95_posterior_cropped/num_spatial_bins;
    median_HPD95_posterior_cropped_normalized = median_HPD95_posterior_cropped/num_spatial_bins;

    mean_std_posterior_full = mean(std_posterior_full);
    median_std_posterior_full = median(std_posterior_full);
    mean_std_posterior_full_normalized = mean_std_posterior_full/num_spatial_bins;
    median_std_posterior_full_normalized = median_std_posterior_full/num_spatial_bins;
    mean_CI95_posterior_full = mean(abs(CI95_posterior_full(:,2)-CI95_posterior_full(:,1)));
    median_CI95_posterior_full = median(abs(CI95_posterior_full(:,2)-CI95_posterior_full(:,1)));
    mean_CI95_posterior_full_normalized = mean_CI95_posterior_full/num_spatial_bins;
    median_CI95_posterior_full_normalized = median_CI95_posterior_full/num_spatial_bins;
    mean_HPD95_posterior_full = mean(abs(HPD95_posterior_full(:,2)-HPD95_posterior_full(:,1)));
    median_HPD95_posterior_full = median(abs(HPD95_posterior_full(:,2)-HPD95_posterior_full(:,1)));
    mean_HPD95_posterior_full_normalized = mean_HPD95_posterior_full/num_spatial_bins;
    median_HPD95_posterior_full_normalized = median_HPD95_posterior_full/num_spatial_bins;



    weighted_r = weighted_correlation(best_posterior');
    sharpness = max(posterior_at_peak);
    mean_sharpness_at_peak = mean(posterior_at_peak); % This is the average max posterior in each time bin
    median_sharpness_at_peak = median(posterior_at_peak); % This is the average max posterior in each time bin
    mean_jump_distance_peak = mean(abs(diff(replay_sequence_peak)))/num_spatial_bins; % of track
    median_jump_distance_peak = median(abs(diff(replay_sequence_peak)))/num_spatial_bins; % of track
    max_jump_distance_peak = compute_max_jump_distance(replay_sequence_peak,params.posBinWidth,params.posBinWidth*num_spatial_bins); % of track


    posterior_threshold_for_good_bin = 5*1/(2*num_spatial_bins); %5 times greater than chance
    % Coverage, as in Foster and Wu:
    %For each spatial bin, find the max posterior over the event.
    set_of_positions_covered = find(max(best_posterior) > posterior_threshold_for_good_bin);
    if isempty(set_of_positions_covered)
        coverage_wu = nan;
        coverage_wu_normalized = nan;
        coverage_wu_range = nan;
        coverage_wu_range_normalized = nan;
    else
    coverage_wu = length(set_of_positions_covered);
    coverage_wu_normalized = coverage_wu/num_spatial_bins;
    coverage_wu_range = (max(set_of_positions_covered)-min(set_of_positions_covered));
    coverage_wu_range_normalized = coverage_wu_range/num_spatial_bins;
    end
    % add up the total travel (sign included). In bins
    cum_coverage = abs(sum(diff(replay_sequence_peak)));
    cum_coverage_normalized = cum_coverage/num_spatial_bins;

    posterior_range = max(replay_sequence_peak)-min(replay_sequence_peak);
    posterior_range_normalized = posterior_range/num_spatial_bins;
    best_map_distance_over_time = (posterior_range*params.posBinWidth/100)/replay_duration; % meters/sec
    % replay_score = replay_score_davidson(best_posterior_cropped_normalized,params.davidson_distance_from_line_to_search ,params.posBinWidth,params.decodingWindowShift,1);
    %     % 1 = score, 2 = intercept, 3 = slope, 4 = continuity

    replay_score = [nan nan];

    % check whether the replay starts locally (near the rat's current location)
    [startsLocal, start_distance_from_rat, endsLocal, end_distance_from_rat] = check_replay_for_local_start_or_stop(replay_sequence_peak, ratPos_scaled, 1, params.num_time_bins_to_check_for_local_start, params.distance_thr_for_local_start);
    start_distance_from_rat_normalized = start_distance_from_rat/num_spatial_bins;
    end_distance_from_rat_normalized = end_distance_from_rat/num_spatial_bins;

    start_distance_from_rat_normalized(start_distance_from_rat_normalized > 1) = 1;
    end_distance_from_rat_normalized(end_distance_from_rat_normalized > 1) = 1;


    dispersion = sqrt( nanmean( (replay_sequence_com(:,1)-nanmean(replay_sequence_com(:,1))).^2 ));


    % add forward/reverse designation
    % 1 = forward; 2 = reverse; nan means the replay was
    % static, had a weighted correlation of 0. we should remove
    % these events.
    if (best_map == 1 && weighted_r < 0) || (best_map == 2 && weighted_r > 0)
        direction = 1;
    elseif (best_map == 1 && weighted_r > 0) || (best_map == 2 && weighted_r < 0)
        direction = 2;
    else
        direction = nan;
    end


    % add congruent designation
    if (best_map == 1 && weighted_r < 0 && reward_zone == 2) || (best_map == 2 && weighted_r > 0 && reward_zone == 1)
        forward_congruent_with_rat_location = 1;
        reverse_congruent_with_rat_location = 0;
        congruent_with_rat_location=1;
    elseif (best_map == 1 && weighted_r > 0 && reward_zone == 1) || (best_map == 2 && weighted_r < 0 && reward_zone == 2)
        forward_congruent_with_rat_location = 0;
        reverse_congruent_with_rat_location = 1;
        congruent_with_rat_location=1;
    else
        forward_congruent_with_rat_location = 0;
        reverse_congruent_with_rat_location = 0;
        congruent_with_rat_location=0;
    end

    % add incongruent designation
    if (best_map == 1 && weighted_r < 0 && reward_zone == 1) || (best_map == 2 && weighted_r > 0 && reward_zone == 2)
        forward_incongruent_with_rat_location = 1;
        reverse_incongruent_with_rat_location = 0;
        incongruent_with_rat_location=1;
    elseif (best_map == 1 && weighted_r > 0 && reward_zone == 2) || (best_map == 2 && weighted_r < 0 && reward_zone == 1)
        forward_incongruent_with_rat_location = 0;
        reverse_incongruent_with_rat_location = 1;
        incongruent_with_rat_location=1;
    else
        forward_incongruent_with_rat_location = 0;
        reverse_incongruent_with_rat_location = 0;
        incongruent_with_rat_location=0;
    end





end

eventMetrics.dispersion = dispersion;
eventMetrics.posterior_left = posterior_left;
eventMetrics.posterior_right = posterior_right;
eventMetrics.weighted_r = weighted_r;
eventMetrics.replay_score = replay_score(1);
eventMetrics.replay_slope = replay_score(2);
eventMetrics.com = replay_sequence_com;
eventMetrics.peak = replay_sequence_peak;
eventMetrics.sharpness = sharpness;
eventMetrics.mean_sharpness_at_peak = mean_sharpness_at_peak;
eventMetrics.median_sharpness_at_peak = median_sharpness_at_peak;
eventMetrics.max_jump_distance = max_jump_distance_peak;
eventMetrics.mean_jump_distance = mean_jump_distance_peak;
eventMetrics.median_jump_distance = median_jump_distance_peak;
eventMetrics.total_posterior_left = total_posterior_left;
eventMetrics.total_posterior_right = total_posterior_right;

eventMetrics.best_map = best_map;
eventMetrics.range = posterior_range;
eventMetrics.range_normalized = posterior_range_normalized;
eventMetrics.range_over_time = best_map_distance_over_time;
eventMetrics.coverage_wu = coverage_wu;
eventMetrics.coverage_wu_normalized = coverage_wu_normalized;
eventMetrics.coverage_wu_range = coverage_wu_range;
eventMetrics.coverage_wu_range_normalized = coverage_wu_range_normalized;
eventMetrics.cum_coverage = cum_coverage;
eventMetrics.cum_coverage_normalized = cum_coverage_normalized;

eventMetrics.percent_posterior_in_left_map = percent_posterior_in_left_map;
eventMetrics.percent_posterior_in_right_map = percent_posterior_in_right_map;
eventMetrics.percent_non_local_posterior_in_left_map = percent_non_local_posterior_in_left_map;
eventMetrics.percent_non_local_posterior_in_right_map = percent_non_local_posterior_in_right_map;
eventMetrics.percent_posterior_ends_masked_in_left_map = percent_posterior_ends_masked_in_left_map;
eventMetrics.percent_posterior_ends_masked_in_right_map = percent_posterior_ends_masked_in_right_map;

eventMetrics.percent_posterior_in_left_map_full = percent_posterior_in_left_map_full;
eventMetrics.percent_posterior_in_right_map_full = percent_posterior_in_right_map_full;
eventMetrics.percent_non_local_posterior_in_left_map_full = percent_non_local_posterior_in_left_map_full;
eventMetrics.percent_non_local_posterior_in_right_map_full = percent_non_local_posterior_in_right_map_full;
eventMetrics.percent_posterior_ends_masked_in_left_map_full = percent_posterior_ends_masked_in_left_map_full;
eventMetrics.percent_posterior_ends_masked_in_right_map_full = percent_posterior_ends_masked_in_right_map_full;







eventMetrics.startsLocal = startsLocal;
eventMetrics.start_distance_from_rat = start_distance_from_rat;
eventMetrics.start_distance_from_rat_normalized = start_distance_from_rat_normalized;
eventMetrics.endsLocal = endsLocal;
eventMetrics.end_distance_from_rat = end_distance_from_rat;
eventMetrics.end_distance_from_rat_normalized = end_distance_from_rat_normalized; 
eventMetrics.direction = direction; %1=forward; 2=reverse;
eventMetrics.forward_congruent_with_rat_location = forward_congruent_with_rat_location;
eventMetrics.reverse_congruent_with_rat_location = reverse_congruent_with_rat_location;
eventMetrics.congruent_with_rat_location=congruent_with_rat_location;
eventMetrics.forward_incongruent_with_rat_location = forward_incongruent_with_rat_location;
eventMetrics.reverse_incongruent_with_rat_location = reverse_incongruent_with_rat_location;
eventMetrics.incongruent_with_rat_location=incongruent_with_rat_location;
eventMetrics.posterior_spread_full = posterior_spread_full;
eventMetrics.posterior_spread_cropped = posterior_spread_cropped;

eventMetrics.mean_std_posterior_cropped = mean_std_posterior_cropped;
eventMetrics.median_std_posterior_cropped =median_std_posterior_cropped;
eventMetrics.mean_std_posterior_cropped_normalized = mean_std_posterior_cropped_normalized;
eventMetrics.median_std_posterior_cropped_normalized = median_std_posterior_cropped_normalized;
eventMetrics.mean_CI95_posterior_cropped = mean_CI95_posterior_cropped;
eventMetrics.median_CI95_posterior_cropped = median_CI95_posterior_cropped;
eventMetrics.mean_CI95_posterior_cropped_normalized = mean_CI95_posterior_cropped_normalized;
eventMetrics.median_CI95_posterior_cropped_normalized = median_CI95_posterior_cropped_normalized;
eventMetrics.mean_HPD95_posterior_cropped = mean_HPD95_posterior_cropped;
eventMetrics.median_HPD95_posterior_cropped = median_HPD95_posterior_cropped;
eventMetrics.mean_HPD95_posterior_cropped_normalized = mean_HPD95_posterior_cropped_normalized;
eventMetrics.median_HPD95_posterior_cropped_normalized = median_HPD95_posterior_cropped_normalized;


eventMetrics.mean_std_posterior_full = mean_std_posterior_full;
eventMetrics.median_std_posterior_full = median_std_posterior_full;
eventMetrics.mean_std_posterior_full_normalized = mean_std_posterior_full_normalized;
eventMetrics.median_std_posterior_full_normalized = median_std_posterior_full_normalized;
eventMetrics.mean_CI95_posterior_full = mean_CI95_posterior_full;
eventMetrics.median_CI95_posterior_full = median_CI95_posterior_full;
eventMetrics.mean_CI95_posterior_full_normalized = mean_CI95_posterior_full_normalized;
eventMetrics.median_CI95_posterior_full_normalized = median_CI95_posterior_full_normalized;
eventMetrics.mean_HPD95_posterior_full = mean_HPD95_posterior_full;
eventMetrics.median_HPD95_posterior_full = median_HPD95_posterior_full;
eventMetrics.mean_HPD95_posterior_full_normalized = mean_HPD95_posterior_full_normalized;
eventMetrics.median_HPD95_posterior_full_normalized = median_HPD95_posterior_full_normalized;
