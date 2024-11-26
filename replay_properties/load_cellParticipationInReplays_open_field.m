function cell_participation_metrics = load_cellParticipationInReplays_open_field(params,clusters,event_timePoints)

%% Criterion for cells:
% params.restrict_to_well_isolated_cells = 0;
% params.num_fields_thr = 0;
% params.min_peak_fr_thr = 0;
keyboard
cell_participation_metrics = struct();
    for i = 1:length(clusters)
        clusters(i).peak_fr = clusters(i).runs(1).maxSpatialFiringRate;
        clusters(i).num_subFields = clusters(i).runs(1).num_subFields;
    end


    if params.restrict_to_well_isolated_cells == 0
        clusters(i).meets_well_isolated_requirement = 1;
    else
        clusters(i).meets_well_isolated_requirement = clusters(i).well_isolated;
    end

    excitatory = [clusters.Excitatory]' == 1 & [clusters.meets_well_isolated_requirement]' == 1  & [clusters.firing_rate]' > params.min_peak_fr_thr;
    inhibitory = [clusters.Excitatory]' == 0 & [clusters.meets_well_isolated_requirement]' == 1  & [clusters.firing_rate]' > params.min_peak_fr_thr;
%     unimodal = [clusters.Modality]' == 1 & excitatory & [clusters.num_subFields]' >= params.num_fields_thr;
%     bimodal = [clusters.Modality]' == 2 & excitatory >= params.num_fields_thr;

    %%
  %  num_spikes = load_numSpks_timeBins(event_timeBins,clusters,(1:length(clusters)),params.decodingWindowShift*params.spikeSampRate,params.decodingWindowSize*params.spikeSampRate);
      cell_participation_metrics.num_spikes_per_decoding_bin = mean(sum(num_spikes));
      num_spikes = nan(length(clusters),1);
     for i = 1:length(clusters)
           num_spikes(i) = compute_dataTemporalConcatenation(clusters(i).spkTime,event_timePoints)
     end
    cell_participation_metrics.num_spikes = sum(sum(num_spikes)); % total number of spikes
    cell_participation_metrics.num_excitatory_spikes = sum(sum(num_spikes(excitatory,:)));
    cell_participation_metrics.num_inhibitory_spikes = sum(sum(num_spikes(inhibitory,:)));
    cell_participation_metrics.ratio_excitatory_inhibitory_spikes = cell_participation_metrics.num_excitatory_spikes/(cell_participation_metrics.num_excitatory_spikes + cell_participation_metrics.num_inhibitory_spikes);

    % find the percentage of cells that fired at least 1 spike in the replay:
    num_spikes_sum = sum(num_spikes,2); % total num of spikes emitted by each cell
    cell_participation_metrics.fraction_of_cells_participating = sum(num_spikes_sum>= params.min_num_spikes_to_participate_in_replay)/size(num_spikes,1);

    % add in the num of cells that fired at least 1
    % spike in the replay:
    cell_participation_metrics.num_of_cells_participating = sum(num_spikes_sum>params.min_num_spikes_to_participate_in_replay);

    % add in the average num of spikes emiited by each
    % participating excitatory cell, and the average firing
    % rate (spikes emitted/replay length)
    cell_participation_metrics.num_of_excitatory_cells_participating = sum(num_spikes_sum(excitatory)>params.min_num_spikes_to_participate_in_replay);
    cell_participation_metrics.mean_num_spikes_all_excitatory_cells = nanmean(num_spikes_sum(excitatory));
    cell_participation_metrics.mean_fr_all_excitatory_cells = nanmean(num_spikes_sum(excitatory))/(size(num_spikes,2)*params.decodingWindowShift);
    cell_participation_metrics.num_of_inhibitory_cells_participating = sum(num_spikes_sum(inhibitory)>params.min_num_spikes_to_participate_in_replay);
    cell_participation_metrics.mean_num_spikes_all_inhibitory_cells = nanmean(num_spikes_sum(inhibitory));
    cell_participation_metrics.mean_fr_all_inhibitory_cells = nanmean(num_spikes_sum(inhibitory))/(size(num_spikes,2)*params.decodingWindowShift);
    cell_participation_metrics.mean_num_spikes_all_bimodal_cells = nanmean(num_spikes_sum(bimodal));
    cell_participation_metrics.mean_fr_all_bimodal_cells = nanmean(num_spikes_sum(bimodal))/(size(num_spikes,2)*params.decodingWindowShift);
    cell_participation_metrics.mean_num_spikes_all_unimodal_cells = nanmean(num_spikes_sum(unimodal));
    cell_participation_metrics.mean_fr_all_unimodal_cells = nanmean(num_spikes_sum(unimodal))/(size(num_spikes,2)*params.decodingWindowShift);


    cell_participation_metrics.mean_num_spikes_participating_excitatory_cells = nanmean(num_spikes_sum(num_spikes_sum>0 & excitatory));
    cell_participation_metrics.mean_fr_participating_excitatory_cells = nanmean(num_spikes_sum(num_spikes_sum>0 & excitatory))/(size(num_spikes,2)*params.decodingWindowShift);
    cell_participation_metrics.mean_num_spikes_participating_inhibitory_cells = nanmean(num_spikes_sum(num_spikes_sum>0 & inhibitory));
    cell_participation_metrics.mean_fr_participating_inhibitory_cells = nanmean(num_spikes_sum(num_spikes_sum>0 & inhibitory))/(size(num_spikes,2)*params.decodingWindowShift);
    cell_participation_metrics.mean_num_spikes_participating_bimodal_cells = nanmean(num_spikes_sum(num_spikes_sum>0 & bimodal));
    cell_participation_metrics.mean_fr_participating_bimodal_cells = nanmean(num_spikes_sum(num_spikes_sum>0 & bimodal))/(size(num_spikes,2)*params.decodingWindowShift);
    cell_participation_metrics.mean_num_spikes_participating_unimodal_cells = nanmean(num_spikes_sum(num_spikes_sum>0 & unimodal));
    cell_participation_metrics.mean_fr_participating_unimodal_cells = nanmean(num_spikes_sum(num_spikes_sum >0 & unimodal))/(size(num_spikes,2)*params.decodingWindowShift);



    % add in the num of cells in the session
    cell_participation_metrics.total_num_cells_in_session = size(num_spikes,1);

    participating_clusters = num_spikes_sum > 0;
    participating_bimodal_cells = participating_clusters(bimodal);
    participating_unimodal_cells = participating_clusters(unimodal);
    participating_excitatory_cells = participating_clusters(excitatory);
    participating_inhibitory_cells = participating_clusters(inhibitory);

    cell_participation_metrics.num_of_excitatory_cells_in_session = sum(excitatory);
    if sum(excitatory)==0
        cell_participation_metrics.fraction_of_excitatory_cells_participating = nan;
        cell_participation_metrics.num_of_excitatory_cells_participating = nan;
    else
        cell_participation_metrics.fraction_of_excitatory_cells_participating = sum(participating_excitatory_cells)/sum(excitatory);
        cell_participation_metrics.num_of_excitatory_cells_participating = sum(participating_excitatory_cells);
    end

    cell_participation_metrics.num_of_bimodal_cells_in_session = sum(bimodal);
    if sum(bimodal)==0
        cell_participation_metrics.fraction_of_bimodal_cells_participating = nan;
    else
        cell_participation_metrics.fraction_of_bimodal_cells_participating = sum(participating_bimodal_cells)/sum(bimodal);
    end

    cell_participation_metrics.num_of_unimodal_cells_in_session = sum(unimodal);
    if sum(unimodal)==0
        cell_participation_metrics.fraction_of_unimodal_cells_participating = nan;
    else
        cell_participation_metrics.fraction_of_unimodal_cells_participating = sum(participating_unimodal_cells)/sum(unimodal);
    end

    cell_participation_metrics.num_of_inhibitory_cells_in_session = sum(inhibitory);
    if sum(inhibitory)==0
        cell_participation_metrics.fraction_of_inhibitory_cells_participating = nan;
    else
        cell_participation_metrics.fraction_of_inhibitory_cells_participating = sum(participating_inhibitory_cells)/sum(inhibitory);
    end

else
    cell_participation_metrics.num_spikes = nan;
    cell_participation_metrics.num_spikes_per_decoding_bin = nan;
    cell_participation_metrics.num_excitatory_spikes = nan;
    cell_participation_metrics.num_inhibitory_spikes = nan;
    cell_participation_metrics.ratio_excitatory_inhibitory_spikes = nan;


    cell_participation_metrics.fraction_of_cells_participating = nan;
    cell_participation_metrics.num_of_cells_participating = nan;
    cell_participation_metrics.num_of_excitatory_cells_participating = nan;
    cell_participation_metrics.mean_num_spikes_all_excitatory_cells = nan;
    cell_participation_metrics.mean_fr_all_excitatory_cells = nan;
    cell_participation_metrics.num_of_inhibitory_cells_participating = nan;
    cell_participation_metrics.mean_num_spikes_all_inhibitory_cells = nan;
    cell_participation_metrics.mean_fr_all_inhibitory_cells = nan;
    cell_participation_metrics.mean_num_spikes_all_bimodal_cells = nan;
    cell_participation_metrics.mean_fr_all_bimodal_cells = nan;
    cell_participation_metrics.mean_num_spikes_all_unimodal_cells = nan;
    cell_participation_metrics.mean_fr_all_unimodal_cells = nan;
    cell_participation_metrics.mean_num_spikes_participating_excitatory_cells = nan;
    cell_participation_metrics.mean_fr_participating_excitatory_cells = nan;
    cell_participation_metrics.mean_num_spikes_participating_inhibitory_cells = nan;
    cell_participation_metrics.mean_fr_participating_inhibitory_cells = nan;
    cell_participation_metrics.mean_num_spikes_participating_bimodal_cells = nan;
    cell_participation_metrics.mean_fr_participating_bimodal_cells = nan;
    cell_participation_metrics.mean_num_spikes_participating_unimodal_cells = nan;
    cell_participation_metrics.mean_fr_participating_unimodal_cells = nan;
    cell_participation_metrics.total_num_cells_in_session = nan;
    cell_participation_metrics.num_of_excitatory_cells_in_session = nan;
    cell_participation_metrics.fraction_of_excitatory_cells_participating = nan;
    cell_participation_metrics.num_of_excitatory_cells_participating = nan;
    cell_participation_metrics.num_of_bimodal_cells_in_session = nan;
    cell_participation_metrics.fraction_of_bimodal_cells_participating = nan;
    cell_participation_metrics.num_of_unimodal_cells_in_session = nan;
    cell_participation_metrics.fraction_of_unimodal_cells_participating = nan;
    cell_participation_metrics.num_of_inhibitory_cells_in_session = nan;
    cell_participation_metrics.fraction_of_inhibitory_cells_participating = nan;

end
