function [corr_value, p_value] = spearman_median(clusters,event_times,sort_idx_fields,plot_figure)
median_spikeTime = [];
cell_IDs = [];
c = 1;
for cell = 1:length(sort_idx_fields)
    cell_num = sort_idx_fields(cell);
    cell_spikes = compute_dataTemporalConcatenation(clusters(cell_num).spkTime,event_times);
    if ~isempty(cell_spikes)
    median_spikeTime(c) = median(cell_spikes);
    cell_IDs(c) = cell;
    c = c+1;
    end
end

 if isempty(cell_IDs) || isempty(median_spikeTime)
     corr_value  = NaN;
     p_value = NaN;
 else
     
 [corr_value, p_value] = corr(cell_IDs',median_spikeTime','type','Spearman');
 % corr_value = abs(corr_value);  %just care about magnitude of correlation
 end

 if plot_figure==1
     plot(cell_IDs,median_spikeTime,'ok')
 end
