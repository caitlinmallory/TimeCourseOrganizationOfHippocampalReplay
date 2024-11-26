function event_metrics = load_replay_ripple_power_spike_density(full_event_time_points,cropped_event_time_points,zscored_LFP,zscored_LFP_15,zscored_spikeDensity_pyr_125,zscored_spikeDensity_pyr_15,zscored_spikeDensity_inhib_15)


plot_fig = 0;
lfpSampRate = 1500;
event_metrics = table();
traces_to_examine = struct();
traces_to_examine.zscored_ripple_power_15 = [zscored_LFP_15(:,1) zscored_LFP_15(:,4) zscored_LFP_15(:,7)];
traces_to_examine.zscored_ripple_power = [zscored_LFP(:,1) zscored_LFP(:,4) zscored_LFP(:,7)]; % [time filtered_lfp amplitude]
traces_to_examine.zscored_sd_pyr_125 = [zscored_spikeDensity_pyr_125(:,1) nan(length(zscored_spikeDensity_pyr_125),1) zscored_spikeDensity_pyr_125(:,3)]; % [time nan amplitude]
traces_to_examine.zscored_sd_pyr_15 = [zscored_spikeDensity_pyr_15(:,1) nan(length(zscored_spikeDensity_pyr_15),1) zscored_spikeDensity_pyr_15(:,3)]; % [time nan amplitude]
traces_to_examine.zscored_sd_inhib_15 = [zscored_spikeDensity_inhib_15(:,1) nan(length(zscored_spikeDensity_inhib_15),1) zscored_spikeDensity_inhib_15(:,3)]; % [time nan amplitude]


%% 1) divide the full event into 3 parts (fraction of event) or 20 ms time bins
% first entry = first third of the event, second entry = second third of
% the event, third entry = third third of the event, fourth entry = entire
% event
time_bin_edges= linspace(full_event_time_points(1),full_event_time_points(2),4)';
full_time_bins_fraction = [time_bin_edges(1:end-1) time_bin_edges(2:end)];
full_time_bins_fraction = [full_time_bins_fraction; [full_event_time_points(1) full_event_time_points(2)]];

% 2) Also save from -0.1 to 1 second of data surrounding the start of the
% event.
min_time = full_event_time_points(1)-30000.*0.1;
max_time = full_event_time_points(1)+30000.*1;
full_time_bins = [min_time max_time];


%% 1) divide the cropped event into 3 parts (fraction of event) or 20 ms time bins
if ~isempty(cropped_event_time_points)
    time_bin_edges= linspace(full_event_time_points(1),full_event_time_points(2),4)';
    cropped_time_bins_fraction = [time_bin_edges(1:end-1) time_bin_edges(2:end)];
    cropped_time_bins_fraction = [cropped_time_bins_fraction; [cropped_event_time_points(1) cropped_event_time_points(2)]];

    min_time = cropped_event_time_points(1)-30000.*0.1;
    max_time = cropped_event_time_points(1)+30000.*1;
    cropped_time_bins = [min_time max_time];
    % full_time_bins = [time_bin_edges(1:end-1) time_bin_edges(2:end)];
else
    cropped_time_bins_fraction = nan(4,2);
    cropped_time_bins = nan(length(full_time_bins),2);
end

trace_names = fieldnames(traces_to_examine);
event_types = {'full_event'; 'cropped_event'};
time_bins = [{full_time_bins_fraction} {full_time_bins}; {cropped_time_bins_fraction} {cropped_time_bins}];

for i = 1:length(trace_names)
    trace_name = trace_names{i};

    data_tbl = table();

    for j = 1 % full or cropped
        
        data_tbl_sub = table();
        data_tbl_sub.('mean_power_over_fraction_of_event') = nan(1,4);
        data_tbl_sub.('peak_power_over_fraction_of_event') = nan(1,4);
        data_tbl_sub.('peak_power_over_time_in_event') = nan(1,round(lfpSampRate*((full_time_bins(2)-full_time_bins(1))/30000)));        
        event_type = event_types{j};


        times_of_interest = time_bins{j,1}; % fractional bins
        if unique(isnan(times_of_interest)) ~= 1
            for bin = 1:length(times_of_interest)
                x = compute_dataTemporalConcatenation(traces_to_examine.(trace_name),times_of_interest(bin,:));
                %inds = traces_to_examine.(trace_name)(:,1)>=times_of_interest(bin,1) & traces_to_examine.(trace_name)(:,1) < times_of_interest(bin,2);
                if ~isempty(x)
                    data_tbl_sub.('mean_power_over_fraction_of_event')(bin) = nanmean(x(:,3));
                    data_tbl_sub.('peak_power_over_fraction_of_event')(bin) = nanmax(x(:,3));
                end
            end
        end

        times_of_interest = time_bins{j,2}; % time bins
        if unique(isnan(times_of_interest)) ~= 1
                x = compute_dataTemporalConcatenation(traces_to_examine.(trace_name),times_of_interest);
                if length(x) ~= round(lfpSampRate*((full_time_bins(2)-full_time_bins(1))/30000)) % rare case where there isn't enough data before or after the event bounds
                data_tbl_sub.('peak_power_over_time_in_event') = nan(1,round(lfpSampRate*((full_time_bins(2)-full_time_bins(1))/30000)));        
                else
                data_tbl_sub.('peak_power_over_time_in_event') = x(:,3)';
                end
        end

        data_tbl.(event_type) = data_tbl_sub;
        event_metrics.(trace_name) = data_tbl_sub;

    end
end


if plot_fig == 1
    figure('Position',[896 297 792 648])
    plot(zscored_LFP_sub(:,1)./30000 - zscored_LFP_sub(1,1)./30000,zscored_LFP_sub(:,3).*(1/10),'color',[0 0 0]); hold on; % filtered ripple
    plot(zscored_LFP_sub(:,1)./30000  - zscored_LFP_sub(1,1)./30000,zscored_LFP_sub(:,end))

    plot(zscored_spikeDensity_sub_pyr(:,1)./30000  - zscored_LFP_sub(1,1)./30000,zscored_spikeDensity_sub_pyr(:,end)),'y'; hold on
    plot(zscored_spikeDensity_pyr_15(:,1)./30000  - zscored_LFP_sub(1,1)./30000,zscored_spikeDensity_pyr_15(:,3),'r')
    legend({'filtered lfp','z-scored R.P.','z-scored SD','z-scored SD less smoothing'})
    drawnow
    keyboard
end

end
