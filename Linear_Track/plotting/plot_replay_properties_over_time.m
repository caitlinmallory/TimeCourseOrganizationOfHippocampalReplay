test_type = 'mean';
plot_all_individual_properties = 0;
control_experimental_animals_to_include = [0 1 2];

t_congruent_replay = t_replay(t_replay.congruent_with_rat_location ==1 & ismember(t_replay.control_experimental_flag,control_experimental_animals_to_include),:);
t_forward_congruent_replay = t_replay(t_replay.direction==1 & t_replay.congruent_with_rat_location ==1 & ismember(t_replay.control_experimental_flag,control_experimental_animals_to_include),:);
t_reverse_congruent_replay = t_replay(t_replay.direction==2 & t_replay.congruent_with_rat_location ==1  & ismember(t_replay.control_experimental_flag,control_experimental_animals_to_include),:);

t_forward_replay = t_replay(t_replay.direction==1 & ismember(t_replay.control_experimental_flag,control_experimental_animals_to_include),:);
t_reverse_replay = t_replay(t_replay.direction==2 & ismember(t_replay.control_experimental_flag,control_experimental_animals_to_include),:);

t_future_replay = t_replay(t_replay.future_map_replay==1 & ismember(t_replay.control_experimental_flag,control_experimental_animals_to_include),:);
t_past_replay = t_replay(t_replay.past_map_replay==1 & ismember(t_replay.control_experimental_flag,control_experimental_animals_to_include),:);

events = {t_forward_congruent_replay,t_reverse_congruent_replay,t_congruent_replay};
event_types = {'forward_congruent_replays','reverse_congruent_replays','t_congruent_replay'};

plot_laser_on_data = 0; % set to one if you want to show both laser on and off data

windowSize = 2;
windowShift = 0.5;
start_time = 0;
end_time = 14;

t_forward_congruent_replay(t_forward_congruent_replay.laser_state==0 & t_forward_congruent_replay.time_since_reward_zone_entry<10,:)
t_reverse_congruent_replay(t_reverse_congruent_replay.laser_state==0 & t_reverse_congruent_replay.time_since_reward_zone_entry<10,:)

bin_start = (start_time:windowShift:(end_time-windowSize))';
bin_end = bin_start + windowSize;

timeBins = [bin_start bin_end];
bin_centers = mean(timeBins,2);

load('/home/caitlin/Data/Processed_Data/properties_table');

properties_to_plot = {...
 'fraction_of_excitatory_cells_participating';...
 'range';...
 'mean_HPD95_posterior_cropped_normalized'};

for event = 1:length(event_types)
    for i = 1:length(timeBins)
        for property = 1:length(properties_to_plot)

            properties_table_row = find(strcmp(properties_to_plot{property},properties.names)==1);

            events_sub = events{event};
            % Look at the property during all candidate events laser off or laser on
            inds = find(events_sub.time_since_reward_zone_entry>timeBins(i,1) & events_sub.time_since_reward_zone_entry <= timeBins(i,2) );
            inds_0 = find(events_sub.time_since_reward_zone_entry>timeBins(i,1) & events_sub.time_since_reward_zone_entry <= timeBins(i,2) & events_sub.laser_state == 0);
            inds_1 = find(events_sub.time_since_reward_zone_entry>timeBins(i,1) & events_sub.time_since_reward_zone_entry <= timeBins(i,2) & events_sub.laser_state == 1);

            data_sub = events_sub{inds,properties_to_plot{property}}; data_sub(isnan(data_sub)) = [];
            data_sub_0 = events_sub{inds_0,properties_to_plot{property}}; data_sub_0(isnan(data_sub_0)) = [];
            data_sub_1 = events_sub{inds_1,properties_to_plot{property}}; data_sub_1(isnan(data_sub_1)) = [];

            properties_over_time.(event_types{event}).(properties_to_plot{property}).group_summary_n(i) = sum(~isnan(data_sub));
            properties_over_time.(event_types{event}).(properties_to_plot{property}).off_summary_n(i) = sum(~isnan(data_sub_0));
            properties_over_time.(event_types{event}).(properties_to_plot{property}).on_summary_n(i) = sum(~isnan(data_sub_1));
            properties_over_time.(event_types{event}).(properties_to_plot{property}).group_summary(i) = group_summary_test(data_sub,test_type);
            properties_over_time.(event_types{event}).(properties_to_plot{property}).off_summary(i) = group_summary_test(data_sub_0,test_type);
            properties_over_time.(event_types{event}).(properties_to_plot{property}).on_summary(i) = group_summary_test(data_sub_1,test_type);
            properties_over_time.(event_types{event}).(properties_to_plot{property}).group_err(i,:) = conf_test(data_sub,test_type);
            properties_over_time.(event_types{event}).(properties_to_plot{property}).off_err(i,:) = conf_test(data_sub_0,test_type);
            properties_over_time.(event_types{event}).(properties_to_plot{property}).on_err(i,:) = conf_test(data_sub_1,test_type);
            properties_over_time.(event_types{event}).(properties_to_plot{property}).pval_on_off(i) = stat_test(data_sub_0,data_sub_1,test_type);
        end
    end
end

% Compute the forwards-reverse difference in each time bin, as well as the 95% confidence interval of that difference
for property = 1:length(properties_to_plot)
    for i = 1:length(timeBins)
        inds_f0 = find(events{1}.time_since_reward_zone_entry>timeBins(i,1) & events{1}.time_since_reward_zone_entry <= timeBins(i,2) & events{1}.laser_state ==0);
        inds_r0 = find(events{2}.time_since_reward_zone_entry>timeBins(i,1) & events{2}.time_since_reward_zone_entry <= timeBins(i,2) & events{2}.laser_state ==0);

        data_sub_f0 = events{1}.(properties_to_plot{property})(inds_f0); data_sub_f0(isnan(data_sub_f0)) = [];
        data_sub_r0 = events{2}.(properties_to_plot{property})(inds_r0); data_sub_r0(isnan(data_sub_r0)) = [];
        [properties_over_time.forward_minus_reverse_replays.(properties_to_plot{property}).for_rev_pval(i),~] = ranksum(data_sub_f0,data_sub_r0);
    end
end

%% Make a nice publication figure-- reverse versus forward replays (Fig 2E)
use_constant_ylims = 0;
figure('Position', [2593 688 90 275]);
tiledlayout(3,1)

properties_to_plot = {'fraction_of_excitatory_cells_participating','range','mean_HPD95_posterior_cropped_normalized'};
forward_color = [.4660 0.6740 0.1880];
reverse_color = [0.4940 0.1840 0.5560];
transparancy_off = 0.2;
transparancy_on = 0.1;

for property = 1:length(properties_to_plot)

    properties_table_row = find(strcmp(properties_to_plot{property},properties.names)==1);

    % plot forward data over the stopping period
    ax = nexttile(property);
    data = properties_over_time.(event_types{1}).(properties_to_plot{property});
    x = data.off_summary';
    err_low = x-data.off_err(:,1);
    err_high = data.off_err(:,2)-x;
    H=shadedErrorBar(bin_centers,x,[err_low err_high],'lineprops','k'); hold on
    H.patch.FaceColor = forward_color;
    H.mainLine.Color = forward_color;
    H.patch.FaceAlpha = transparancy_off;
    H.edge(1).Color = 'none'; H.edge(2).Color = 'none';
    H.mainLine.LineWidth = 1;

    if use_constant_ylims == 1
        ylim(properties.ylims{properties_table_row}(1,1:2))
    end
    if strcmp(properties_to_plot{property},'weighted_r')
        ylim([0.7 0.88]);
    elseif strcmp(properties_to_plot{property},'mean_HPD95_posterior_cropped_normalized')
        ylim([15 25])
    elseif strcmp(properties_to_plot{property},'range')
        ylim([75 180])
        yticks([100,180]);
        yticklabels({'10' '18'})
    elseif strcmp(properties_to_plot{property},'fraction_of_excitatory_cells_participating')
    end
    a = events{1}.time_since_reward_zone_entry(events{1}.laser_state==0 & events{1}.time_since_reward_zone_entry <= 10);
    b = events{1}.(properties_to_plot{property})(events{1}.laser_state==0 & events{1}.time_since_reward_zone_entry <= 10);
    [rho_unbinned_for,p_unbinned_for] = nancorr(a,b);

    % plot reverse data over the stopping period
    data = properties_over_time.(event_types{2}).(properties_to_plot{property});
    x = data.off_summary';
    err_low = x-data.off_err(:,1);
    err_high = data.off_err(:,2)-x;
    H=shadedErrorBar(bin_centers,x,[err_low err_high],'lineprops','k');
    H.patch.FaceColor = reverse_color;
    H.mainLine.Color = reverse_color;
    H.patch.FaceAlpha = transparancy_off;
    H.edge(1).Color = 'none'; H.edge(2).Color = 'none';
    H.mainLine.LineWidth = 1;

    if use_constant_ylims == 1
        ylim(properties.ylims{properties_table_row}(1,1:2))
    end

    a = events{2}.time_since_reward_zone_entry(events{2}.laser_state==0 & events{2}.time_since_reward_zone_entry <= 10);
    b = events{2}.(properties_to_plot{property})(events{2}.laser_state==0 & events{2}.time_since_reward_zone_entry <= 10);
    [rho_unbinned_rev,p_unbinned_rev] = nancorr(a,b);

    data = properties_over_time.forward_minus_reverse_replays.(properties_to_plot{property});
    ylimit = gca().YLim;

    plot(bin_centers(data.for_rev_pval<0.05),repmat(ylimit(2),[sum(data.for_rev_pval<0.05),1]),'.k')
    xtickangle(0);
    xticks(0:2:10);
    if property==3
        xlabel('Time (s)')
    end

    xlim([0 10])
    xticks([0:2:10]);
    xtickangle(0);
    yticklabels({})
end

set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures';
saveas(gcf,fullfile(fig_path,'summary_for_rev_overlaid'),'jpg')
saveas(gcf,fullfile(fig_path,'summary_for_rev_overlaid'),'pdf')

for_n = [min(properties_over_time.(event_types{1}).(properties_to_plot{property}).off_summary_n(bin_centers<=10)) ...
    max(properties_over_time.(event_types{1}).(properties_to_plot{property}).off_summary_n(bin_centers<=10))]

rev_n = [min(properties_over_time.(event_types{2}).(properties_to_plot{property}).off_summary_n(bin_centers<=10)) ...
    max(properties_over_time.(event_types{2}).(properties_to_plot{property}).off_summary_n(bin_centers<=10))];

function summary_x = group_summary_test(x,test_type)
if strcmp(test_type,'median')==1
    summary_x = median(x);
end
if strcmp(test_type,'mean')==1
    summary_x = mean(x);
end
end

function CI = conf_test(x,test_type)
if strcmp(test_type,'median')==1
    CI = bootci(1000,@(x)median(x),x);
end
if strcmp(test_type,'mean')==1
    err = std(x)/(sqrt(length(x)));
    CI = [mean(x)-err mean(x) + err];
end
end

function p = stat_test(x,y,test_type)
% make sure there is data in both x and y
if ~isempty(x) && ~isempty(y)
    [p] = ranksum(x,y);
else
    p = nan;
end
end

function CI = bootstrap_differences(data1,data2,test_type)
% bootstrap:
boot_diffs = nan(1000,1);
for i = 1:1000
    d1 = randsample(data1,length(data1),'true');
    d2 = randsample(data2,length(data2),'true');
    if strcmp(test_type,'mean')
        boot_diffs(i) = mean(d1)-mean(d2);
    elseif strcmp(test_type,'median')
        boot_diffs(i) = median(d1)-median(d2);
    end
end
CI = [prctile(boot_diffs,2.5) prctile(boot_diffs,97.5)];
end