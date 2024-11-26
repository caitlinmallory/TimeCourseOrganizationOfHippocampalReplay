set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
sig_test = 'medians';

properties = {'dispersion','meanAngDisplacement_futPath','meanAngDisplacement_pastPath', 'percent_participation_excitatory','start_to_end_distance','posteriorSpread'};
plot_all_properties = 0;
%%
windowSize = 2;
windowShift = 0.5;
start_time = 0;
end_time = 10;
bin_start = (start_time:windowShift:(end_time-windowSize))';
bin_end = bin_start + windowSize;
bin_edges = [bin_start bin_end];
bin_centers = mean(bin_edges,2);

binned_replay_properties = struct();
drink_times = [0 10];

replay_sub = replay(replay.time_since_real_drink_onset>=drink_times(1) & replay.time_since_real_drink_onset <= drink_times(2),:);
replay_sub.ripple_power(replay_sub.ripple_power==-inf) = nan; % Occasionally a session didn't have good LFP- but we want to include these anyway.

% For Figure 2L, use the following (but then only use laser off data)
group1 = 'off';
group2 = 'on';
replay_sub_group1 = replay_sub(replay_sub.laser_state_binary==0,:);
replay_sub_group2 = replay_sub(replay_sub.laser_state_binary==1,:);
group_1_n_total = height(replay_sub_group1);
group_2_n_total = height(replay_sub_group2);

for i = 1:length(bin_edges)
    replay_sub_sub_group1 = replay_sub_group1(replay_sub_group1.time_since_real_drink_onset>=bin_edges(i,1) & replay_sub_group1.time_since_real_drink_onset<=bin_edges(i,2),:);
    replay_sub_sub_group2 = replay_sub_group2(replay_sub_group2.time_since_real_drink_onset>=bin_edges(i,1) & replay_sub_group2.time_since_real_drink_onset<=bin_edges(i,2),:);
    retro_group1_sub = replay_sub_sub_group1(replay_sub_sub_group1.past == 1,:);
    retro_group2_sub  = replay_sub_sub_group2(replay_sub_sub_group2.past == 1,:);
    prospective_group1_sub = replay_sub_sub_group1(replay_sub_sub_group1.future == 1,:);
    prospective_group2_sub  = replay_sub_sub_group2(replay_sub_sub_group2.future == 1,:);

    for j = 1:length(properties)

        binned_replay_properties.(properties{j}).replay_sub_sub_group1.data{i} = replay_sub_sub_group1.(properties{j});
        binned_replay_properties.(properties{j}).replay_sub_sub_group1.mean(i) = nanmean(replay_sub_sub_group1.(properties{j}));
        binned_replay_properties.(properties{j}).replay_sub_sub_group1.sem(i) = nanstd(replay_sub_sub_group1.(properties{j}))./sqrt(sum(~isnan(replay_sub_sub_group1.(properties{j}))));
        binned_replay_properties.(properties{j}).replay_sub_sub_group1.n(i) = sum(~isnan(replay_sub_sub_group1.(properties{j})));

        binned_replay_properties.(properties{j}).replay_sub_sub_group2.data{i} = replay_sub_sub_group2.(properties{j});
        binned_replay_properties.(properties{j}).replay_sub_sub_group2.mean(i) = nanmean(replay_sub_sub_group2.(properties{j}));
        binned_replay_properties.(properties{j}).replay_sub_sub_group2.sem(i) = nanstd(replay_sub_sub_group2.(properties{j}))./sqrt(sum(~isnan(replay_sub_sub_group2.(properties{j}))));
        binned_replay_properties.(properties{j}).replay_sub_sub_group2.n(i) = sum(~isnan(replay_sub_sub_group2.(properties{j})));

        binned_replay_properties.(properties{j}).retro_group1_sub.data{i} = retro_group1_sub.(properties{j});
        binned_replay_properties.(properties{j}).retro_group1_sub.mean(i) = nanmean(retro_group1_sub.(properties{j}));
        binned_replay_properties.(properties{j}).retro_group1_sub.sem(i) = nanstd(retro_group1_sub.(properties{j}))./sqrt(sum(~isnan(retro_group1_sub.(properties{j}))));
        binned_replay_properties.(properties{j}).retro_group1_sub.n(i) = sum(~isnan(retro_group1_sub.(properties{j})));

        binned_replay_properties.(properties{j}).retro_group2_sub.data{i} = retro_group2_sub.(properties{j});
        binned_replay_properties.(properties{j}).retro_group2_sub.mean(i) = nanmean(retro_group2_sub.(properties{j}));
        binned_replay_properties.(properties{j}).retro_group2_sub.sem(i) = nanstd(retro_group2_sub.(properties{j}))./sqrt(sum(~isnan(retro_group2_sub.(properties{j}))));
        binned_replay_properties.(properties{j}).retro_group2_sub.n(i) = sum(~isnan(retro_group2_sub.(properties{j})));

        binned_replay_properties.(properties{j}).prospective_group1_sub.data{i} = prospective_group1_sub.(properties{j});
        binned_replay_properties.(properties{j}).prospective_group1_sub.mean(i) = nanmean(prospective_group1_sub.(properties{j}));
        binned_replay_properties.(properties{j}).prospective_group1_sub.sem(i) = nanstd(prospective_group1_sub.(properties{j}))./sqrt(sum(~isnan(prospective_group1_sub.(properties{j}))));
        binned_replay_properties.(properties{j}).prospective_group1_sub.n(i) = sum(~isnan(prospective_group1_sub.(properties{j})));

        binned_replay_properties.(properties{j}).prospective_group2_sub.data{i} = prospective_group2_sub.(properties{j});
        binned_replay_properties.(properties{j}).prospective_group2_sub.mean(i) = nanmean(prospective_group2_sub.(properties{j}));
        binned_replay_properties.(properties{j}).prospective_group2_sub.sem(i) = nanstd(prospective_group2_sub.(properties{j}))./sqrt(sum(~isnan(prospective_group2_sub.(properties{j}))));
        binned_replay_properties.(properties{j}).prospective_group2_sub.n(i) = sum(~isnan(prospective_group2_sub.(properties{j})));
    end
    % special case:
    a = replay_sub_sub_group1.('meanAngDisplacement_pastPath');
    b = replay_sub_sub_group1.('meanAngDisplacement_futPath');
    c = a - b;
    binned_replay_properties.('angle_difference').replay_sub_sub_group1.data{i} = c;
    binned_replay_properties.('angle_difference').replay_sub_sub_group1.mean(i) = nanmean(c);
    binned_replay_properties.('angle_difference').replay_sub_sub_group1.sem(i) =  nanstd(c./sqrt(sum(~isnan(c))));

    % special case:
    a = replay_sub_sub_group2.('meanAngDisplacement_pastPath');
    b = replay_sub_sub_group2.('meanAngDisplacement_futPath');
    c = a - b;
    binned_replay_properties.('angle_difference').replay_sub_sub_group2.data{i} = c;
    binned_replay_properties.('angle_difference').replay_sub_sub_group2.mean(i) = nanmean(c);
    binned_replay_properties.('angle_difference').replay_sub_sub_group2.sem(i) =  nanstd(c./sqrt(sum(~isnan(c))));

end

% Pool ALL data
pooled_data = struct ;
aField     = fields(binned_replay_properties);
for i = 1:length(aField)
    for j = 1:length(binned_replay_properties.(aField{i}).replay_sub_sub_group1.data)
        pooled_data.(aField{i}).data{j} = [binned_replay_properties.(aField{i}).replay_sub_sub_group1.data{j};  binned_replay_properties.(aField{i}).replay_sub_sub_group2.data{j}];
    end
end


shuffles_group1_versus_group2 = struct();
num_time_bins = length(binned_replay_properties.meanAngDisplacement_futPath.replay_sub_sub_group1.data);
properties = fields(binned_replay_properties);
if strcmp(sig_test,'shuffle')
    for property = 1:length(properties)
        shuffle_diff_future_past = nan(1000,num_time_bins);
        for i = 1:1000
            for j = 1:num_time_bins
                inds1 = randperm(size(pooled_data.(properties{property}).data{j},1), size(binned_replay_properties.(properties{property}).replay_sub_sub_group1.data{j},1));
                inds2 = setdiff(1:size(pooled_data.(properties{property}).data{j},1),inds1);
                shuffle_diff_future_past(i,j) = nanmean(pooled_data.(properties{property}).data{j}(inds1,:))-nanmean(pooled_data.(properties{property}).data{j}(inds2,:));
            end
        end
        shuffles_group1_versus_group2.(properties{property}) = shuffle_diff_future_past;
    end
else
    shuffle_diff_future_past = nan(1000,num_time_bins);
    for i = 1:1000
        for j = 1:num_time_bins
            inds1 = randperm(size(pooled_data.('angle_difference').data{j},1), size(binned_replay_properties.('angle_difference').replay_sub_sub_group1.data{j},1));
            inds2 = setdiff(1:size(pooled_data.('angle_difference').data{j},1),inds1);
            shuffle_diff_future_past(i,j) = abs(nanmean(pooled_data.('angle_difference').data{j}(inds1,:))-nanmean(pooled_data.('angle_difference').data{j}(inds2,:)));
        end
    end
    shuffles_group1_versus_group2.('angle_difference') = shuffle_diff_future_past;
end

%%
% figure('Position',[1921 909 940 300])
% subplot(1,3,1)
% y_mean = binned_replay_properties.meanAngDisplacement_futPath.replay_sub_sub_group1.mean;
% y_sem = binned_replay_properties.meanAngDisplacement_futPath.replay_sub_sub_group1.sem;
% H = shadedErrorBar(bin_centers,y_mean,y_sem,'lineprops','k'); hold on;
% y_mean = binned_replay_properties.meanAngDisplacement_futPath.replay_sub_sub_group2.mean;
% y_sem = binned_replay_properties.meanAngDisplacement_futPath.replay_sub_sub_group2.sem;
% H = shadedErrorBar(bin_centers,y_mean,y_sem,'lineprops','r'); hold on;
% ylim([45 135])
% xlabel('Time since arival (s)');
% title('dist from future path')
% subplot(1,3,2)
% y_mean = binned_replay_properties.meanAngDisplacement_pastPath.replay_sub_sub_group1.mean;
% y_sem = binned_replay_properties.meanAngDisplacement_pastPath.replay_sub_sub_group1.sem;
% H = shadedErrorBar(bin_centers,y_mean,y_sem,'lineprops','k'); hold on;
% y_mean = binned_replay_properties.meanAngDisplacement_pastPath.replay_sub_sub_group2.mean;
% y_sem = binned_replay_properties.meanAngDisplacement_pastPath.replay_sub_sub_group2.sem;
% H = shadedErrorBar(bin_centers,y_mean,y_sem,'lineprops','r'); hold on;
% ylim([45 135])
% xlabel('Time since arival (s)');
% title('dist from past path')
% subplot(1,3,3) % differences
% y_mean = binned_replay_properties.angle_difference.replay_sub_sub_group1.mean;
% y_sem = binned_replay_properties.angle_difference.replay_sub_sub_group1.sem;
% H = shadedErrorBar(bin_centers,y_mean,y_sem,'lineprops','k'); hold on;
% y_mean = binned_replay_properties.angle_difference.replay_sub_sub_group2.mean;
% y_sem = binned_replay_properties.angle_difference.replay_sub_sub_group2.sem;
% H = shadedErrorBar(bin_centers,y_mean,y_sem,'lineprops','r'); hold on;
% xlim([0 10])
% xticks(0:2:10)
% xlabel('Time since arival (s)');
% title('angle_difference_all_replays','Interpreter','none')
% hline(0)
% hold on;
% ylimit = gca().YLim;
%
%
% real_difference = abs(binned_replay_properties.angle_difference.replay_sub_sub_group1.mean - binned_replay_properties.angle_difference.replay_sub_sub_group2.mean);
% sig_bins_shuffle = real_difference>quantile(shuffles_group1_versus_group2.angle_difference,0.95);
% sig_bins_means = nan(length(binned_replay_properties.angle_difference.replay_sub_sub_group1.data),1);
% sig_bins_medians = nan(length(binned_replay_properties.angle_difference.replay_sub_sub_group1.data),1);
%
% if length(inds2)>0
%     for bin = 1:length(binned_replay_properties.angle_difference.replay_sub_sub_group1.data)
%         [sig_bins_means(bin)]  = ttest2(binned_replay_properties.angle_difference.replay_sub_sub_group1.data{bin},binned_replay_properties.angle_difference.replay_sub_sub_group2.data{bin});
%         [~,sig_bins_medians(bin)]  = ranksum(binned_replay_properties.angle_difference.replay_sub_sub_group1.data{bin},binned_replay_properties.angle_difference.replay_sub_sub_group2.data{bin});
%     end
%     if strcmp(sig_test,'means')
%         sig_bins = logical(sig_bins_means);
%     elseif strcmp(sig_test,'medians')
%         sig_bins = logical(sig_bins_medians);
%     end
% end
%
% set(gcf,'PaperPositionMode','auto')
% saveas(gcf,['angle_difference_' group1 '_' group2 '_' num2str(thrForCategorization_include)],'jpeg')
% saveas(gcf,['angle_difference_' group2 '_' group2 '_' num2str(thrForCategorization_include)],'pdf')

%% Plot qualities of prospective and retrospective replays overlaid, and shuffle to see which bins differ significantly.
figure('Position', [2593 688 90 275]);
tiledlayout(3,1)
prospective_color = [.4660 0.6740 0.1880];
retrospective_color = [0.4940 0.1840 0.5560];
transparancy = 0.2;
properties = {'percent_participation_excitatory','start_to_end_distance','posteriorSpread'};
n_prospective_total = height(replay(replay.future==1&replay.time_since_real_drink_onset>0 & replay.time_since_real_drink_onset<=10,:));
n_retrospective_total = height(replay(replay.past==1&replay.time_since_real_drink_onset>0 & replay.time_since_real_drink_onset<=10,:));

if strcmp(sig_test,'shuffle')
    % shuffle for significance
    pooled_data = struct ;
    for property = 1:length(properties)
        for j = 1:length(binned_replay_properties.(properties{property}).retro_group1_sub.data) % loop across time bins
            pooled_data.(properties{property}).data{j} = [binned_replay_properties.(properties{property}).retro_group1_sub.data{j};  binned_replay_properties.(properties{property}).prospective_group1_sub.data{j}];
        end
    end
    shuffles_future_versus_past = struct();
    num_time_bins = length(binned_replay_properties.(properties{property}).retro_group1_sub.data);
    for property = 1:length(properties)
        shuffle_diff_future_past = nan(1000,17);
        for i = 1:1000
            for j = 1:num_time_bins
                inds1 = randperm(size(pooled_data.(properties{property}).data{j},1), size(binned_replay_properties.(properties{property}).retro_group1_sub.data{j},1));
                inds2 = setdiff(1:size(pooled_data.(properties{property}).data{j},1),inds1);
                shuffle_diff_future_past(i,j) = nanmean(pooled_data.(properties{property}).data{j}(inds1,:))-nanmean(pooled_data.(properties{property}).data{j}(inds2,:));
            end
        end
        shuffles_future_versus_past.(properties{property}) = shuffle_diff_future_past;
    end
end

for i = 1:length(properties)
    % plot forward data over the stopping period
    ax = nexttile(i);
    x = binned_replay_properties.(properties{i}).prospective_group1_sub.mean;
    err = binned_replay_properties.(properties{i}).prospective_group1_sub.sem;
    H=shadedErrorBar(bin_centers,x,err,'lineprops','k'); hold on
    H.patch.FaceColor = prospective_color;
    H.mainLine.Color = prospective_color;
    H.patch.FaceAlpha = transparancy;
    H.edge(1).Color = 'none'; H.edge(2).Color = 'none';
    H.mainLine.LineWidth = 1;

    % plot reverse data over the stopping period
    x2 = binned_replay_properties.(properties{i}).retro_group1_sub.mean;
    err = binned_replay_properties.(properties{i}).retro_group1_sub.sem;
    H=shadedErrorBar(bin_centers,x2,err,'lineprops','k'); hold on
    H.patch.FaceColor = retrospective_color;
    H.mainLine.Color = retrospective_color;
    H.patch.FaceAlpha = transparancy;
    H.edge(1).Color = 'none'; H.edge(2).Color = 'none';
    H.mainLine.LineWidth = 1;

    real_difference = binned_replay_properties.(properties{i}).retro_group1_sub.mean - binned_replay_properties.(properties{i}).prospective_group1_sub.mean;

    sig_bins_means = nan(length(binned_replay_properties.(properties{i}).retro_group1_sub.data),1);
    sig_bins_medians = nan(length(binned_replay_properties.(properties{i}).retro_group1_sub.data),1);
    ns = nan(length(binned_replay_properties.(properties{i}).retro_group1_sub.data),2);
    for bin = 1:length(binned_replay_properties.(properties{i}).retro_group1_sub.data)
        [sig_bins_means(bin)]  = ttest2(binned_replay_properties.(properties{i}).retro_group1_sub.data{bin},binned_replay_properties.(properties{i}).prospective_group1_sub.data{bin});
        [~,sig_bins_medians(bin)]  = ranksum(binned_replay_properties.(properties{i}).retro_group1_sub.data{bin},binned_replay_properties.(properties{i}).prospective_group1_sub.data{bin});
        ns(bin,1) = sum(~isnan(binned_replay_properties.(properties{i}).retro_group1_sub.data{bin}))
        ns(bin,2) = sum(~isnan(binned_replay_properties.(properties{i}).prospective_group1_sub.data{bin}))
    end

    n_forward = [min(ns(bin_centers<=10,2)) max(ns(bin_centers<=10,2))]
    n_reverse = [min(ns(bin_centers<=10,1)) max(ns(bin_centers<=10,1))]

    if strcmp(sig_test,'shuffle')
        sig_bins = real_difference < quantile(shuffles_future_versus_past.(properties{i}),0.025) | real_difference > quantile(shuffles_future_versus_past.(properties{i}),0.975);
    elseif strcmp(sig_test,'means')
        sig_bins = logical(sig_bins_means);
    elseif strcmp(sig_test,'medians')
        sig_bins = logical(sig_bins_medians);
    end

    if strcmp(properties{i},'percent_participation_excitatory')
        yticks([6 8 10])
        ylim([6 11])
    end
    if strcmp(properties{i},'posteriorSpread')
        yticks([2 3])
        ylim([1.5 3.5])
    end
    if strcmp(properties{i},'start_to_end_distance')
        yticks([15 20 25])
        ylim([15 25])
    end
    if i == 3
        xlabel('Time (s)')
    end
    xticks([0:2:10]);

    xtickangle(0);
    yticklabels({})
    ylimit = gca().YLim;
    plot(bin_centers(sig_bins),ylimit(2).*ones(sum(sig_bins),1),'.k')

end

set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
fig_path = '/home/caitlin/Data//Processed_Data/Manuscript_Figures';
saveas(gcf,fullfile(fig_path,'summary_for_rev_overlaid'),'jpg')
saveas(gcf,fullfile(fig_path,'summary_for_rev_overlaid'),'pdf')

%%
if ~isempty(inds2)
    % Plot angle between replay and future path, angle between replay and past path over time into stopping period.
    color_future = [.4660 0.6740 0.1880];
    color_past =   [0.4940 0.1840 0.556];
    color_shuffle_future = [0 0 0];
    color_shuffle_past = [0.3 0.3 0.3];
    figure('Position',[1921 560 800 300])
    tiledlayout(2,4)
    nexttile(1)

    colors = [.4660 0.6740 0.1880; 0.4940 0.1840 0.5560;];
    transparency_pcnt = 1;
    colors_2 = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
        [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];

    ylabels = {'|Ang. displacement|'};

    nexttile(1)
    x = binned_replay_properties.meanAngDisplacement_futPath.replay_sub_sub_group1.mean;
    err = binned_replay_properties.meanAngDisplacement_futPath.replay_sub_sub_group1.sem;
    h = shadedErrorBar(bin_centers,x,err,'lineprops','k'); hold on;
    h.patch.FaceColor = colors(1,:);
    h.patch.FaceAlpha = 0.4;
    h.mainLine.Color = colors(1,:);
    h.mainLine.LineWidth = 1;
    h.edge(1).Color = 'none'; h.edge(2).Color = 'none';

    x = binned_replay_properties.meanAngDisplacement_futPath.replay_sub_sub_group2.mean;
    err = binned_replay_properties.meanAngDisplacement_futPath.replay_sub_sub_group2.sem;
    h = shadedErrorBar(bin_centers,x,err,'lineprops','--k'); hold on;
    h.patch.FaceColor = colors_2(1,:);
    h.patch.FaceAlpha = 0.4;
    h.mainLine.Color = colors(1,:);
    h.mainLine.LineWidth = 1;
    h.edge(1).Color = 'none'; h.edge(2).Color = 'none';

    xticks(0:2:10)
    xtickangle(0)
    ylim([45 135])
    xlabel('Time since arrival (s)')
    ylabel('|Displacement|')

    nexttile(2)
    x = binned_replay_properties.meanAngDisplacement_pastPath.replay_sub_sub_group1.mean;
    err = binned_replay_properties.meanAngDisplacement_pastPath.replay_sub_sub_group1.sem;
    h = shadedErrorBar(bin_centers,x,err,'lineprops','k'); hold on;
    h.patch.FaceColor = colors(2,:);
    h.patch.FaceAlpha = 0.4;
    h.mainLine.Color = colors(2,:);
    h.mainLine.LineWidth = 1;
    h.edge(1).Color = 'none'; h.edge(2).Color = 'none';

    x = binned_replay_properties.meanAngDisplacement_pastPath.replay_sub_sub_group2.mean;
    err = binned_replay_properties.meanAngDisplacement_pastPath.replay_sub_sub_group2.sem;
    h = shadedErrorBar(bin_centers,x,err,'lineprops','--k'); hold on;
    h.patch.FaceColor = colors_2(2,:);
    h.patch.FaceAlpha = 0.4;
    h.mainLine.Color = colors(2,:);
    h.mainLine.LineWidth = 1;
    h.edge(1).Color = 'none'; h.edge(2).Color = 'none';

    xticks(0:2:10)
    xtickangle(0)
    ylim([45 135])
    xlabel('Time since arrival (s)')
    ylabel('|Displacement|')

    nexttile(5)
    x = binned_replay_properties.meanAngDisplacement_futPath.replay_sub_sub_group1.mean;
    err = binned_replay_properties.meanAngDisplacement_futPath.replay_sub_sub_group1.sem;
    h = shadedErrorBar(bin_centers,x,err,'lineprops','k'); hold on;
    h.patch.FaceColor = colors(1,:);
    h.patch.FaceAlpha = 0.4;
    h.mainLine.Color = colors(1,:);
    h.mainLine.LineWidth = 1;
    h.edge(1).Color = 'none'; h.edge(2).Color = 'none';

    x = binned_replay_properties.meanAngDisplacement_pastPath.replay_sub_sub_group1.mean;
    err = binned_replay_properties.meanAngDisplacement_pastPath.replay_sub_sub_group1.sem;
    h = shadedErrorBar(bin_centers,x,err,'lineprops','k'); hold on;
    h.patch.FaceColor = colors(2,:);
    h.patch.FaceAlpha = 0.4;
    h.mainLine.Color = colors(2,:);
    h.mainLine.LineWidth = 1;
    h.edge(1).Color = 'none'; h.edge(2).Color = 'none';

    xticks(0:2:10)
    xtickangle(0)
    ylim([45 135])
    xlabel('Time since arrival (s)')
    ylabel('|Displacement|')


    nexttile(6)
    x = binned_replay_properties.meanAngDisplacement_futPath.replay_sub_sub_group2.mean;
    err = binned_replay_properties.meanAngDisplacement_futPath.replay_sub_sub_group2.sem;
    h = shadedErrorBar(bin_centers,x,err,'lineprops','--k'); hold on;
    h.patch.FaceColor = colors(1,:);
    h.patch.FaceAlpha = 0.4;
    h.mainLine.Color = colors(1,:);
    h.mainLine.LineWidth = 1;
    h.edge(1).Color = 'none'; h.edge(2).Color = 'none';

    x = binned_replay_properties.meanAngDisplacement_pastPath.replay_sub_sub_group2.mean;
    err = binned_replay_properties.meanAngDisplacement_pastPath.replay_sub_sub_group2.sem;
    h = shadedErrorBar(bin_centers,x,err,'lineprops','--k'); hold on;
    h.patch.FaceColor = colors(2,:);
    h.patch.FaceAlpha = 0.4;
    h.mainLine.Color = colors(2,:);
    h.mainLine.LineWidth = 1;
    h.edge(1).Color = 'none'; h.edge(2).Color = 'none';

    xticks(0:2:10)
    xtickangle(0)
    ylim([45 135])
    xlabel('Time since arrival (s)')
    ylabel('|Displacement|')

    nexttile(7)
    x = binned_replay_properties.meanAngDisplacement_futPath.replay_sub_sub_group1.mean;
    err = binned_replay_properties.meanAngDisplacement_futPath.replay_sub_sub_group1.sem;
    h = shadedErrorBar(bin_centers,x,err,'lineprops','k'); hold on;
    h.patch.FaceColor = colors(1,:);
    h.patch.FaceAlpha = 0.4;
    h.mainLine.Color = colors(1,:);
    h.mainLine.LineWidth = 1;
    h.edge(1).Color = 'none'; h.edge(2).Color = 'none';

    x = binned_replay_properties.meanAngDisplacement_pastPath.replay_sub_sub_group1.mean;
    err = binned_replay_properties.meanAngDisplacement_pastPath.replay_sub_sub_group1.sem;
    h = shadedErrorBar(bin_centers,x,err,'lineprops','k'); hold on;
    h.patch.FaceColor = colors(2,:);
    h.patch.FaceAlpha = 0.4;
    h.mainLine.Color = colors(2,:);
    h.mainLine.LineWidth = 1;
    h.edge(1).Color = 'none'; h.edge(2).Color = 'none';


    x = binned_replay_properties.meanAngDisplacement_futPath.replay_sub_sub_group2.mean;
    err = binned_replay_properties.meanAngDisplacement_futPath.replay_sub_sub_group2.sem;
    h = shadedErrorBar(bin_centers,x,err,'lineprops','--k'); hold on;
    h.patch.FaceColor = colors(1,:);
    h.patch.FaceAlpha = 0.4;
    h.mainLine.Color = colors(1,:);
    h.mainLine.LineWidth = 1;
    h.edge(1).Color = 'none'; h.edge(2).Color = 'none';

    x = binned_replay_properties.meanAngDisplacement_pastPath.replay_sub_sub_group2.mean;
    err = binned_replay_properties.meanAngDisplacement_pastPath.replay_sub_sub_group2.sem;
    h = shadedErrorBar(bin_centers,x,err,'lineprops','--k'); hold on;
    h.patch.FaceColor = colors(2,:);
    h.patch.FaceAlpha = 0.4;
    h.mainLine.Color = colors(2,:);
    h.mainLine.LineWidth = 1;
    h.edge(1).Color = 'none'; h.edge(2).Color = 'none';

    xticks(0:2:10)
    xtickangle(0)
    ylim([45 135])
    xlabel('Time since arrival (s)')
    ylabel('|Displacement|')

    nexttile(3)
    % Difference of differences, real versus shuffle:
    colors = [0 0 0; 0.5 0.5 0.5];

    x = binned_replay_properties.angle_difference.replay_sub_sub_group1.mean;
    err = binned_replay_properties.angle_difference.replay_sub_sub_group1.sem;
    h = shadedErrorBar(bin_centers,x,err,'lineprops','k'); hold on;
    h.patch.FaceColor = colors(1,:);
    h.patch.FaceAlpha = 0.4;
    h.mainLine.Color = colors(1,:);
    h.mainLine.LineWidth = 1;
    h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
    hold on;
    x = binned_replay_properties.angle_difference.replay_sub_sub_group2.mean;
    err = binned_replay_properties.angle_difference.replay_sub_sub_group2.sem;
    h = shadedErrorBar(bin_centers,x,err,'lineprops','--k'); hold on;
    h.patch.FaceColor = colors(2,:);
    h.patch.FaceAlpha = 0.4;
    h.mainLine.Color = colors(1,:);
    h.mainLine.LineWidth = 1;
    h.edge(1).Color = 'none'; h.edge(2).Color = 'none';

    yline(0)

    pvals = nan(length(bin_centers),1);
    for i = 1:length(bin_centers)
        pvals(i) = ranksum(binned_replay_properties.angle_difference.replay_sub_sub_group1.data{i},binned_replay_properties.angle_difference.replay_sub_sub_group2.data{i});
    end
    hold on
    ylimit = gca().YLim(2);
    plot(bin_centers(pvals<0.05),ylimit*ones(sum(pvals<0.05)),'.k')

    xticks(0:2:10)
    xtickangle(0)
    xlabel('Time since arrival (s)')
    ylabel('Past-future')
    % Difference of differences, real versus shuffle:
    nexttile(8)
    colors = [0 0 0; 0.5 0.5 0.5];
    real_difference1 = binned_replay_properties.meanAngDisplacement_pastPath.replay_sub_sub_group1.mean;
    real_difference2 = binned_replay_properties.meanAngDisplacement_pastPath.replay_sub_sub_group2.mean;

    x = binned_replay_properties.angle_difference.replay_sub_sub_group1.mean;
    err = binned_replay_properties.angle_difference.replay_sub_sub_group1.sem;
    h = shadedErrorBar(bin_centers,x,err,'lineprops','k'); hold on;
    h.patch.FaceColor = colors(1,:);
    h.patch.FaceAlpha = 0.4;
    h.mainLine.Color = colors(1,:);
    h.mainLine.LineWidth = 1;
    h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
    hold on;
    x = binned_replay_properties.angle_difference.replay_sub_sub_group2.mean;
    err = binned_replay_properties.angle_difference.replay_sub_sub_group2.sem;
    h = shadedErrorBar(bin_centers,x,err,'lineprops','--k'); hold on;
    h.patch.FaceColor = colors(2,:);
    h.patch.FaceAlpha = 0.4;
    h.mainLine.Color = colors(1,:);
    h.mainLine.LineWidth = 1;
    h.edge(1).Color = 'none'; h.edge(2).Color = 'none';

    yline(0)

    pvals = nan(length(bin_centers),1);
    for i = 1:length(bin_centers)
        pvals(i) = ranksum(binned_replay_properties.angle_difference.replay_sub_sub_group1.data{i},binned_replay_properties.angle_difference.replay_sub_sub_group2.data{i});
    end
    hold on
    ylimit = gca().YLim(2);
    plot(bin_centers(pvals<0.05),ylimit*ones(sum(pvals<0.05)),'.k')

    xticks(0:2:10)
    xtickangle(0)
    xlabel('Time since arrival (s)')
    ylabel('Past-future')
    set(gcf, 'Color', 'white','PaperPositionMode','auto');
    fig_title = [group1 'v' group2];

    saveas(gcf,fullfile(fig_path,[fig_title '_disp' num2str(replay_dispersionThr,2) '_rats_' num2str(rats)]),'jpg')

    %% Plot qualities of replays with MEC active or inactive (Figure 4D or 4K)
    use_constant_ylims = 0;
    fig = figure();
    fig.Units = 'pixels';
    fig.Position =  [952 557 125 175];
    tiledlayout(2,3,'TileSpacing','tight')
    min_time = 0;
    max_time = inf;
    prospective_color = [.4660 0.6740 0.1880];
    retrospective_color = [0.4940 0.1840 0.5560];
    transparancy = 0.2;
    transparancy_on = 0.1;
    properties = {'start_to_end_distance','percent_participation_excitatory',};

    if strcmp(sig_test,'shuffle')
        % shuffle for significance
        pooled_data = struct ;
        for property = 1:length(properties)
            for j = 1:length(binned_replay_properties.(properties{property}).replay_sub_sub_group1.data) % loop across time bins
                pooled_data.(properties{property}).data{j} = [binned_replay_properties.(properties{property}).replay_sub_sub_group1.data{j};  binned_replay_properties.(properties{property}).replay_sub_sub_group2.data{j}];
            end
        end
        shuffles_group1_v_group2 = struct();
        num_time_bins = length(binned_replay_properties.(properties{property}).replay_sub_sub_group1.data);
        for property = 1:length(properties)
            shuffle_difprospective_group1_group2 = nan(1000,17);
            for i = 1:1000
                for j = 1:num_time_bins
                    inds1 = randperm(size(pooled_data.(properties{property}).data{j},1), size(binned_replay_properties.(properties{property}).replay_sub_sub_group1.data{j},1));
                    inds2 = setdiff(1:size(pooled_data.(properties{property}).data{j},1),inds1);
                    shuffle_difprospective_group1_group2(i,j) = abs(nanmean(pooled_data.(properties{property}).data{j}(inds1,:))-nanmean(pooled_data.(properties{property}).data{j}(inds2,:)));
                end
            end
            shuffles_group1_v_group2.(properties{property}) = shuffle_difprospective_group1_group2;
        end
    end

    for i = 1:length(properties)
        % plot group1 data over the stopping period
        nexttile(3*(i-1)+1,[1,2]);
        x = binned_replay_properties.(properties{i}).replay_sub_sub_group1.mean;
        err = binned_replay_properties.(properties{i}).replay_sub_sub_group1.sem;
        H=shadedErrorBar(bin_centers,x,err,'lineprops','k'); hold on
        H.patch.FaceAlpha = transparancy;
        H.edge(1).Color = 'none'; H.edge(2).Color = 'none';
        H.mainLine.LineWidth = 1;

        % plot group2 data over the stopping period
        x2 = binned_replay_properties.(properties{i}).replay_sub_sub_group2.mean;
        err = binned_replay_properties.(properties{i}).replay_sub_sub_group2.sem;
        H=shadedErrorBar(bin_centers,x2,err,'lineprops','--k'); hold on
        H.patch.FaceAlpha = transparancy;
        H.edge(1).Color = 'none'; H.edge(2).Color = 'none';
        H.mainLine.LineWidth = 1;

        group_1_n = [min(binned_replay_properties.(properties{i}).replay_sub_sub_group1.n(bin_centers<=10)) max(binned_replay_properties.(properties{i}).replay_sub_sub_group1.n(bin_centers<=10))]
        group_2_n = [min(binned_replay_properties.(properties{i}).replay_sub_sub_group2.n(bin_centers<=10)) max(binned_replay_properties.(properties{i}).replay_sub_sub_group2.n(bin_centers<=10))]

        if strcmp(properties{i},'start_to_end_distance')
            ylim([15 25]);
        end
        yticklabels({})
        real_difference = abs(binned_replay_properties.(properties{i}).replay_sub_sub_group1.mean - binned_replay_properties.(properties{i}).replay_sub_sub_group2.mean);

        sig_bins_means = nan(length(binned_replay_properties.(properties{i}).replay_sub_sub_group1.data),1);
        sig_bins_medians = nan(length(binned_replay_properties.(properties{i}).replay_sub_sub_group1.data),1) ;
        for bin = 1:length(binned_replay_properties.(properties{i}).replay_sub_sub_group1.data)
            [sig_bins_means(bin)]  = ttest2(binned_replay_properties.(properties{i}).replay_sub_sub_group1.data{bin},binned_replay_properties.(properties{i}).replay_sub_sub_group2.data{bin});
            [~,sig_bins_medians(bin)]  = ranksum(binned_replay_properties.(properties{i}).replay_sub_sub_group1.data{bin},binned_replay_properties.(properties{i}).replay_sub_sub_group2.data{bin});
        end

        if strcmp(sig_test,'shuffle')
            sig_bins = real_difference > quantile(shuffles_group1_v_group2.(properties{i}),0.975);
        elseif strcmp(sig_test,'means')
            sig_bins = logical(sig_bins_means);
        elseif strcmp(sig_test,'medians')
            sig_bins = logical(sig_bins_medians);
        end

        ylimit = gca().YLim;
        plot(bin_centers(sig_bins),ylimit(2).*ones(sum(sig_bins),1),'.k')
        if i~=1
            xlabel('Time (s)')
        end
        xticks([0:5:10]);
        xtickangle(0);

        nexttile(3*(i-1)+3);
        x1 = replay_sub_group1.(properties{i})(replay_sub_group1.time_since_real_drink_onset > min_time & replay_sub_group1.time_since_real_drink_onset <= max_time);
        group1_mean = nanmean(x1);
        group1_sem = nanstd(x1)./sqrt(sum(~isnan(x1)));
        x2 = replay_sub_group2.(properties{i})(replay_sub_group2.time_since_real_drink_onset > min_time & replay_sub_group2.time_since_real_drink_onset <= max_time);
        group2_mean = nanmean(x2);
        group2_sem = nanstd(x2)./sqrt(sum(~isnan(x2)));

        x = [1 2];
        data = [group1_mean group2_mean];
        sem = [group1_sem group2_sem];

        x_boxplot = [ones(size(x1)); 2*ones(size(x2))];
        y_boxplot = [x1; x2];

        b = bar(x,data,'FaceColor','flat'); hold on
        b(1).CData(1,:) = [0 0 0];
        b(1).CData(2,:) = [0.5 0.5 0.5]
        er = errorbar(x,data,sem,sem);
        er.Color = [0 0 0];
        er.LineStyle = 'none';
        b.EdgeColor = 'none';
        %title(num2str(pval))
        xticklabels({''})

        box off
        hold on

        ax = gca;
        ax.FontSize = 8;
        set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
        fig_path = '/home/caitlin/Data//Processed_Data/Manuscript_Figures';
        saveas(fig,fullfile(fig_path,['summary_' group1 '_' group2 '_overlaid']),'jpg')
        saveas(fig,fullfile(fig_path,['summary_' group1 '_' group2 '_overlaid']),'pdf')
    end

    % stats for fig 4D or 4K
    num_shuffles = 10000;
    use_means_for_shuffle = 1;
    time_window_1 = [0 3];
    time_window_2 = [3 10];
    property = 'percent_participation_excitatory';
    %property = 'start_to_end_distance';

    real_off = replay_sub.(property)(replay_sub.time_since_real_drink_onset >= 0 & ...
        replay_sub.time_since_real_drink_onset < 10 & replay_sub.laser_state_binary==0);
    real_on = replay_sub.(property)(replay_sub.time_since_real_drink_onset >= 0 & ...
        replay_sub.time_since_real_drink_onset < 10 & replay_sub.laser_state_binary==1);

    real_early_off = replay_sub.(property)(replay_sub.time_since_real_drink_onset >= time_window_1(1) & ...
        replay_sub.time_since_real_drink_onset < time_window_1(2) & replay_sub.laser_state_binary==0);
    real_late_off = replay_sub.(property)(replay_sub.time_since_real_drink_onset >= time_window_2(1) & ...
        replay_sub.time_since_real_drink_onset < time_window_2(2) & replay_sub.laser_state_binary==0);

    real_early_on = replay_sub.(property)(replay_sub.time_since_real_drink_onset >= time_window_1(1) & ...
        replay_sub.time_since_real_drink_onset < time_window_1(2) & replay_sub.laser_state_binary==1);
    real_late_on = replay_sub.(property)(replay_sub.time_since_real_drink_onset >= time_window_2(1) & ...
        replay_sub.time_since_real_drink_onset < time_window_2(2) & replay_sub.laser_state_binary==1);

    if use_means_for_shuffle==1
        real_late_off_minus_real_late_on = mean(real_late_off) - mean(real_late_on);
        real_early_off_minus_real_early_on = mean(real_early_off) - mean(real_early_on);
        real_interaction = real_late_off_minus_real_late_on - real_early_off_minus_real_early_on;

        early = [real_early_on; real_early_off];
        late = [real_late_on; real_late_off];
        shuffled_early_diffs = nan(num_shuffles,1);
        shuffled_late_diffs = nan(num_shuffles,1);
        for n = 1:num_shuffles
            shuffled_early_off_inds = randperm(length(early),length(real_early_off));
            shuffled_early_on_inds = setdiff(1:length(early),shuffled_early_off_inds);
            shuffled_early_diffs(n) = mean(early(shuffled_early_off_inds))-mean(early(shuffled_early_on_inds));

            shuffled_late_off_inds = randperm(length(late),length(real_late_off));
            shuffled_late_on_inds = setdiff(1:length(late),shuffled_late_off_inds);
            shuffled_late_diffs(n) = mean(late(shuffled_late_off_inds))-mean(late(shuffled_late_on_inds));
        end

        shuffled_interactions = shuffled_late_diffs-shuffled_early_diffs;
        sig_interaction = (1+sum(shuffled_interactions>abs(real_interaction))+sum(shuffled_interactions<(-1*abs(real_interaction))))/(num_shuffles+1);
        sig_late = (1+sum(shuffled_late_diffs > abs(real_late_off_minus_real_late_on))+sum(shuffled_late_diffs < (-1*abs(real_late_off_minus_real_late_on))))/(num_shuffles+1);
        sig_early = (1+sum(shuffled_early_diffs > abs(real_early_off_minus_real_early_on)) +sum(shuffled_early_diffs < (-1*abs(real_early_off_minus_real_early_on))))/(num_shuffles+1);
    else
        real_late_off_minus_real_late_on = median(real_late_off) - median(real_late_on);
        real_early_off_minus_real_early_on = median(real_early_off) - median(real_early_on);
        real_interaction = real_late_off_minus_real_late_on - real_early_off_minus_real_early_on;

        num_shuffles = 10000;
        early = [real_early_on; real_early_off];
        late = [real_late_on; real_late_off];
        shuffled_early_diffs = nan(num_shuffles,1);
        shuffled_late_diffs = nan(num_shuffles,1);
        for n = 1:num_shuffles
            shuffled_early_off_inds = randperm(length(early),length(real_early_off));
            shuffled_early_on_inds = setdiff(1:length(early),shuffled_early_off_inds);
            shuffled_early_diffs(n) = median(early(shuffled_early_off_inds))-median(early(shuffled_early_on_inds));

            shuffled_late_off_inds = randperm(length(late),length(real_late_off));
            shuffled_late_on_inds = setdiff(1:length(late),shuffled_late_off_inds);
            shuffled_late_diffs(n) = median(late(shuffled_late_off_inds))-median(late(shuffled_late_on_inds));
        end

        shuffled_interactions = shuffled_late_diffs-shuffled_early_diffs;
        sig_interaction = (1+sum(shuffled_interactions>abs(real_interaction))+sum(shuffled_interactions<(-1*abs(real_interaction))))/(num_shuffles+1);
        sig_late = (1+sum(shuffled_late_diffs > abs(real_late_off_minus_real_late_on))+sum(shuffled_late_diffs < (-1*abs(real_late_off_minus_real_late_on))))/(num_shuffles+1);
        sig_early = (1+sum(shuffled_early_diffs > abs(real_early_off_minus_real_early_on)) +sum(shuffled_early_diffs < (-1*abs(real_early_off_minus_real_early_on))))/(num_shuffles+1);
    end

    % compare to two way anova
    % data_table_sub = table();
    % data_table_sub.property = [real_early_off; real_late_off; real_early_on; real_late_on];
    % data_table_sub.group = [zeros(height(real_early_off),1); zeros(height(real_late_off),1); ones(height(real_early_on),1); ones(height(real_late_on),1)];
    % data_table_sub.time = [zeros(height(real_early_off),1); ones(height(real_late_off),1); zeros(height(real_early_on),1); ones(height(real_late_on),1)];
    % groupOrder = ["1","2"];
    % namedGroup = categorical(data_table_sub.group,0:1,groupOrder);
    % directionOrder = ["Early","Late"];
    % namedDirection = categorical(data_table_sub.time,[0,1],directionOrder);
    % data_sub = data_table_sub.property;
    % [p,tbl,stats] = anovan(data_sub',{namedGroup,namedDirection},'model','interaction','varnames',{'namedGroup','namedDirection'});
    %  [results,~,~,gnames] = multcompare(stats,"Dimension",[1 2])

    [p,h,z] = ranksum(real_early_off,real_early_on)
    [p,h,z] = ranksum(real_late_off,real_late_on)
    data_early_off = nanmean(real_early_off);
    sem_early_off = nanstd(real_early_off)./sqrt(sum(~isnan(real_early_off)));
    data_early_on = nanmean(real_early_on);
    sem_early_on = nanstd(real_early_on)./sqrt(sum(~isnan(real_early_on)));
    data_late_off = nanmean(real_late_off);
    sem_late_off = nanstd(real_late_off)./sqrt(sum(~isnan(real_late_off)));
    data_late_on = nanmean(real_late_on);
    sem_late_on = nanstd(real_late_on)./sqrt(sum(~isnan(real_late_on)));

    data = [data_early_off, data_early_on, data_late_off, data_late_on];
    err = [sem_early_off, sem_early_on, sem_late_off, sem_late_on];

    colors = [.4660 0.6740 0.1880; 0.4940 0.1840 0.5560;];
    transparency_pcnt = 0.5;
    colors_2 = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
        [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];
    b1 = bar(1,data(1)); hold on;
    e1 = errorbar(1,data(1),err(1),'k','linestyle','none');
    e1.CapSize = 4;
    hold on
    b2 = bar(2,data(2)); hold on;
    e2 = errorbar(2,data(2),err(2),'k','linestyle','none');
    e2.CapSize = 4;
    b3 = bar(4,data(3)); hold on;
    e3 = errorbar(4,data(3),err(3),'k','linestyle','none');
    e3.CapSize = 4;
    b4 = bar(5,data(4)); hold on;
    e4 = errorbar(5,data(4),err(4),'k','linestyle','none');
    e4.CapSize = 4;
    b1.FaceColor = colors(1,:);
    b1.EdgeColor = 'none';
    b2.FaceColor = colors_2(1,:);
    b2.EdgeColor = 'none';
    b3.FaceColor = colors(2,:);
    b3.EdgeColor = 'none';
    b4.FaceColor = colors_2(2,:);
    b4.EdgeColor = 'none';

    if strcmp(property,'dispersion')
        ylim([5 9])
    end
    if strcmp(property,'percent_participation_excitatory')
        ylim([10 16])
    end

    early_off_n = length(real_early_off);
    early_on_n = length(real_early_on);
    late_off_n = length(real_late_off);
    late_on_n = length(real_late_on);
end
