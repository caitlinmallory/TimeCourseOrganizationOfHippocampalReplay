%%

grouping = 'rat';
grouping = 'session';

if strcmp(grouping,'rat')
    data_tbl.unique_n = data_tbl.rat_label;
elseif strcmp(grouping,'session')
    data_tbl.unique_n = data_tbl.unique_session_id;
end
unique_n = unique(data_tbl.unique_n);

fig1 = figure();
fig1.Units = 'centimeters';
fig1.Position = [50.8000 26.2467 8 6.5];

tiledlayout(2,2)
nexttile()

group_summary = struct();

for n = 1:length(unique_n)
replay = data_tbl(data_tbl.unique_n == unique_n(n) & data_tbl.dispersion > replay_dispersionThr & data_tbl.duration >= replay_durationThr & ...
    data_tbl.drink_period_time >= min_stopping_period_duration & data_tbl.drink_period_time <= max_stopping_period_duration,:);

% convert from radians to degrees
replay.mean_abs_angular_displacement_ratLoc = rad2deg(replay.mean_abs_angular_displacement_ratLoc);
replay.meanAngDisplacement_futPath = rad2deg(replay.meanAngDisplacement_futPath);
replay.meanAngDisplacement_pastPath = rad2deg(replay.meanAngDisplacement_pastPath);
replay.angle_between_past_future_trajectory = rad2deg(replay.angle_between_past_future_trajectory);
replay.angDisplacement_futPath = cellfun(@(x) rad2deg(x),replay.angDisplacement_futPath, 'UniformOutput', false);
replay.angDisplacement_pastPath = cellfun(@(x) rad2deg(x),replay.angDisplacement_pastPath, 'UniformOutput', false);
replay.replay_shuffle_past = cellfun(@(x) rad2deg(x),replay.replay_shuffle_past, 'UniformOutput', false);
replay.replay_shuffle_future = cellfun(@(x) rad2deg(x),replay.replay_shuffle_future, 'UniformOutput', false);
replay.all_angularDisplacement_past_shuffled = cellfun(@(x) rad2deg(x),replay.all_angularDisplacement_past_shuffled, 'UniformOutput', false);
replay.all_angularDisplacement_future_shuffled = cellfun(@(x) rad2deg(x),replay.all_angularDisplacement_future_shuffled, 'UniformOutput', false);

% remove shuffles that came BEFORE or AFTER the current trial
modifyCell = @(shuffle,num_trials) shuffle(num_trials>= num_trials_away_low & num_trials <= num_trials_away_high);
replay.replay_shuffle_past = cellfun(modifyCell, replay.replay_shuffle_past, replay.num_trials_between_current_stopping_period_and_shuffle, 'UniformOutput', false);
replay.replay_shuffle_future = cellfun(modifyCell, replay.replay_shuffle_future, replay.num_trials_between_current_stopping_period_and_shuffle, 'UniformOutput', false);

replay(replay.angle_between_past_future_trajectory < angle_between_past_future_trajectory_Thr(1) | ...
    replay.angle_between_past_future_trajectory > angle_between_past_future_trajectory_Thr(2),:) = [];


if home_trials_only == 1
    replay(replay.home_event~=1,:) = [];
end
if away_trials_only == 1
    replay(replay.away_event~=1,:) = [];
end
if align_to_start_of_anticipatory_licking == 1
    remove_trials_with_anticipatory_licking = 0;
    replay.time_since_real_drink_onset(replay.time_between_anticipatory_licking_and_drink_start > 0.5) = replay.time_since_real_drink_onset(replay.time_between_anticipatory_licking_and_drink_start > 0.5) + replay.time_between_anticipatory_licking_and_drink_start(replay.time_between_anticipatory_licking_and_drink_start > 0.5);
end
if remove_trials_with_anticipatory_licking==1
    replay(replay.time_between_anticipatory_licking_and_drink_start > 0.5,:) = [];
end

if align_to_drink_onset == 1
replay.time_into_stopping_period = replay.time_since_real_drink_onset;
elseif align_to_drink_offset == 1
replay.time_into_stopping_period = replay.time_till_real_drink_offset;
end


% replay.time_into_stopping_period(replay.rat_label==6|replay.rat_label==4) = replay.time_since_real_drink_onset(replay.rat_label==6|replay.rat_label==4);
% replay.time_into_stopping_period(replay.rat_label~=6|replay.rat_label~=4) = replay.time_since_stopping_period_onset(replay.rat_label~=6|replay.rat_label~=4); 

% mean_shuffle_future = nanmean(replay.replay_shuffle_future);
% mean_shuffle_past = nanmean(replay.replay_shuffle_past);
if plot_shuffles == 1
    replay.mean_shuffle_future = cellfun(@nanmean,replay.replay_shuffle_future);
    replay.mean_shuffle_past = cellfun(@nanmean,replay.replay_shuffle_past);
    replay.replay_all = cellfun(@(x,y) horzcat(x,y), replay.replay_shuffle_future, replay.replay_shuffle_past,'UniformOutput',false);
    replay.mean_shuffle_all = cellfun(@nanmean,replay.replay_all);
%     replay.mean_shuffle_future = replay.replay_shuffle_future;
%     replay.mean_shuffle_past = replay.replay_shuffle_past;

    properties = {'meanAngDisplacement_futPath'; 'meanAngDisplacement_pastPath';
        'startDistFromRat'; 'endDistFromRat'; 'mean_abs_angular_displacement_ratLoc';'distance';'dispersion';'duration';'mean_shuffle_past';'mean_shuffle_future';'mean_shuffle_all'};
else
    properties = {'meanAngDisplacement_futPath'; 'meanAngDisplacement_pastPath';
        'startDistFromRat'; 'endDistFromRat'; 'mean_abs_angular_displacement_ratLoc';'distance';'dispersion';'duration'};
end

windowSize = 1;
windowShift = 0.5;
start_time = 0;
end_time = 10;
bin_start = (start_time:windowShift:(end_time-windowSize))';
bin_end = bin_start + windowSize;

bin_edges = [bin_start bin_end];
bin_centers = mean(bin_edges,2);

binned_properties = table();
for property = 1:length(properties)
    binned_data = table();
    for time_bin = 1:length(bin_edges)
        binned_data.data{time_bin} = replay.(properties{property})(replay.time_into_stopping_period>= bin_edges(time_bin,1) & replay.time_into_stopping_period<bin_edges(time_bin,2));
        binned_data.mean(time_bin) = nanmean(binned_data.data{time_bin});
        binned_data.sem(time_bin) = nanstd(binned_data.data{time_bin})/sqrt(sum(~isnan(binned_data.data{time_bin})));
        binned_data.low(time_bin) = quantile(binned_data.data{time_bin},0.025);
        binned_data.high(time_bin) = quantile(binned_data.data{time_bin},0.975);
    end
    binned_properties.(properties{property}) = binned_data;
    group_summary(n).(properties{property}) = binned_data;
end
binned_data = table();
for time_bin = 1:length(bin_centers)
    a = group_summary(n).meanAngDisplacement_futPath(time_bin,:).data{:};
    b = group_summary(n).meanAngDisplacement_pastPath(time_bin,:).data{:};
    c = b-a;
    binned_data.data{time_bin} = c;
    binned_data.mean(time_bin) = nanmean(c);
    binned_data.sem(time_bin) = nanstd(c)/sqrt(sum(~isnan(c)));
    binned_data.low(time_bin) = quantile(c,0.025);
    binned_data.high(time_bin) = quantile(c,0.975);
end
group_summary(n).past_minus_future_angle = binned_data;



graph_x_low = 0;
graph_x_high = 10;

if align_to_drink_offset == 1
    bin_centers = -1*(bin_centers);
    graph_x_low = -10;
    graph_x_high = 0;
end
% Plot angle between replay and future path, angle between replay and past path over time into stopping period.
properties_to_plot = {'meanAngDisplacement_futPath' 'meanAngDisplacement_pastPath'};
ylabels = {};
fig_title = 'replay_angle_past_future_over_timed';
%color_future = [ 0.3922    0.8314    0.0745];
color_future = [ 62 150 81]./255;
color_past =   [0.6392    0.0118    0.6392];
color_shuffle_future = [0 0 0];
color_shuffle_past = [0.3 0.3 0.3];

fig1()

nexttile(1)
plot(bin_centers,binned_properties.(properties_to_plot{1}).mean); hold on;
if align_to_drink_onset==1
    xlabel('Time since stopping (s)')
elseif align_to_drink_offset == 1
    xlabel('Time till run (s)')
    rho1 = -1*rho1; rho2 = -1*rho2;
end
ylabel(ylabels);

plot([graph_x_low graph_x_high],[(90) (90)],'--','color',[0.5 0.5 0.5],'LineWidth',2); hold on;
xticks(graph_x_low:2:graph_x_high)
xtickangle(0)
if strcmp(grouping,'rat')
    ylim_deg = [45 135];
elseif strcmp(grouping,'session')
    ylim_deg = [0 180];
end
ylim ([(ylim_deg(1)) (ylim_deg(2))]);
yticks([(ylim_deg(1)) (90) (ylim_deg(2))])
yticklabels({num2str(ylim_deg(1)),num2str(90),num2str(ylim_deg(2))})
ylimit = gca().YLim; y_range = ylimit(2)-ylimit(1);
hold on;
box off

nexttile(2)
plot(bin_centers,binned_properties.(properties_to_plot{2}).mean); hold on;
if align_to_drink_onset==1
    xlabel('Time since stopping (s)')
elseif align_to_drink_offset == 1
    xlabel('Time till run (s)')
    rho1 = -1*rho1; rho2 = -1*rho2;
end
ylabel(ylabels);

plot([graph_x_low graph_x_high],[(90) (90)],'--','color',[0.5 0.5 0.5],'LineWidth',2); hold on;
xticks(graph_x_low:2:graph_x_high)
xtickangle(0)
if strcmp(grouping,'rat')
    ylim_deg = [45 135];
elseif strcmp(grouping,'session')
    ylim_deg = [0 180];
end
ylim ([(ylim_deg(1)) (ylim_deg(2))]);
yticks([(ylim_deg(1)) (90) (ylim_deg(2))])
yticklabels({num2str(ylim_deg(1)),num2str(90),num2str(ylim_deg(2))})
ylimit = gca().YLim; y_range = ylimit(2)-ylimit(1);
hold on;
box off


nexttile(3)
plot(bin_centers,binned_properties.(properties_to_plot{2}).mean - binned_properties.(properties_to_plot{1}).mean); hold on;
if align_to_drink_onset==1
    xlabel('Time since stopping (s)')
elseif align_to_drink_offset == 1
    xlabel('Time till run (s)')
    rho1 = -1*rho1; rho2 = -1*rho2;
end
ylabel(ylabels);
if strcmp(grouping,'rat')
    ylim_deg = [-30 60];
elseif strcmp(grouping,'session')
    ylim_deg = [-90 90];
end
ylim ([(ylim_deg(1)) (ylim_deg(2))]);
yticks([(ylim_deg(1)) 0 (ylim_deg(2))])
yticklabels({num2str(ylim_deg(1)),num2str(0),num2str(ylim_deg(2))})

plot([graph_x_low graph_x_high],[(0) (0)],'--','color',[0.5 0.5 0.5],'LineWidth',2); hold on;
xticks(graph_x_low:2:graph_x_high)
xtickangle(0)
box off
end

%     text(0.25,ylimit(1)+1.15*y_range,['R=' num2str(rho1,2) ' P=' num2str(p1,2)],'color',color_future,'FontSize',8)
%     text(0.25,ylimit(1)+1.08*y_range,['R=' num2str(rho2,2) ' P=' num2str(p2,2)],'color',color_past,'FontSize',8)
    set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
    export_fig(fullfile(fig_path,[fig_title '_disp' num2str(replay_dispersionThr,2) '_rats_' num2str(rats)]),'-jpeg')
    saveas(gcf,fullfile(fig_path,[fig_title '_disp' num2str(replay_dispersionThr,2) '_rats_' num2str(rats)]),'pdf')   

    %%
    Variables = {'meanAngDisplacement_futPath','meanAngDisplacement_pastPath'};
    group = struct();

    for n = 1:length(unique_n)
        for i = 1:length(Variables)
            group.(Variables{i})(n,:) = group_summary(n).(Variables{i}).mean;
        end
    end

    nexttile(1)
    plot(bin_centers,nanmean(group.(Variables{1})),'k','LineWidth',2)
    nexttile(2)
    plot(bin_centers,nanmean(group.(Variables{2})),'k','LineWidth',2)
    nexttile(3)
    plot(bin_centers,nanmean(group.(Variables{2}) - group.(Variables{1})),'k','LineWidth',2)
    
    data = group.(Variables{2}) - group.(Variables{1});
    pval = ttest(data,0);
     plot(bin_centers(pval==1),repmat((50),[sum(pval==1),1]),'.k') 

    set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
    saveas(gcf,fullfile(fig_path,[fig_title '_disp' num2str(replay_dispersionThr,2) grouping ]),'pdf')   
   saveas(gcf,fullfile(fig_path,[fig_title '_disp' num2str(replay_dispersionThr,2) grouping ]),'jpeg')   
    %%
    fig2 = figure();
    fig2.Units = 'centimeters';
    fig2.Position = [50.8000 26.2467 8 6.5];
    tiledlayout(2,2)
    nexttile(1)
    shadedErrorBar(bin_centers,nanmean(group.(Variables{1})),nanstd(group.(Variables{1}))./sqrt(sum(~isnan(group.(Variables{1})))));
    hold on;
    plot([graph_x_low graph_x_high],[(90) (90)],'--','color',[0.5 0.5 0.5],'LineWidth',2); hold on;
    xticks(graph_x_low:2:graph_x_high)
    xtickangle(0)
    ylim_deg = [45 135];
    ylim ([(ylim_deg(1)) (ylim_deg(2))]);
    yticks([(ylim_deg(1)) (90) (ylim_deg(2))])
    yticklabels({num2str(ylim_deg(1)),num2str(90),num2str(ylim_deg(2))})
    ylimit = gca().YLim; y_range = ylimit(2)-ylimit(1);
    xlabel('Time since arrival')
    hold on;
    box off


    nexttile(2)
    shadedErrorBar(bin_centers,nanmean(group.(Variables{2})),nanstd(group.(Variables{2}))./sqrt(sum(~isnan(group.(Variables{2})))))
    hold on;
    plot([graph_x_low graph_x_high],[(90) (90)],'--','color',[0.5 0.5 0.5],'LineWidth',2); hold on;
    xticks(graph_x_low:2:graph_x_high)
    xtickangle(0)
    ylim_deg = [45 135];
    ylim ([(ylim_deg(1)) (ylim_deg(2))]);
    yticks([(ylim_deg(1)) (90) (ylim_deg(2))])
    yticklabels({num2str(ylim_deg(1)),num2str(90),num2str(ylim_deg(2))})
    ylimit = gca().YLim; y_range = ylimit(2)-ylimit(1);
    xlabel('Time since arrival')
    hold on;
    box off



    nexttile(3)
    data = group.(Variables{2}) - group.(Variables{1});
    pval = ttest(data,0);
    shadedErrorBar(bin_centers,nanmean(data),nanstd(data)./sqrt(sum(~isnan(data))))
    hold on;
    plot([graph_x_low graph_x_high],[(0) (0)],'--','color',[0.5 0.5 0.5],'LineWidth',2); hold on;
    xticks(graph_x_low:2:graph_x_high)
    xtickangle(0)
    ylim_deg = [-10 40];
    ylim ([(ylim_deg(1)) (ylim_deg(2))]);
    yticks([(ylim_deg(1)) (0) (ylim_deg(2))])
    yticklabels({num2str(ylim_deg(1)),num2str(0),num2str(ylim_deg(2))})
    ylimit = gca().YLim; y_range = ylimit(2)-ylimit(1);
    hold on;
    box off
    hold on;
    xlabel('Time since arrival')
    plot(bin_centers(pval==1),repmat((ylim_deg(2)),[sum(pval==1),1]),'.k')


    set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
    saveas(gcf,fullfile(fig_path,[fig_title '_disp' num2str(replay_dispersionThr,2) grouping '_sem' ]),'pdf')  
    saveas(gcf,fullfile(fig_path,[fig_title '_disp' num2str(replay_dispersionThr,2) grouping '_sem' ]),'jpeg')  
%%

