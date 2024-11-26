
%%

replay = replay(~isnan(replay.meanAngDisplacement_futPath) & ~isnan(replay.meanAngDisplacement_pastPath),:);

grouping = 'rat';
%grouping = 'session';
if strcmp(grouping,'rat')
    ylim_deg_difference_plot = [-45 45];
elseif strcmp(grouping,'session')
    ylim_deg_difference_plot = [-90 90];
end
ylim (ylim_deg_difference_plot);
if strcmp(grouping,'rat')
    replay.unique_n = replay.rat_label;
elseif strcmp(grouping,'session')
    replay.unique_n = replay.unique_session_id;
end
unique_n = unique(replay.unique_n);
time_of_crossing = nan(length(unique_n),1);
all_coefficients = nan(length(unique_n),3);

windowSize = 2;
windowShift = 0.5;
start_time = 0;
end_time = 14;
bin_start = (start_time:windowShift:(end_time-windowSize))';
bin_end = bin_start + windowSize;
bin_edges = [bin_start bin_end];
bin_centers = mean(bin_edges,2);
properties_to_plot = {'angDisplacement_futPath' 'angDisplacement_pastPath'};
%% First do a grand binning of all data:

binned_properties = table();
for property = 1:length(properties_to_plot)
    binned_data = table();
    for time_bin = 1:length(bin_edges)
        binned_data.data{time_bin} = abs(cell2mat(replay.(properties_to_plot{property})(replay.time_into_stopping_period>= bin_edges(time_bin,1) & replay.time_into_stopping_period<bin_edges(time_bin,2))));
        binned_data.mean(time_bin) = nanmean(binned_data.data{time_bin});
        binned_data.median(time_bin) = nanmedian(binned_data.data{time_bin});
        binned_data.sem(time_bin) = nanstd(binned_data.data{time_bin})/sqrt(sum(~isnan(binned_data.data{time_bin})));
        binned_data.n(time_bin) = sum(~isnan(binned_data.data{time_bin}));
        binned_data.std(time_bin) = nanstd(binned_data.data{time_bin});
        binned_data.low(time_bin) = quantile(binned_data.data{time_bin},0.025);
        binned_data.high(time_bin) = quantile(binned_data.data{time_bin},0.975);
        binned_data.sum(time_bin) = nansum(binned_data.data{time_bin});
        binned_data.pcnt(time_bin) = nansum(binned_data.data{time_bin})./length(binned_data.data{time_bin});
    end
    binned_properties.(properties_to_plot{property}) = binned_data;
end

binned_data = table();
for time_bin = 1:length(bin_centers)
    a = binned_properties.angDisplacement_futPath(time_bin,:).data{:};
    b = binned_properties.angDisplacement_pastPath(time_bin,:).data{:};
    c = a-b;
    binned_data.data{time_bin} = c;
    binned_data.mean(time_bin) = nanmean(c);
    binned_data.n(time_bin) = sum(~isnan(c));
    binned_data.sem(time_bin) = nanstd(c)/sqrt(sum(~isnan(c)));
    binned_data.std(time_bin) = nanstd(c);
    binned_data.low(time_bin) = quantile(c,0.025);
    binned_data.high(time_bin) = quantile(c,0.975);
end
binned_properties.future_minus_past_angle = binned_data;

%%

fig1 = figure();
fig1.Position = [600 600 300 250];
tiledlayout(2,2)
nexttile()

group_summary = struct();

for n = 1:length(unique_n)
    replay_sub = replay(replay.unique_n == unique_n(n),:);


    binned_properties_sub = table();
    for property = 1:length(properties_to_plot)
        binned_data = table();
        for time_bin = 1:length(bin_edges)
            binned_data.data{time_bin} = abs(cell2mat(replay_sub.(properties_to_plot{property})(replay_sub.time_into_stopping_period>= bin_edges(time_bin,1) & replay_sub.time_into_stopping_period<bin_edges(time_bin,2))));
            binned_data.mean(time_bin) = nanmean(binned_data.data{time_bin});
            binned_data.sem(time_bin) = nanstd(binned_data.data{time_bin})/sqrt(sum(~isnan(binned_data.data{time_bin})));
            binned_data.low(time_bin) = quantile(binned_data.data{time_bin},0.025);
            binned_data.high(time_bin) = quantile(binned_data.data{time_bin},0.975);
        end
        binned_properties_sub.(properties_to_plot{property}) = binned_data;
        group_summary(n).(properties_to_plot{property}) = binned_data;
    end
    binned_data = table();
    for time_bin = 1:length(bin_edges)
        a = group_summary(n).angDisplacement_futPath(time_bin,:).data{:};
        b = group_summary(n).angDisplacement_pastPath(time_bin,:).data{:};
        c = a-b;
        binned_data.data{time_bin} = c;
        binned_data.mean(time_bin) = nanmean(c);
        binned_data.sem(time_bin) = nanstd(c)/sqrt(sum(~isnan(c)));
        binned_data.low(time_bin) = quantile(c,0.025);
        binned_data.high(time_bin) = quantile(c,0.975);
    end
    binned_properties_sub.future_minus_past_angle = binned_data;
    group_summary(n).future_minus_past_angle = binned_data;


    graph_x_low = 0;
    graph_x_high = 10;

    if align_to_drink_offset == 1
        bin_centers = -1*(bin_centers);
        graph_x_low = -10;
        graph_x_high = 0;
    end
    % Plot angle between replay and future path, angle between replay and past path over time into stopping period.
    properties_to_plot = {'angDisplacement_futPath' 'angDisplacement_pastPath'};
    ylabels = {};
    fig_title = 'replay_angle_past_future_over_timed';
    %color_future = [ 0.3922    0.8314    0.0745];
    color_future = [ 62 150 81]./255;
    color_past =   [0.6392    0.0118    0.6392];
    color_shuffle_future = [0 0 0];
    color_shuffle_past = [0.3 0.3 0.3];

    fig1()
    nexttile(1)
    plot(bin_centers,binned_properties_sub.(properties_to_plot{1}).mean); hold on;
    % if align_to_drink_onset==1
    %     xlabel('Time since arrival (s)')
    % elseif align_to_drink_offset == 1
    %     xlabel('Time till run (s)')
    %     rho1 = -1*rho1; rho2 = -1*rho2;
    % end
    ylabel(ylabels);

    plot([graph_x_low graph_x_high],[(90) (90)],'--','color',[0.5 0.5 0.5],'LineWidth',2); hold on;
    xticks(graph_x_low:2:graph_x_high)
    xlim([graph_x_low graph_x_high])
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
    ylabel('|Displacement|')
    box off

    nexttile(2)
    plot(bin_centers,binned_properties_sub.(properties_to_plot{2}).mean); hold on;
    % if align_to_drink_onset==1
    %     xlabel('Time since arrival (s)')
    % elseif align_to_drink_offset == 1
    %     xlabel('Time till run (s)')
    %     rho1 = -1*rho1; rho2 = -1*rho2;
    % end
    ylabel(ylabels);

    % plot([graph_x_low graph_x_high],[(90) (90)],'--','color',[0.5 0.5 0.5],'LineWidth',2); hold on;
    xticks(graph_x_low:2:graph_x_high)
    xlim([graph_x_low graph_x_high])
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
    ylabel('|Displacement|')

    % [exp_tr,gof_tr,z] = fit(bin_centers(bin_centers <= 10),binned_properties_sub.meanAngDisplacement_pastPath.mean(bin_centers<=10),"exp2")
    % all_coefficients(n,:) = coeffvalues(exp_tr);
    % hold on;
    % plot(exp_tr,bin_centers(bin_centers<=10),binned_properties_sub.meanAngDisplacement_pastPath.mean(bin_centers<=10))

    % figure()
    % fit_exponential_curve
    % all_coefficients(n,:) = coefficients;

    nexttile(3)
    % plot(bin_centers,binned_properties_sub.(properties_to_plot{1}).mean - binned_properties_sub.(properties_to_plot{2}).mean); hold on;
    plot(bin_centers,binned_properties_sub.future_minus_past_angle.mean); hold on;
    % exp_tr(x) = a*exp(b*x)
    %plot(exp_tr.a)


    ind_of_crossing = find(binned_properties_sub.future_minus_past_angle.mean>0,1,'first');
    time_of_crossing(n) = bin_centers(ind_of_crossing);


    if align_to_drink_onset==1
        xlabel('Time since arrival (s)')
    elseif align_to_drink_offset == 1
        xlabel('Time till run (s)')
        rho1 = -1*rho1; rho2 = -1*rho2;
    end
    ylabel(ylabels);
    plot([graph_x_low graph_x_high],[(0) (0)],'--','color',[0.5 0.5 0.5],'LineWidth',2); hold on;
    xticks(graph_x_low:2:graph_x_high)
    xlim([graph_x_low graph_x_high])
    xtickangle(0)
    box off
    ylabel('Difference')
    ylim([-40 40])
end


%%
Variables = {'angDisplacement_futPath','angDisplacement_pastPath','future_minus_past_angle'};
group = struct();

for n = 1:length(unique_n)
    for i = 1:length(Variables)
        group.(Variables{i})(n,:) = group_summary(n).(Variables{i}).mean;
    end
end

nexttile(1);
plot(bin_centers,nanmean(group.(Variables{1})),'k','LineWidth',2)
nexttile(2);
plot(bin_centers,nanmean(group.(Variables{2})),'k','LineWidth',2)
nexttile(3);
plot(bin_centers,nanmean(group.(Variables{3})),'k','LineWidth',2)
pvals_subject_is_n = nan(length(bin_centers),1);
for bin = 1:length(bin_centers)
    [pvals_subject_is_n(bin)] = signrank(group.future_minus_past_angle(:,bin));
end
plot(bin_centers(pvals_subject_is_n<0.05),ylim_deg_difference_plot(2)*ones(sum(pvals_subject_is_n<0.05),1),'.k')


ylim(ylim_deg_difference_plot)
ylim
yticks(-80:20:80)
%      plot(bin_centers(pval==1),repmat((50),[sum(pval==1),1]),'.k')
nexttile(4);
% This is plotting the GRAND mean and sem of ALL events, NOT the average
% of each rat:
%     h = shadedErrorBar(bin_centers,binned_properties.future_minus_past_angle.mean,...
%     binned_properties.future_minus_past_angle.sem,'lineprops','k'); hold on;
%     h.mainLine.LineWidth = 2;
%     h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
%     hold on
pvals_datapoint_is_n = nan(length(bin_centers),1);
for bin = 1:length(bin_centers)
    [pvals_datapoint_is_n(bin)] = signrank(binned_properties.future_minus_past_angle.data{bin});
end

h = shadedErrorBar(bin_centers,binned_properties.future_minus_past_angle.mean,...
    binned_properties.future_minus_past_angle.sem,'lineprops','k'); hold on;
h.mainLine.LineWidth = 2;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
hold on
ylim(ylim_deg_difference_plot)
plot(bin_centers(pvals_datapoint_is_n<0.05),ylim_deg_difference_plot(2)*ones(sum(pvals_datapoint_is_n<0.05),1),'.k')
yticks(-80:20:80)
xlim([0 10])
plot([graph_x_low graph_x_high],[(0) (0)],'--','color',[0.5 0.5 0.5],'LineWidth',2); hold on;


set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
set(gcf, 'Color', 'white','PaperPositionMode','auto');
saveas(gcf,fullfile(fig_path,[fig_title '_disp' num2str(replay_dispersionThr,2) grouping ]),'pdf')
saveas(gcf,fullfile(fig_path,[fig_title '_disp' num2str(replay_dispersionThr,2) grouping ]),'jpeg')
%%
fig2 = figure();
fig2.Position = [600 600 300 250];
tiledlayout(2,2)

nexttile(1)
shadedErrorBar(bin_centers,nanmean(group.(Variables{1})),nanstd(group.(Variables{1}))./sqrt(sum(~isnan(group.(Variables{1})))));
hold on;
xticks(graph_x_low:2:graph_x_high)
xlim([graph_x_low graph_x_high])
xtickangle(0)
ylim_deg = [45 135];
ylim ([(ylim_deg(1)) (ylim_deg(2))]);
yticks([(ylim_deg(1)) (90) (ylim_deg(2))])
yticklabels({num2str(ylim_deg(1)),num2str(90),num2str(ylim_deg(2))})
ylimit = gca().YLim; y_range = ylimit(2)-ylimit(1);
xlabel('Time since arrival')
hold on;
box off
ylabel('|Displacement|')

nexttile(2)
shadedErrorBar(bin_centers,nanmean(group.(Variables{2})),nanstd(group.(Variables{2}))./sqrt(sum(~isnan(group.(Variables{2})))))
hold on;
plot([graph_x_low graph_x_high],[(90) (90)],'--','color',[0.5 0.5 0.5],'LineWidth',2); hold on;
xticks(graph_x_low:2:graph_x_high)
xlim([graph_x_low graph_x_high])
xtickangle(0)
ylim_deg = [45 135];
ylim ([(ylim_deg(1)) (ylim_deg(2))]);
yticks([(ylim_deg(1)) (90) (ylim_deg(2))])
yticklabels({num2str(ylim_deg(1)),num2str(90),num2str(ylim_deg(2))})
ylimit = gca().YLim; y_range = ylimit(2)-ylimit(1);
xlabel('Time since arrival')
hold on;
box off
ylabel('|Displacement|')


nexttile(3)
shadedErrorBar(bin_centers,nanmean(group.future_minus_past_angle),nanstd(group.future_minus_past_angle)./sqrt(sum(~isnan(group.future_minus_past_angle))))
hold on;
h.mainLine.LineWidth = 2;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
hold on
% data = group.(Variables{2}) - group.(Variables{1});
ylim_deg_difference_plot = [-50 20];
ylim(ylim_deg_difference_plot)
plot(bin_centers(pvals_subject_is_n<0.05),ylim_deg_difference_plot(2)*ones(sum(pvals_subject_is_n<0.05),1),'.k')
yticks(-40:20:20)
xlim([0 10])

xlabel('Time since arrival')
ylabel('difference')

set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
set(gcf, 'Color', 'white','PaperPositionMode','auto');
saveas(gcf,fullfile(fig_path,[fig_title '_disp' num2str(replay_dispersionThr,2) grouping '_sem' ]),'pdf')
saveas(gcf,fullfile(fig_path,[fig_title '_disp' num2str(replay_dispersionThr,2) grouping '_sem' ]),'jpeg')
%%
figure('Position',[2274 1006 50 10])
jitter = rand(1,6);
jitter = jitter-0.5;
plot(ones(length(unique_n),1)+jitter',time_of_crossing,'ok');
ylim([0 10])
xlim([0 2])
hold on
plot([0.5 1.5],[median(time_of_crossing) median(time_of_crossing)],'-k')
set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
set(gcf, 'Color', 'white','PaperPositionMode','auto');
saveas(gcf,fullfile(fig_path,['time_of_crossing_'  grouping ]),'pdf')
