set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'}) % ensure that figure text and axes are black.

unit_choice = 'Sessions';
%unit_choice = 'Subjects';
connect_dataPoints = 1;
min_number_events = 5;

events = t_replay(t_replay.congruent_with_rat_location==1 & t_replay.time_since_reward_zone_entry > 0 & t_replay.time_since_reward_zone_entry < 10,:);
if strcmp(unit_choice,'Subjects')
    units= unique(t_replay.rat_label);
    events.unit_id = events.rat_label;
elseif strcmp(unit_choice,'Sessions')
    units = unique(t_replay.unique_session_id);
    events.unit_id = events.unique_session_id;
end

events.group(events.laser_state==0 & events.direction==2) = 1;
events.group(events.laser_state==1 & events.direction==2) = 2;
events.group(events.laser_state==0 & events.direction==1) = 3;
events.group(events.laser_state==1 & events.direction==1) = 4;

all_events = struct();
all_events.reverse_off = events(events.group==1,:);
all_events.reverse_on = events(events.group==2,:);
all_events.forward_off = events(events.group==3,:);
all_events.forward_on = events(events.group==4,:);

events_to_compare = {'forward_off','reverse_off'};

reverse_color = [0.4940 0.1840 0.5560];
forward_color = [.4660 0.6740 0.1880];


%%
rev_min_for_time = nan(length(units),1);
rev_time_mean = nan(length(units),1);
rev_time_median = nan(length(units),1);
rev_n = nan(length(units),1);
for_time_mean = nan(length(units),1);
for_time_median = nan(length(units),1);
for_n = nan(length(units),1);

for i = 1:length(units)
    data1_rat_inds = find(all_events.(events_to_compare{1}).unit_id == units(i));
    data1=all_events.(events_to_compare{1}).time_since_reward_zone_entry(data1_rat_inds);

    data2_rat_inds = find(all_events.(events_to_compare{2}).unit_id == units(i));    
    data2=all_events.(events_to_compare{2}).time_since_reward_zone_entry(data2_rat_inds);
%     figure(fig1);
%     subplot(4,4,i)
%     ecdf(data1);
%     hold on
%     ecdf(data2);
%     drawnow();
    rev_n(i) = length(data2)
    for_n(i) = length(data1)

    if length(data2)>=min_number_events & length(data1)>=min_number_events
    rev_min_for_time(i) = mean(data2)-mean(data1);
    rev_time_mean(i) = mean(data2);
    rev_time_median(i) = median(data2);
    
    for_time_mean(i) = mean(data1);
    for_time_median(i) = median(data1);
    end
%     figure(fig2);
%     subplot(4,4,i)
%     bins = [0:1:10]
%     [counts,bins] = hist(data1,bins)
%     b1 = bar(bins,counts); hold on;
%     b1.FaceColor = 'g';
%     b1.FaceAlpha = 0.5;
% 
%     bins = [0:1:10]
%     [counts,bins] = hist(data2,bins)
%     b2 = bar(bins,counts); hold on;    
%     b2.FaceColor = 'm';
%     b2.FaceAlpha = 0.5
end

fig3 = figure();
fig3.Position = [600 600 80 100];
offset = 0.5;
jitter = 0.15*rand(length(units),1); jitter = jitter-mean(jitter);
if connect_dataPoints==1
for data_point = 1:length(for_time_median)
line([jitter(data_point)-offset jitter(data_point)],[for_time_median(data_point),rev_time_median(data_point)],'color',[0.5 0.5 0.5], 'LineWidth', 0.25); hold on;
end
end
scatter(jitter,rev_time_median,20,'MarkerFaceColor',reverse_color,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','none'); hold on;
line([-0.2 0.2],[nanmean(rev_time_median) nanmean(rev_time_median)],'color',reverse_color,'LineWidth', 1.5); hold on;
scatter(jitter-offset,for_time_median,20,'MarkerFaceColor',forward_color,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');

line([0-offset-0.2 0-offset+0.2],[nanmean(for_time_median) nanmean(for_time_median)],'color',forward_color,'LineWidth', 1.5); hold on;
xlim([(0-offset-0.3) 0.3])
ylim([0 10])
yticks(0:2:10)
xticks((gca().XLim(1) + gca().XLim(2))./2)
xticklabels({})
ylabel('Time since arrival (s)')
xlabel(unit_choice)
box off
set(gcf, 'PaperPositionMode', 'auto');
fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures/';
saveas(fig3,fullfile([fig_path,'median_time_forwards_reverse_replays_per_' unit_choice]),'pdf')

fig4 = figure();
fig4.Position = [600 600 80 100];
offset = 0.5;
jitter = 0.15*rand(length(units),1); jitter = jitter-mean(jitter);
if connect_dataPoints==1
for data_point = 1:length(for_time_mean)
line([jitter(data_point)-offset jitter(data_point)],[for_time_mean(data_point),rev_time_mean(data_point)],'color',[0.5 0.5 0.5], 'LineWidth', 0.25); hold on;
end
end
scatter(jitter,rev_time_mean,20,'MarkerFaceColor',reverse_color,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','none'); hold on;
line([-0.2 0.2],[nanmean(rev_time_mean) nanmean(rev_time_mean)],'color',reverse_color,'LineWidth', 1.5); hold on;
scatter(jitter-offset,for_time_mean,20,'MarkerFaceColor',forward_color,'MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');

line([0-offset-0.2 0-offset+0.2],[nanmean(for_time_mean) nanmean(for_time_mean)],'color',forward_color,'LineWidth', 1.5); hold on;
xlim([(0-offset-0.3) 0.3])
ylim([0 10])
yticks(0:2:10)
xticks((gca().XLim(1) + gca().XLim(2))./2)
xticklabels({})
ylabel('Time since arrival (s)')
xlabel(unit_choice)
box off
set(gcf, 'PaperPositionMode', 'auto');
fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures/';
saveas(fig4,fullfile([fig_path,'mean_time_forwards_reverse_replays_per_' unit_choice]),'pdf')


[p,h,z] = signrank(rev_time_median,for_time_median);

