% Good examples:
% CM1: 20230228 ~25 min in

shiftSizeDecoding = 0.005;

% User will define the example stopping period.
load Behavior_Data
load Analysis_Information.mat
load Experiment_Information.mat
load binDecoding_02
load Position_Data
load('laser_state')

windowSizeDecoding = decoder_binDecoding.windowSizeDecoding;
sessionNum_decoder = 1;

decoding_timeBin_centers = mean(decoder_binDecoding.timeBins,2);
% truncate everything to the run
clear posterior_color
posterior = decoder_binDecoding(sessionNum_decoder).posterior;

% limit posterior to Segment 1
good_inds = find(decoding_timeBin_centers >= Experiment_Information.Segments(1).Times(1) & decoding_timeBin_centers <= Experiment_Information.Segments(1).Times(2));
posterior = posterior(good_inds,:);
decoding_timeBin_centers = decoding_timeBin_centers(good_inds);

posterior_size = size(decoder_binDecoding(sessionNum_decoder).posterior,2);
posterior_direction_1 = posterior(:,1:posterior_size/2)';
posterior_direction_2 =  posterior(:,posterior_size/2+1:end)';


%%
speedThr = 5;
sessionNum=1;
plot_spikeDensityPeak=1;

% load the replays from this session;
replay_weighted_r_thr = 0.6;
replay_coverage_thr = 0.2;
replay_posterior_diff_thr = 0.33;
min_time_into_stopping_period=0;
max_time_into_stopping_period = inf;

% all sdes:
events = 'spike_filtered';
sde_thr = 3;
ripple_thr = -inf;
use_duration_og = 1;
posterior_diff_thr = 0;
coverage_thr = -inf;
weighted_r_thr = -inf;

t_events.('sde') = load_replays_from_individual_session(events,use_duration_og,coverage_thr,...
    weighted_r_thr,posterior_diff_thr,sde_thr,ripple_thr,min_time_into_stopping_period,max_time_into_stopping_period);

% all ripples:
events = 'ripple_filtered';
sde_thr = -inf;
ripple_thr = 3;
use_duration_og = 1;
posterior_diff_thr = 0;
coverage_thr = -inf;
weighted_r_thr = -inf;

t_events.('ripples') = load_replays_from_individual_session(events,use_duration_og,coverage_thr,...
    weighted_r_thr,posterior_diff_thr,sde_thr,ripple_thr,min_time_into_stopping_period,max_time_into_stopping_period);
% spike filtered replays:
events = 'spike_filtered';
sde_thr = 3;
ripple_thr = -inf;
use_duration_og = 0;
posterior_diff_thr = replay_posterior_diff_thr;
coverage_thr = replay_coverage_thr;
weighted_r_thr = replay_weighted_r_thr;

sub_events = load_replays_from_individual_session(events,use_duration_og,coverage_thr,...
    weighted_r_thr,posterior_diff_thr,sde_thr,ripple_thr,min_time_into_stopping_period,max_time_into_stopping_period);
t_events.('reverse_replays_spike') = sub_events(sub_events.reverse_replay==1,:);
t_events.('forward_replays_spike') = sub_events(sub_events.forward_replay==1,:);
t_events.('congruent_replays_spike') = sub_events(sub_events.congruent_replay==1,:);
t_events.('incongruent_replays_spike') = sub_events(sub_events.incongruent_replay==1,:);
t_events.('forward_congruent_replays_spike') = sub_events(sub_events.forward_congruent_replay==1,:);
t_events.('forward_incongruent_replays_spike') = sub_events(sub_events.forward_incongruent_replay==1,:);
t_events.('reverse_congruent_replays_spike') = sub_events(sub_events.reverse_congruent_replay==1,:);
t_events.('reverse_incongruent_replays_spike') = sub_events(sub_events.reverse_incongruent_replay==1,:);


%load spike density
load zscored_sd_pyr
%load LFPs
load zscored_ripple_power

posterior_color(:,:,1) = (1-posterior_direction_1 + ones(size(posterior_direction_1)))/2;
posterior_color(:,:,2) = (1-posterior_direction_1 + 1-posterior_direction_2)/2;
posterior_color(:,:,3) = (ones(size(posterior_direction_1)) + 1-posterior_direction_2)/2;

%change saturation
HSV = rgb2hsv(posterior_color);
HSV(:,:,2) = HSV(:,:,2)*25;
posterior_color = hsv2rgb(HSV);


laser_state_sub = compute_dataTemporalConcatenation(laser_state,[decoding_timeBin_centers(1) decoding_timeBin_centers(end)]);
laser_state_sub = compute_dataInterpolation(laser_state_sub,decoding_timeBin_centers,[]);
Position_Data_sub = compute_dataTemporalConcatenation(Position_Data,[decoding_timeBin_centers(1) decoding_timeBin_centers(end)]);
Position_Data_sub = compute_dataInterpolation(Position_Data_sub,decoding_timeBin_centers,4);
zscored_ripple_power = compute_dataTemporalConcatenation(zscored_ripple_power,[decoding_timeBin_centers(1) decoding_timeBin_centers(end)]);
zscored_ripple_power = compute_dataInterpolation(zscored_ripple_power ,decoding_timeBin_centers,[]);
zscored_sd_sub = compute_dataTemporalConcatenation(zscored_sd,[decoding_timeBin_centers(1) decoding_timeBin_centers(end)]);
zscored_sd_sub = compute_dataInterpolation(zscored_sd_sub,decoding_timeBin_centers,[]);

% times_to_label = round(linspace(1,length(times),300));
% times = 1:length(Position_Data_sub);
% time_labels = Position_Data_sub(:,1)./30000;
% time_labels = seconds(time_labels-time_labels(1));
% time_labels.Format = 'mm:ss';
% time_labels = cellstr(time_labels);

x_inds = 1:length(Position_Data_sub);
x_times = Position_Data_sub./30000 - Position_Data_sub(1)./30000;

x_inds_to_label = round(linspace(1,length(x_inds),300));
time_labels = seconds(x_times(x_inds_to_label));
time_labels.Format = 'mm:ss';
time_labels = cellstr(time_labels);

%%
ripple_color = [115 82 68]./260;
spike_density_color = [194 150 130]./260;

figure('Position',[398 320 1027 603])
ax1 = subplot(2,1,1);
yyaxis left
L1 = plot(x_inds,Position_Data_sub(:,2),'Color','k','LineStyle','-','LineWidth',2);
yyaxis right
L2 = plot(x_inds,zscored_sd_sub(:,3),'Color',spike_density_color,'LineStyle','-');
hold on
L3 = plot(x_inds,zscored_ripple_power(:,7),'Color',ripple_color,'LineStyle','-','Marker','none');
ax1.YAxis(1).Limits = [0 size(posterior,2)];
ax1.YAxis(2).Limits = [-2 12];
yyaxis left
%ax1.YTick = [0 50 100 150 200 250];
xticks(x_inds_to_label);

xticklabels(time_labels)
ylabel('Position (cm)')
box off

ax1.YColor = [0.1500 0.1500 0.1500];
yyaxis right
% ax1.YTick = [0 2 4 6 8 10 12 14 16];
ylabel('Z-score')
% [~,hobj,~,~] = legend({'','spike density','ripple power'});
h = legend({'','spike density','ripple power'});
% hl = findobj(hobj,'type','line');
% set(hl,'LineWidth',2);
yyaxis left
set(gca,'FontSize',12)

ax2 = subplot(2,1,2);
x_inds_to_label = round(linspace(1,length(decoding_timeBin_centers),300));
imagesc(posterior_color);
set(gca,'YDir','normal')
xticks(x_inds_to_label);
xticklabels(time_labels)
ax2.YLim = [0 size(posterior_color,1)];
% ax2.YTick = [0 25 50 75 100 125];
ax2.YTick = ax1.YAxis(1).TickValues./2;
ax2.YTickLabel = char([cellstr(num2str(ax1.YAxis(1).TickValues'))]);
xlabel('Time since reward onset (s)')
ylabel('Position (cm)')
set(gca,'FontSize',12)
linkaxes([ax1,ax2],'x','y')
ZoomHandle = zoom(gcf);
set(ZoomHandle,'Motion','horizontal')


%% find replays within the time being shown

% ToDo: show time on x axis as time since drink onset
% inds_plotted = gca().XLim;
% inds_plotted = round(inds_plotted); % has to be an index
inds_plotted = [296093 299978];
start_time_plotted = Position_Data_sub(inds_plotted(1),1);
stop_time_plotted = Position_Data_sub(inds_plotted(2),1);
% Find the nearest stopping time
stopping_period_starts = Reward_Epoch_Time_Boundaries_speed_thresholded(:,1);
[~,stopping_period] = min(abs(stopping_period_starts-start_time_plotted(1)));

% Plot the stopping period, +/- 2 seconds on either end
num_secs_flanking_stopping_period = 3;
stopping_period_start_time = Reward_Epoch_Time_Boundaries_speed_thresholded(stopping_period,1);
stopping_period_end_time = Reward_Epoch_Time_Boundaries_speed_thresholded(stopping_period,2);


figure_start_time = stopping_period_start_time - num_secs_flanking_stopping_period*spikeSampRate;
figure_end_time = stopping_period_end_time + num_secs_flanking_stopping_period*spikeSampRate;

[~,figure_start_ind] = min(abs(figure_start_time-decoding_timeBin_centers));
[~,figure_end_ind] = min(abs(figure_end_time-decoding_timeBin_centers));
[~,figure_drink_start_ind] = min(abs(stopping_period_start_time-decoding_timeBin_centers));

% new_times_to_label
% Frame rate = 0.005;
diff_x_labels_in_sec = 2;
plot_frame_rate = mean(diff(Position_Data_sub));

new_x_inds_to_label = figure_drink_start_ind - (diff_x_labels_in_sec/0.005):(diff_x_labels_in_sec/0.005):figure_end_ind;

xStart = -1*diff_x_labels_in_sec;
dx = 2;
N = length(new_x_inds_to_label);
x = xStart+ (0:N-1)*dx;

new_time_labels = seconds(x);
new_time_labels.Format = 'mm:ss';
new_time_labels = cellstr(new_time_labels);

% New Figure
%%
% figure('Position',[350 474 1000 350])
figure('Position',[350 474 680 220])


tiledlayout(2,8,'Padding','compact')
% ax1 = nexttile(1,[1,7]);
ax1 = nexttile(1,[1,8]);

yyaxis right
figure_x_inds = figure_start_ind:figure_end_ind;
L1 = plot(figure_x_inds,Position_Data_sub(figure_x_inds,2),'Color','k','LineStyle','-','LineWidth',2);
yyaxis left
L2 = plot(figure_x_inds,zscored_sd_sub(figure_x_inds,2),'Color',[115 82 68]./260,'LineStyle','-');
hold on
L3 = plot(figure_x_inds,zscored_ripple_power(figure_x_inds,6),'Color',[194 150 130]./260,'LineStyle','-','Marker','none');
[val,reward_start_ind] = min(abs(Position_Data_sub(:,1) - Reward_Epoch_Time_Boundaries_speed_thresholded(stopping_period,1)));
hold on;
h = xline(reward_start_ind,':k');
h.LineWidth = 2;
[val,reward_stop_ind] = min(abs(Position_Data_sub(:,1) - Reward_Epoch_Time_Boundaries_speed_thresholded(stopping_period,2)));
h = xline(reward_stop_ind,':k');
h.LineWidth = 2;

ax1.XLim = [figure_x_inds(1) figure_x_inds(end)];
ax1.YAxis(2).Limits = [0 size(posterior,2)];
ax1.YAxis(1).Limits = [-2 12];

ax1.XTick = new_x_inds_to_label;
ax1.XTickLabel = new_time_labels;
yyaxis right
ylabel('Position (cm)')

ax1.YColor = [0.1500 0.1500 0.1500];
yyaxis left
ylabel('Z-score')
% [~,hobj,~,~] = legend({'','spike density','ripple power'},'box','off','FontSize',9);
% h = legend({'','spike density','ripple power'});
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',2);
yyaxis left
set(gca,'FontSize',8)
box off

ax2 = nexttile(9,[1,8]);
imagesc(posterior_color(:,figure_x_inds,:));
set(gca,'YDir','normal')
new_x_inds_to_label = find(ismember(figure_x_inds,new_x_inds_to_label)==1);
ax2.XTick = new_x_inds_to_label;
ax2.XTickLabel = new_time_labels;
ax2.YLim = [0 size(posterior_color,1)];
% ax2.YTick = [0 25 50 75 100 125];
ax2.YTick = ax1.YAxis(2).TickValues./2;
ax2.YTickLabel = char([cellstr(num2str(ax1.YAxis(2).TickValues'))]);
xlabel('Time since reward onset (s)')
ylabel('Position (cm)')
set(gca,'FontSize',8)

%%

% Find reverse replays in this stopping period
event_types = {'forward_replays_spike','reverse_replays_spike'};
for i = 1:2
    event_times = t_events.(event_types{i}).('timePoints_og')(:,1);
    stopping_period_events.(event_types{i}) = t_events.(event_types{i})(event_times > Reward_Epoch_Time_Boundaries_speed_thresholded(stopping_period,1) & event_times < Reward_Epoch_Time_Boundaries_speed_thresholded(stopping_period,2),:);
end

nexttile(1)
yyaxis left
y_high = ax1.YLim(2);
y_low = ax1.YLim(2)-30;
for i = 1:2
    for j = 1:height(stopping_period_events.(event_types{i}))
        event  = stopping_period_events.(event_types{i})(j,:);
        event_time = mean(event.timePoints);
        [val,event_ind] = min(abs(event_time-decoding_timeBin_centers));
        [normx, normy] = coord2norm(ax1, [event_ind event_ind], [y_high y_low]);
        ar = annotation('arrow',normx,normy); hold on;
        if i == 1
            ar.Color = [.4660 0.6740 0.1880];
        elseif i == 2
            ar.Color = [0.4940 0.1840 0.5560];
        end
    end
end

ytick_values = ax1.YAxis(1).TickValues;

% replay_num = 0;
% event_types = {'forward_replays_spike','reverse_replays_spike'};
% while replay_num<2
%     for i = 1:2
%         for j = 1:height(stopping_period_events.(event_types{i}))
%             replay_num = replay_num+1;
%             event  = stopping_period_events.(event_types{i})(j,:);
% 
%             posterior_left = cell2mat(event.posterior_full_left);
%             posterior_right = cell2mat(event.posterior_full_right);
% 
%             numTimeBins = size(posterior_left,1);
%             numSpatialBins = size(posterior_left,2);
% 
%             % Pull out replay metrics to add to plot:
% 
%             weighted_r_sub = event.weighted_r;
%             max_jump_distance_sub = event.max_jump_distance;
%             ave_jump_distance_sub = event.mean_jump_distance;
%             range_normalized_sub = event.range_normalized;
%             coverage_wu_sub = event.coverage_wu;
%             ave_sharpness_at_peak_sub = event.mean_sharpness_at_peak;
%             best_map_sub = event.best_map;
%             duration_sub = event.duration;
%             ratPos_sub = event.ratPos(1);
%             startsLocal_sub = event.startsLocal;
%             total_posterior_right_sub = event.total_posterior_right;
%             total_posterior_left_sub = event.total_posterior_left;
%             time_since_drink_onset_sub = event.time_since_drink_onset;
%             cum_coverage_sub = event.cum_coverage;
%             replay_ripple_power_sub = event.replay_ripple_power;
%             replay_spikeDensity_power_sub = event.replay_spikeDensity_power;
%             laser_state = event.laser_state;
% 
%             if i == 1
%             info_string0 = 'For';
%             elseif i ==2
%                 info_string0 = 'Rev';
%             end
% 
%             info_string1 = ['R ' num2str(weighted_r_sub,2)];
%             info_string2 = ['Cov ' num2str(coverage_wu_sub,2)];
% 
%             if ratPos_sub(1) < 1
%                 ratPos_sub(1) = 1;
%             end
%             if ratPos_sub(1) > size(posterior_left,2)
%                 ratPos_sub(1) = size(posterior_left,2);
%             end
% 
%             clear posterior_color_ex
%             posterior_color_ex(:,:,1) = (1-posterior_left' + ones(size(posterior_left')))/2;
%             posterior_color_ex(:,:,2) = (1-posterior_left' + 1-posterior_right')/2;
%             posterior_color_ex(:,:,3) = (ones(size(posterior_left')) + 1-posterior_right')/2;
% 
%             %change saturation
%             HSV = rgb2hsv(posterior_color_ex);
%             HSV(:,:,2) = HSV(:,:,2)*10;
%             posterior_color_ex = hsv2rgb(HSV);
% 
% 
%             if replay_num == 1
%                 ax = nexttile(8,[1,1]);
%             elseif replay_num == 2
%                 ax = nexttile(16,[1,1]);
%             end
%             imagesc(posterior_color_ex)
%             set(ax,'YDir','normal');
%             ax.YTick = ax1.YAxis(1).TickValues./2;
%             ax.YTickLabel = [];
%             textColor = 'k';
%             hold on
% 
%             ax.LineWidth = 1;
%             if i == 1
%                 ax.XColor = [.4660 0.6740 0.1880];
%                 ax.YColor = [.4660 0.6740 0.1880];
%             elseif i == 2
%                 ax.XColor = [0.4940 0.1840 0.5560];
%                 ax.YColor = [0.4940 0.1840 0.5560];
%             end
%             xticks([0.5 size(posterior_right,1)+0.5])
%             xticklabels({num2str(0) num2str(windowSizeDecoding + (size(posterior_right,1)-1)*shiftSizeDecoding,2)})
%             xlabel('Replay time (s)','FontSize',8)
%             ylabel('Position (cm)','FontSize',8)
% 
%             y_axis_range = ax.YLim(2) - ax.YLim(1);
% %             text(2,ax.YLim(1)+ 0.98*y_axis_range, info_string0, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Color',textColor,'fontsize',10)
% %             text(2,ax.YLim(1)+ 0.90*y_axis_range, info_string1, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Color',textColor,'fontsize',10)
% %             text(2,ax.YLim(1)+ 0.82*y_axis_range, info_string2, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Color',textColor,'fontsize',10)
% %             
%             %print(['stopping period' num2str(stopping_period) ' event' num2str(replay_num)],'-djpeg')
%         end
%     end
%end
set(gcf, 'Color', 'white','Renderer','painters');

export_fig(['Example_posterior_over_stopping_period_' num2str(stopping_period)],'-jpeg');
export_fig(['Example_posterior_over_stopping_period_' num2str(stopping_period)],'-pdf');

fig_path = '/home/caitlin/Insync/caitlinmallory@berkeley.edu/Google Drive/Processed_Data/Manuscript_Figures/';
set(gcf,'Renderer','painters','PaperPositionMode', 'auto');
print(['Example_posterior_over_stopping_period_' num2str(stopping_period)],'-dpdf','-bestfit')
saveas(gcf,fullfile(fig_path,['Example_posterior_over_stopping_period_' num2str(stopping_period)]),'pdf')

%saveas(gcf,['Example_posterior_over_stopping_period_' num2str(stopping_period)],'jpg')