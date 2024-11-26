align_to_drink_offset = 1;


colorlimits = [0 3]; %supp
%colorlimits = [0 8]; %main

% this code does not work for overlapping bins
sig_test = 'means'; %means, medians or shuffle

replay(isnan(replay.drink_period_number),:) = [];
[C,ia, ic] = unique(replay(:,{'unique_session_id','drink_period_number','home_event','rat_label','drink_period_time','laser_state_binary'}));
% [C,ia, ic] = unique(replay(:,{'unique_session_id','drink_period_number','home_event','rat_label','drink_period_time','laser_state_binary', 'angle_between_past_future_trajectory'}));


% category1 = 'home';
% category2 = 'away';
% inds1 = find(C.home_event==1);
% inds2 = find(C.home_event==0);

category1 = 'off';
category2 = 'on';
inds1 = find(C.laser_state_binary==0);
inds2 = find(C.laser_state_binary==1);

% category1 = 'tight';
% category2 = 'wide';
% inds1 = find(C.angle_between_past_future_trajectoyr >= 0 & C.angle_between_past_future_trajectory<=60);
% inds2 = find(C.angle_between_past_future_trajectory>=120 & C.angle_between_past_future_trajectory<=180);


plot_anti_past = 0;
smooth_rate_plots = 1;
smoothing_sigma = 1; % bins
filter_length = smoothing_sigma*6;
if(mod(filter_length,2)==0) % iseven
    filter_length = filter_length+1;
end


num_stopping_periods = height(C);
event_types = {'future','past','antiPast','replay','home_ending_future'};
bin_edges = start_time:windowSize:end_time;
bin_centers = (bin_edges(1:end-1)' + bin_edges(2:end)')./2;

combined_hists = struct();
combined_hists_rates = struct();

for i =1:length(event_types)
    combined_hists.(event_types{i}) = nan(num_stopping_periods,length(bin_edges)-1);
    combined_hists_rates.(event_types{i}) = nan(num_stopping_periods,length(bin_edges)-1);
end

binned_data_by_stopping_period.replay = nan(num_stopping_periods,length(bin_edges)-1);

for i = 1:num_stopping_periods
    % find the data in this unique stopping period:
    replay_sub = replay(replay.unique_session_id == C.unique_session_id(i) & replay.drink_period_number == C.drink_period_number(i),:);

    for n = 1:length(event_types)
        event_inds = find(replay_sub.(event_types{n})==1);
        if align_to_drink_onset==1
            event_hists =  histcounts(replay_sub.time_since_real_drink_onset(event_inds),bin_edges);
        elseif align_to_drink_offset==1
            event_hists =  histcounts(replay_sub.time_till_real_drink_offset(event_inds),bin_edges);
        end
        event_hists(replay_sub.drink_period_time(1) < bin_centers(:,1)) = nan;

        % length of stopping period:
        combined_hists.(event_types{n})(i,:) = event_hists;
    end
end

for i = 1:length(event_types)
    combined_hists_rates.(event_types{i}) = combined_hists.(event_types{i})./windowSize;
end

w = setUp_gaussFilt([1 filter_length ],smoothing_sigma);
%w = setUp_gaussFilt_sigma(smoothing_sigma,bin_width);
% w2 = setUp_gaussFilt_sigma(1,0.5)
% smooth the rates on each lap, if requested
if smooth_rate_plots==1
    for i = 1:length(event_types)
        for j = 1:size(combined_hists_rates.future,1)

            rate_sub = combined_hists_rates.(event_types{i})(j,:);
            no_nan_inds = find(~isnan(rate_sub));
            rate_sub(isnan(rate_sub)) = [];
            combined_hists_rates.(event_types{i})(j,no_nan_inds) = conv(rate_sub,w,'same');
            %combined_hists_rates.(event_types{i})(j,no_nan_inds) = conv(rate_sub,w,'valid');
        end
    end
end

%% Plot rates over the stopping period

color_replay = [0 0 0];
color_future = [.4660 0.6740 0.1880];
color_past = [0.4940 0.1840 0.5560];
colors = [color_future; color_past];

transparency_pcnt = 0.99;
colors_2 = [1 - transparency_pcnt*(1- colors(:,1)) 1 - transparency_pcnt*(1- colors(:,2)) 1 - transparency_pcnt*(1-colors(:,3))];

fig = figure();
fig.Units = 'centimeters';
fig.Position = [50.8000 26.2467 4.5 10];
tiledlayout(3,1,'TileSpacing','compact')

ax1 = nexttile();
ax1.FontSize = 8;

x_mean  = nanmean(combined_hists_rates.future);
x_err = nanstd(combined_hists_rates.future./sqrt((sum(~isnan(combined_hists_rates.future)))));

if align_to_drink_offset==1
    bin_centers = -1.*bin_centers;
end
h = shadedErrorBar(bin_centers,x_mean,...
    x_err,'lineprops','g'); hold on;
h.patch.FaceColor = colors_2(1,:);
h.mainLine.Color = colors(1,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';

x_mean  = nanmean(combined_hists_rates.past);
x_err = nanstd(combined_hists_rates.past./sqrt((sum(~isnan(combined_hists_rates.past)))));

h = shadedErrorBar(bin_centers,x_mean,...
    x_err,'lineprops','m'); hold on;
h.patch.FaceColor = colors_2(2,:);
h.mainLine.Color = colors(2,:);
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';

color_antiPast = [0 0 1];

if plot_anti_past == 1
    x_mean  = nanmean(combined_hists_rates.antiPast);
    x_err = nanstd(combined_hists_rates.antiPast./sqrt((sum(~isnan(combined_hists_rates.antiPast)))));

    h = shadedErrorBar(bin_centers,x_mean,...
        x_err,'lineprops','m'); hold on;
    h.patch.FaceColor = color_antiPast;
    h.patch.FaceAlpha = 0.4;
    h.mainLine.Color = color_antiPast;
    h.mainLine.LineWidth = 1;
    h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
    % ylim([0 0.4])
end
if align_to_drink_onset==1
    xticks(0:2:10)
    xlabel('Time since arrival (s)')
else
    xticks(-10:2:0)
    xlabel('Time till departure (s)')
end
xlim([0 10])
xtickangle(0);
ylabel('Events/s')


ax2 = nexttile();
ax2.FontSize = 8;
x_mean  = nanmean(combined_hists_rates.future - combined_hists_rates.past);
x_err = nanstd(combined_hists_rates.future - combined_hists_rates.past)./sqrt((sum(~isnan((combined_hists_rates.future - combined_hists_rates.past)))));

h = shadedErrorBar(bin_centers,x_mean,...
    x_err,'lineprops','m'); hold on;
h.patch.FaceColor = [0 0 0];
h.patch.FaceAlpha = 0.4;
h.mainLine.Color = [0 0 0];
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
% ylim([0 0.4])
% ylim([-0.06 0.06]);
if align_to_drink_onset==1
    xticks(0:2:10)
    xlim([0 10])
    xlabel('Time since arrival (s)')
else
    xticks(-10:2:0)
    xlabel('Time till departure (s)')
    xlim([-10 0])
end
xtickangle(0);
ylabel('Events/s')
hold on
yline(0)
ax3 = nexttile();
ax3.FontSize = 8;
color_replay = [0 0 1];
hold on
x_mean  = nanmean(combined_hists_rates.replay);
x_err = nanstd(combined_hists_rates.replay./sqrt((sum(~isnan(combined_hists_rates.replay)))));

h = shadedErrorBar(bin_centers,x_mean,...
    x_err,'lineprops','m'); hold on;
h.patch.FaceColor = color_replay;
h.patch.FaceAlpha = 0.4;
h.mainLine.Color = color_replay;
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
% ylim([0.2 1.2])
if align_to_drink_onset==1
    xticks(0:2:10)
    xlabel('Time since arrival (s)')
else
    xticks(-10:2:0)
    xlabel('Time till departure (s)')
end
xtickangle(0);
ylabel('Events/s')

set(gcf, 'Color', 'white','Renderer','painters');
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,fullfile(fig_path,'future_past_difference_replay_rate_by_stopping_period'),'jpg')

%%
x_mean  = nanmean(combined_hists_rates.home_ending_future);
x_err = nanstd(combined_hists_rates.home_ending_future./sqrt((sum(~isnan(combined_hists_rates.home_ending_future)))));

h = shadedErrorBar(bin_centers,x_mean,...
    x_err,'lineprops','m'); hold on;
%%
num_laps_to_plot = 60;
num_bins_to_plot = 20; % 10 seconds

past_time_in_stopping_period_by_lap = nan(num_laps_to_plot,num_bins_to_plot);
future_time_in_stopping_period_by_lap = nan(num_laps_to_plot,num_bins_to_plot);
green_colormap = customcolormap(linspace(0,1,5), {'#014419','#1c7735','#5aae60','#a6db9d','#d7f1d6'});
purple_colormap = customcolormap(linspace(0,1,5), {'#410149','#762a84','#9b6fac','#c1a5cd','#e7d4e8'});

for lap = 1:num_laps_to_plot
    inds_0 = find(C.drink_period_number==lap);

    future_time_in_stopping_period_by_lap(lap,:) = nansum(combined_hists.('future')(inds_0,1:num_bins_to_plot),1);
    past_time_in_stopping_period_by_lap(lap,:) = nansum(combined_hists.('past')(inds_0,1:num_bins_to_plot),1);

end


fig1 = figure();
fig1.Position = [600 600 300 250];
tiledlayout(2,2);



ax1=nexttile(1);
imagesc(future_time_in_stopping_period_by_lap);
box off
colormap(ax1,green_colormap)
caxis(colorlimits)
xticks (linspace(0.5,20.5,6))
xticklabels(num2cell(0:2:10))
xtickangle(0);
yticks([1,10,20,30,40,50,60]);
ylabel('Trial')
% title('Forward replays')
ax1.FontSize = 8;

ax3=nexttile(3);
imagesc(past_time_in_stopping_period_by_lap);
box off
colormap(ax3,purple_colormap)
caxis(colorlimits)
xticks (linspace(0.5,20.5,6))
xticklabels(num2cell(0:2:10))
%     if align_to_reward_onset==1 || align_to_stopping_period_start==1
%         xticklabels({'0'; '2'; '4'; '6'; '8'; '10'})
%         % xlabel('Time (s)')
%     elseif align_to_stopping_period_departure==1
%         set(gca,'xticklabel',round(-10:2:0))
%         % xlabel('Time (s)')
%     end
xtickangle(0);
yticks([1,10,20,30,40,50,60]);
ylabel('Trial')
%title('Reverse replays')
ax3.FontSize = 8;
xlabel('Time since arrival (s)')
ax2 = nexttile(2);
x_mean  = nanmean(combined_hists_rates.future);
x_err = nanstd(combined_hists_rates.future./sqrt((sum(~isnan(combined_hists_rates.future)))));

if align_to_drink_offset==1
    bin_centers = -1.*bin_centers;
end
h = shadedErrorBar(bin_centers,x_mean,...
    x_err,'lineprops','g'); hold on;
h.patch.FaceColor = color_future;
% h.patch.FaceAlpha = 0.4;
h.mainLine.Color = color_future;
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';

x_mean  = nanmean(combined_hists_rates.past);
x_err = nanstd(combined_hists_rates.past./sqrt((sum(~isnan(combined_hists_rates.past)))));
hold on;

h = shadedErrorBar(bin_centers,x_mean,...
    x_err,'lineprops','m'); hold on;
h.patch.FaceColor = color_past;
%h.patch.FaceAlpha = 0.4;
h.mainLine.Color = color_past;
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
xlim([0 10])
xticks(0:2:10)
xticklabels({''})
ylabel('Events/s')

ax4 = nexttile(4);
x_mean  = nanmean(combined_hists_rates.future - combined_hists_rates.past);
x_err = nanstd(combined_hists_rates.future - combined_hists_rates.past)./sqrt((sum(~isnan((combined_hists_rates.future - combined_hists_rates.past)))));

h = shadedErrorBar(bin_centers,x_mean,...
    x_err,'lineprops','m'); hold on;
h.patch.FaceColor = [0 0 0];
h.mainLine.Color = [0 0 0];
h.mainLine.LineWidth = 1;
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
% ylim([0 0.4])
% ylim([-0.06 0.06]);
if align_to_drink_onset==1
    xticks(0:2:10)
    xlabel('Time since arrival (s)')
else
    xticks(-10:2:0)
    xlabel('Time till departure (s)')
end
sig_bins_means = ttest2(combined_hists_rates.future,combined_hists_rates.past);
sig_bins_means(isnan(sig_bins_means))=0;
sig_bins_means = logical(sig_bins_means);

ylimit = gca().YLim;
plot(bin_centers(sig_bins_means),ylimit(2)*ones(sum(sig_bins_means)),'.k')
xlim([0 10])
xtickangle(0);
ylabel('Events/s')
hold on
yline(0)

set(gcf, 'Color', 'white','Renderer','painters');
set(gcf, 'PaperPositionMode', 'auto');
fig_path = '/home/caitlin/Insync/caitlinmallory@berkeley.edu/Google Drive/Processed_Data/Manuscript_Figures';
saveas(gcf,fullfile(fig_path,'Main_open_field_rate_summary'),'pdf')

%% Look at the overall rate of retrospective or prospective replay on each trial type.
plot_fig2_rate_summary
% Plot home prospective versus home retrospective:
plot_figure2_prospective_retrospective_and_diff

keyboard

% OLD PLOT:
transparency_pcnt = 1;
colors_2 = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
    [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];

num_bins = size(combined_hists_rates.future,2);
figure('Position',[1921 560 800 300])
tiledlayout(2,4)

ax1 = nexttile(1);
data = combined_hists_rates.future;
y_mean = nanmean(data);
y_sem = nanstd(data)./sqrt(sum(~isnan(data)));
H = shadedErrorBar(bin_centers,y_mean,y_sem,'lineprops','k');
H.patch.FaceColor = colors(1,:);
H.mainLine.Color =  colors(1,:);
H.mainLine.LineWidth = 1;
H.patch.FaceColor = colors_2(1,:);
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';
hold on
data = combined_hists_rates.past;
y_mean = nanmean(data);
y_sem = nanstd(data)./sqrt(sum(~isnan(data)));
H = shadedErrorBar(bin_centers,y_mean,y_sem,'lineprops','k');
H.patch.FaceColor = colors(2,:);
H.mainLine.Color =  colors(2,:);
H.mainLine.LineWidth = 1;
H.patch.FaceColor = colors_2(2,:);
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';
% ylim([0 0.2])
xlim([0 10]);
xticks(0:2:10);
xtickangle(0);
xlabel('Time since arrival (s)')
ylabel('Events/s')

ax7 = nexttile(7);
data = combined_hists_rates.future(inds1,:);
y_mean = nanmean(data);
y_sem = nanstd(data)./sqrt(sum(~isnan(data)));
H = shadedErrorBar(bin_centers,y_mean,y_sem,'lineprops','k');
H.patch.FaceColor = colors(1,:);
H.mainLine.Color =  colors(1,:);
H.mainLine.LineWidth = 1;
H.patch.FaceColor = colors_2(1,:);
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';
hold on
data = combined_hists_rates.past(inds1,:);
y_mean = nanmean(data);
y_sem = nanstd(data)./sqrt(sum(~isnan(data)));
H = shadedErrorBar(bin_centers,y_mean,y_sem,'lineprops','k');
H.patch.FaceColor = colors(2,:);
H.mainLine.Color =  colors(2,:);
H.mainLine.LineWidth = 1;
H.patch.FaceColor = colors_2(2,:);
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';
hold on
data = combined_hists_rates.future(inds2,:);
y_mean = nanmean(data);
y_sem = nanstd(data)./sqrt(sum(~isnan(data)));
H = shadedErrorBar(bin_centers,y_mean,y_sem,'lineprops','--k');
H.patch.FaceColor = colors(1,:);
H.mainLine.Color =  colors(1,:);
H.mainLine.LineWidth = 1;
H.patch.FaceColor = colors_2(1,:);
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';
hold on
data = combined_hists_rates.past(inds2,:);
y_mean = nanmean(data);
y_sem = nanstd(data)./sqrt(sum(~isnan(data)));
H = shadedErrorBar(bin_centers,y_mean,y_sem,'lineprops','--k');
H.patch.FaceColor = colors(2,:);
H.mainLine.Color =  colors(2,:);
H.mainLine.LineWidth = 1;
H.patch.FaceColor = colors_2(2,:);
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';
% ylim([0 0.4])
xlim([0 10]);
xticks(0:2:10);
xtickangle(0);
xlabel('Time since arrival (s)')
ylabel('Events/s')

% Plot home prospective versus away prospective
ax2 = nexttile(2);
data = combined_hists_rates.future(inds1,:);
y_mean = nanmean(data);
y_sem = nanstd(data)./sqrt(sum(~isnan(data)));
H = shadedErrorBar(bin_centers,y_mean,y_sem,'lineprops','k');
H.patch.FaceColor = colors(1,:);
H.mainLine.Color = colors(1,:);
H.mainLine.LineWidth = 1;
H.patch.FaceColor = colors_2(1,:);
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';
hold on
data = combined_hists_rates.future(inds2,:);
y_mean = nanmean(data);
y_sem = nanstd(data)./sqrt(sum(~isnan(data)));
H = shadedErrorBar(bin_centers,y_mean,y_sem,'lineprops','--k');
H.patch.FaceColor = colors(1,:);
H.mainLine.Color = colors(1,:);
H.mainLine.LineWidth = 1;
H.patch.FaceColor = colors_2(1,:);
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';
xlabel('Time since arrival (s)')
ylabel('Events/s')
%legend({'Future, off', 'Past, off', 'Future, on', 'Past, on'})
% ylim([0 0.4])
xlim([0 10]);
xticks(0:2:10);
xtickangle(0);
% Plot home retrospective versus away retrospective
ax3 = nexttile(3);
data = combined_hists_rates.past(inds1,:);
y_mean = nanmean(data);
y_sem = nanstd(data)./sqrt(sum(~isnan(data)));
H = shadedErrorBar(bin_centers,y_mean,y_sem,'lineprops','k');
H.patch.FaceColor = colors(2,:);
H.mainLine.Color = colors(2,:);
H.mainLine.LineWidth = 1;
H.patch.FaceColor = colors_2(2,:);
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';
hold on
data = combined_hists_rates.past(inds2,:);
y_mean = nanmean(data);
y_sem = nanstd(data)./sqrt(sum(~isnan(data)));
H = shadedErrorBar(bin_centers,y_mean,y_sem,'lineprops','--k');
H.patch.FaceColor = colors(2,:);
H.mainLine.Color = colors(2,:);
H.mainLine.LineWidth = 1;
H.patch.FaceColor = colors_2(2,:);
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';
xlabel('Time since arrival (s)')
ylabel('Events/s')
%legend({'Future, off', 'Past, off', 'Future, on', 'Past, on'})
% ylim([0 0.2])
xlim([0 10]);
xticks(0:2:10);
xtickangle(0);

nexttile(4)
hold on
data1 = combined_hists_rates.future(inds1,:) - combined_hists_rates.past(inds1,:);
y_mean = nanmean(data1);
y_sem = nanstd(data1)./sqrt(sum(~isnan(data1)));
H = shadedErrorBar(bin_centers,y_mean,y_sem,'lineprops','-k');
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';
H.mainLine.LineWidth = 1;
data2 = combined_hists_rates.future(inds2,:) - combined_hists_rates.past(inds2,:);
y_mean = nanmean(data2);
y_sem = nanstd(data2)./sqrt(sum(~isnan(data2)));
hold on;
H = shadedErrorBar(bin_centers,y_mean,y_sem,'lineprops','--k');
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';
H.mainLine.LineWidth = 1;
H.mainLine.LineWidth = 1;

hold on
xlabel('Time since arrival (s)')
ylabel('Events/s')
yline(0)
% box off
% ylim([-0.3 0.3])
xlim([0 10]);
xticks(0:2:10);
xtickangle(0);

diff_of_diffs = abs(nanmean(data1)-nanmean(data2));
sig_times = nan(size(combined_hists_rates.future,2),1);

if strcmp(sig_test,'shuffle') == 1
    shuffle_diff_of_diffs = nan(1000,size(diff_of_diffs,2));
    for i = 1:1000
        shuffled_inds1 = randperm(size(combined_hists_rates.future,1), size(inds1,1));
        shuffled_inds2 = setdiff(1:size(combined_hists_rates.future,1),shuffled_inds1);
        shuffle1_diff = combined_hists_rates.future(shuffled_inds1,:)-combined_hists_rates.past(shuffled_inds1,:);
        shuffle2_diff = combined_hists_rates.future(shuffled_inds2,:)-combined_hists_rates.past(shuffled_inds2,:);
        shuffle_diff_of_diffs(i,:) = abs(nanmean(shuffle1_diff)-nanmean(shuffle2_diff));
    end
    sig_times = diff_of_diffs>quantile(shuffle_diff_of_diffs,0.95);
elseif strcmp(sig_test,'means')
    for i = 1:size(combined_hists_rates.future,2)
        a = combined_hists_rates.future(inds1,i)-combined_hists_rates.past(inds1,i);
        b = combined_hists_rates.future(inds2,i)-combined_hists_rates.past(inds2,i);
        sig_times(i) = ttest2(a,b);
    end
elseif strcmp(sig_test,'medians')
    for i = 1:size(combined_hists_rates.future,2)
        a = combined_hists_rates.future(inds1,i)-combined_hists_rates.past(inds1,i);
        b = combined_hists_rates.future(inds2,i)-combined_hists_rates.past(inds2,i);
        [~,sig_times(i)] = ranksum(a,b);
    end
end
sig_times = logical(sig_times);


hold on
ylimit = gca().YLim;
plot(bin_centers(sig_times),ylimit(2)*ones(sum(sig_times),1),'.k')

% Option 2
% Plot home retrospective versus retrospective
ax5 = nexttile(5);
data = combined_hists_rates.future(inds1,:);
y_mean = nanmean(data);
y_sem = nanstd(data)./sqrt(sum(~isnan(data)));
H = shadedErrorBar(bin_centers,y_mean,y_sem,'lineprops','k');
H.patch.FaceColor = colors(1,:);
H.mainLine.Color = colors(1,:);
H.mainLine.LineWidth = 1;
H.patch.FaceColor = colors_2(1,:);
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';
hold on
data = combined_hists_rates.past(inds1,:);
y_mean = nanmean(data);
y_sem = nanstd(data)./sqrt(sum(~isnan(data)));
H = shadedErrorBar(bin_centers,y_mean,y_sem,'lineprops','k');
H.patch.FaceColor = colors(2,:);
H.mainLine.Color = colors(2,:);
H.mainLine.LineWidth = 1;
H.patch.FaceColor = colors_2(2,:);
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';
xlim([0 10]);
xticks(0:2:10);
xtickangle(0);
xlabel('Time since arrival (s)')
% Plot away retrospective versus retrospective
ax6 = nexttile(6);
data = combined_hists_rates.future(inds2,:);
y_mean = nanmean(data);
y_sem = nanstd(data)./sqrt(sum(~isnan(data)));
H = shadedErrorBar(bin_centers,y_mean,y_sem,'lineprops','--k');
H.patch.FaceColor = colors(1,:);
H.mainLine.Color = colors(1,:);
H.mainLine.LineWidth = 1;
H.patch.FaceColor = colors_2(1,:);
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';
hold on
data = combined_hists_rates.past(inds2,:);
y_mean = nanmean(data);
y_sem = nanstd(data)./sqrt(sum(~isnan(data)));
H = shadedErrorBar(bin_centers,y_mean,y_sem,'lineprops','--k');
H.patch.FaceColor = colors(2,:);
H.mainLine.Color = colors(2,:);
H.mainLine.LineWidth = 1;
H.patch.FaceColor = colors_2(2,:);
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';
xlim([0 10]);
xticks(0:2:10);
xtickangle(0);
xlabel('Time since arrival (s)')

ax8 = nexttile(8);
y_mean = nanmean(data1);
y_sem = nanstd(data1)./sqrt(sum(~isnan(data1)));
H = shadedErrorBar(bin_centers,y_mean,y_sem,'lineprops','-k');
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';
H.mainLine.LineWidth = 1;
y_mean = nanmean(data2);
y_sem = nanstd(data2)./sqrt(sum(~isnan(data2)));
hold on;
H = shadedErrorBar(bin_centers,y_mean,y_sem,'lineprops','--k');
H.edge(1).Color = 'none';
H.edge(2).Color = 'none';
H.mainLine.LineWidth = 1;
H.mainLine.LineWidth = 1;

hold on
xlabel('Time since arrival (s)')
ylabel('Events/s')
yline(0)
% box off
% ylim([-0.2 0.35])
xlim([0 10]);
xticks(0:2:10);
xtickangle(0);

hold on
ylimit = gca().YLim;
plot(bin_centers(sig_times),ylimit(2)*ones(sum(sig_times),1),'.k')

set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
set(gcf, 'Color', 'white','PaperPositionMode','auto');
saveas(gcf,fullfile(fig_path,['replay_rates_' category1 '_versus_' category2 '_DegThr_' num2str(thrForCategorization_include) 'sdeThr_' num2str(sd_thr) '_disp' num2str(replay_dispersionThr) '_' title_ending]),'jpg')
saveas(gcf,fullfile(fig_path,['replay_rates_' category1 '_versus_' category2 '_DegThr_' num2str(thrForCategorization_include) 'sdeThr_' num2str(sd_thr) '_disp' num2str(replay_dispersionThr) '_' title_ending]),'pdf')

%%
















%%
%%
% rats_to_include = [1 4 6 12 13 14];
% 
% figure()
% tiledlayout(10,2)
% drink_bins = [1 3];
% colors = parula(length(drink_bins)-1);
% subplot(2,2,1)
% for i = 1:length(drink_bins)-1
%     data_sub = combined_hists_rates.future(C.drink_period_number >= drink_bins(i) &  C.drink_period_number < drink_bins(i+1) & C.home_event == 1& ismember(C.rat_label,rats_to_include),:);
%     x_mean = nanmean(data_sub);
%     x_err = nanstd(data_sub)./sqrt(sum(~isnan(data_sub)));
%     plot(bin_centers,x_mean,'color',colors(i,:),'LineWidth',2)
%     hold on
% end
% title('Home, future replay rate')
% ylim([0 0.3])
% subplot(2,2,2)
% for i = 1:length(drink_bins)-1
% 
%     data_sub = combined_hists_rates.future(C.drink_period_number >= drink_bins(i) &  C.drink_period_number < drink_bins(i+1) & C.home_event == 0& ismember(C.rat_label,rats_to_include),:);
%     x_mean = nanmean(data_sub);
%     x_err = nanstd(data_sub)./sqrt(sum(~isnan(data_sub)));
%     plot(bin_centers,x_mean,'color',colors(i,:),'LineWidth',2)
%     hold on
% end
% ylim([0 0.3])
% title('Away, future replay rate')
% subplot(2,2,3)
% for i = 1:length(drink_bins)-1
% 
%     data_sub = combined_hists_rates.past(C.drink_period_number >= drink_bins(i) &  C.drink_period_number < drink_bins(i+1) & C.home_event == 1& ismember(C.rat_label,rats_to_include),:);
%     x_mean = nanmean(data_sub);
%     x_err = nanstd(data_sub)./sqrt(sum(~isnan(data_sub)));
%     plot(bin_centers,x_mean,'color',colors(i,:),'LineWidth',2)
%     hold on
% end
% ylim([0 0.3])
% title('Home, past replay rate')
% subplot(2,2,4)
% for i = 1:length(drink_bins)-1
% 
%     data_sub = combined_hists_rates.past(C.drink_period_number >= drink_bins(i) &  C.drink_period_number < drink_bins(i+1) & C.home_event == 0 & ismember(C.rat_label,rats_to_include),:);
%     x_mean = nanmean(data_sub);
%     x_err = nanstd(data_sub)./sqrt(sum(~isnan(data_sub)));
%     plot(bin_centers,x_mean,'color',colors(i,:),'LineWidth',2)
%     hold on
% end
% ylim([0 0.3])
% title('Away, past replay rate')
% figure()
% subplot(1,2,1)
% for i = 1:length(drink_bins)-1
%   
%     data_sub = combined_hists_rates.future(C.drink_period_number >= drink_bins(i) &  C.drink_period_number < drink_bins(i+1) & C.home_event == 1 & ismember(C.rat_label,rats_to_include),:) - ...
%         combined_hists_rates.past(C.drink_period_number >= drink_bins(i) &  C.drink_period_number < drink_bins(i+1) & C.home_event ==1 & ismember(C.rat_label,rats_to_include),:);
%     x_mean = nanmean(data_sub);
%     x_err = nanstd(data_sub)./sqrt(sum(~isnan(data_sub)));
%     plot(bin_centers,x_mean,'color',colors(i,:),'LineWidth',2)
%     hold on
% end
% ylim([-0.15 0.2])
% title('Home, past replay rate')
% yline(0)
% subplot(1,2,2)
% for i = 1:length(drink_bins)-1
%   
%     data_sub = combined_hists_rates.future(C.drink_period_number >= drink_bins(i) &  C.drink_period_number < drink_bins(i+1) & C.home_event == 0 & ismember(C.rat_label,rats_to_include),:) - ...
%         combined_hists_rates.past(C.drink_period_number >= drink_bins(i) &  C.drink_period_number < drink_bins(i+1) & C.home_event == 0 & ismember(C.rat_label,rats_to_include),:);
%     x_mean = nanmean(data_sub);
%     x_err = nanstd(data_sub)./sqrt(sum(~isnan(data_sub)));
%     plot(bin_centers,x_mean,'color',colors(i,:),'LineWidth',2)
%     hold on
% end
% ylim([-0.15 0.2])
% title('Away, past replay rate')
% yline(0)
% 
% %%
% rats_to_include = [1 4 6 12 13 14];
% 
% figure()
% tiledlayout(10,2)
% drink_bins = [1 21; 21 41]; % [31 33]
% subplot(2,2,1)
% for i = 1:2
%     data_sub = combined_hists_rates.future(C.drink_period_number >= drink_bins(i,1) &  C.drink_period_number < drink_bins(i,2) & C.home_event == 1 & ismember(C.rat_label,rats_to_include),:);
%     x_mean = nanmean(data_sub);
%     x_err = nanstd(data_sub)./sqrt(sum(~isnan(data_sub)));
%        h = shadedErrorBar(bin_centers,x_mean,...
%         x_err,'lineprops','k'); hold on;
%     h.patch.FaceColor = [0 0 0];
%     h.patch.FaceAlpha = 0.4;
%     h.mainLine.Color = [0 0 0];
%     h.mainLine.LineWidth = 1;
%     h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
%     hold on
% end
% subplot(2,2,2)
% for i = 1:2
%     data_sub = combined_hists_rates.past(C.drink_period_number >= drink_bins(i,1) &  C.drink_period_number < drink_bins(i,2) & C.home_event == 1 & ismember(C.rat_label,rats_to_include),:);
%     x_mean = nanmean(data_sub);
%     x_err = nanstd(data_sub)./sqrt(sum(~isnan(data_sub)));
%        h = shadedErrorBar(bin_centers,x_mean,...
%         x_err,'lineprops','k'); hold on;
%     h.patch.FaceColor = [0 0 0];
%     h.patch.FaceAlpha = 0.4;
%     h.mainLine.Color = [0 0 0];
%     h.mainLine.LineWidth = 1;
%     h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
%     hold on
% end
% subplot(2,2,3)
% for i = 1:2
%     data_sub = combined_hists_rates.future(C.drink_period_number >= drink_bins(i,1) &  C.drink_period_number < drink_bins(i,2) & C.home_event == 0 & ismember(C.rat_label,rats_to_include),:);
%     x_mean = nanmean(data_sub);
%     x_err = nanstd(data_sub)./sqrt(sum(~isnan(data_sub)));
%        h = shadedErrorBar(bin_centers,x_mean,...
%         x_err,'lineprops','k'); hold on;
%     h.patch.FaceColor = [0 0 0];
%     h.patch.FaceAlpha = 0.4;
%     h.mainLine.Color = [0 0 0];
%     h.mainLine.LineWidth = 1;
%     h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
%     hold on
% end
% subplot(2,2,4)
% for i = 1:2
%     data_sub = combined_hists_rates.past(C.drink_period_number >= drink_bins(i,1) &  C.drink_period_number < drink_bins(i,2) & C.home_event == 0 & ismember(C.rat_label,rats_to_include),:);
%     x_mean = nanmean(data_sub);
%     x_err = nanstd(data_sub)./sqrt(sum(~isnan(data_sub)));
%        h = shadedErrorBar(bin_centers,x_mean,...
%         x_err,'lineprops','k'); hold on;
%     h.patch.FaceColor = [0 0 0];
%     h.patch.FaceAlpha = 0.4;
%     h.mainLine.Color = [0 0 0];
%     h.mainLine.LineWidth = 1;
%     h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
%     hold on
% end

% tiledlayout(10,2)
% drink_bins = 1:6:60;
% for i = 1:length(drink_bins)-1
%     nexttile()
%     data_sub = combined_hists_rates.future(C.drink_period_number >= drink_bins(i) &  C.drink_period_number < drink_bins(i+1) & C.home_event == 1,:);
%     x_mean = nanmean(data_sub);
%     x_err = nanstd(data_sub)./sqrt(sum(~isnan(data_sub)));
%     h = shadedErrorBar(bin_centers,x_mean,...
%         x_err,'lineprops','b'); hold on;
%     h.patch.FaceColor = [0 0 1];
%     h.patch.FaceAlpha = 0.4;
%     h.mainLine.Color = [0 0 1];
%     h.mainLine.LineWidth = 1;
%     h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
%     hold on
%     data_sub = combined_hists_rates.future(C.drink_period_number >= drink_bins(i) &  C.drink_period_number < drink_bins(i+1) & C.home_event == 0,:);
%     x_mean = nanmean(data_sub);
%     x_err = nanstd(data_sub)./sqrt(sum(~isnan(data_sub)));
%     h = shadedErrorBar(bin_centers,x_mean,...
%         x_err,'lineprops','k'); hold on;
%     h.patch.FaceColor = [0 0 0];
%     h.patch.FaceAlpha = 0.4;
%     h.mainLine.Color = [0 0 0];
%     h.mainLine.LineWidth = 1;
%     h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
%     nexttile()
%     data_sub = combined_hists_rates.past(C.drink_period_number >= drink_bins(i) &  C.drink_period_number < drink_bins(i+1) & C.home_event == 1,:);
%     x_mean = nanmean(data_sub);
%     x_err = nanstd(data_sub)./sqrt(sum(~isnan(data_sub)));
%     h = shadedErrorBar(bin_centers,x_mean,...
%         x_err,'lineprops','b'); hold on;
%     h.patch.FaceColor = [0 0 1];
%     h.patch.FaceAlpha = 0.4;
%     h.mainLine.Color = [0 0 1];
%     h.mainLine.LineWidth = 1;
%     h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
%     hold on
%     data_sub = combined_hists_rates.past(C.drink_period_number >= drink_bins(i) &  C.drink_period_number < drink_bins(i+1) & C.home_event == 0,:);
%     x_mean = nanmean(data_sub);
%     x_err = nanstd(data_sub)./sqrt(sum(~isnan(data_sub)));
%     h = shadedErrorBar(bin_centers,x_mean,...
%         x_err,'lineprops','k'); hold on;
%     h.patch.FaceColor = [0 0 0];
%     h.patch.FaceAlpha = 0.4;
%     h.mainLine.Color = [0 0 0];
%     h.mainLine.LineWidth = 1;
%     h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
% end

%%
