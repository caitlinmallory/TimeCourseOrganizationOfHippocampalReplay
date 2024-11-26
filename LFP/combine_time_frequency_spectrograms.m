


%fig_path = '/home/caitlin/Insync/caitlinmallory@berkeley.edu/Google Drive/Processed_Data/lfp_properties/LFP_over_stopping_period';
 fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures';
time_of_gamma_burst = 0;
time_of_peak_forward_rate = 1.75;
time_of_peak_reverse_rate = 3.75; %

reverse_window = [3 10];
forward_window = [0 3];

duration_thr = 10; %0, 6, or 10s. This is the minimum length of the stopping periods included

plot_laser_on_data = 0;
plot_regions_of_interest = 0;
mark_forward_reverse_times = 1;
num_yticks = 24;
rats = [1:11];
rats = [1 4 6 12 13 14]
reward_size = [0 inf];
novel_range = [10.2 10.5];
flags.flags_to_include = [11 ];
flags.flags_to_exclude = [1002 1003 0 101 14 18];
% 18 = dreadds; 19 = saline; 13 = laser off; 14 = laser on
must_include_all_flags = 1; % 1: requires that all flags in Flags_to_include are present. 0: requires that at least one flag in Flags_to_include are present.

smooth_frequency_band_plots = 0;
smooth_speed_plots = 0;
smoothing_sigma = 0.5; % seconds
min_freq = 1;
max_freq = 251;
num_frex = 250;
linear = 0; % 1 = linear, 0 = log spaced frequencies

delta_band = [1 6];
theta_band = [7 10];
beta_band = [11 30];
low_gamma_band = [30 50];
medium_gamma_band = [50 80];
high_gamma_band = [90 140];
ripple_band = [150 250];



if linear == 1
    frex = linspace(min_freq,max_freq,num_frex);
elseif linear == 0
    frex = logspace(log10(min_freq),log10(max_freq),num_frex);
end
freq_inds.delta = find(frex >=delta_band(1) & frex <= delta_band(2));
freq_inds.theta = find(frex >=theta_band(1) & frex <= theta_band(2));
freq_inds.beta = find(frex >= beta_band(1) & frex <= beta_band(2));
freq_inds.low_gamma = find(frex >= low_gamma_band(1) & frex <= low_gamma_band(2));
freq_inds.medium_gamma = find(frex >= medium_gamma_band(1) & frex <= medium_gamma_band(2));
freq_inds.high_gamma = find(frex >= high_gamma_band(1) & frex <= high_gamma_band(2));
freq_inds.ripple = find(frex >= ripple_band(1) & frex <= ripple_band(2));
frequency_bands = fieldnames(freq_inds);

lfp_sampling_rate = 1500;
all_sessions_speed_off = [];
all_sessions_speed_on = [];
all_sessions_tf = [];
all_sessions_tf_off = [];
all_sessions_tf_on = [];
all_sessions_tf_zscored_within_trial = [];
all_sessions_tf_off_zscored_within_trial = [];
all_sessions_tf_on_zscored_within_trial = [];
all_sessions_tf_zscored_within_session = [];
all_sessions_tf_off_zscored_within_session = [];
all_sessions_tf_on_zscored_within_session = [];


time_freq_band = struct();
time_freq_groups = {'tf','tf_off','tf_on','tf_zscored_within_trial','tf_off_zscored_within_trial','tf_on_zscored_within_trial',...
    'tf_zscored_within_session','tf_off_zscored_within_session','tf_on_zscored_within_session'};
for group = 1:length(time_freq_groups)
    for frequency_band = 1:length(frequency_bands)
        time_freq_band.(time_freq_groups{group}).(frequency_bands{frequency_band}) = [];
    end
end


for rat = rats
    %load_session_list_lfp
load_session_list_lfp_open_field
    for session = 1:length(dayFiles)
        cd(fullfile(directory,dayFiles{session}))

        load Experiment_Information
        if ~(Experiment_Information.reward_size >= reward_size(1) && Experiment_Information.reward_size < reward_size(2))
            continue
        end


        sessions_that_meet_criterion_day = check_for_relevant_sessions(flags,1,Experiment_Information,must_include_all_flags,1,novel_range);

        for segment =  1:length(sessions_that_meet_criterion_day)
            if duration_thr==0
                load(['tf_spectrograms_segment' num2str(sessions_that_meet_criterion_day(segment))]);
            elseif duration_thr==6
                load(['tf_spectrograms_low_duration_thr_segment' num2str(sessions_that_meet_criterion_day(segment))]);
            elseif duration_thr == 10
                load(['tf_spectrograms_high_duration_thr_segment' num2str(sessions_that_meet_criterion_day(segment))]);
            end
            stepSize = mean(diff(times2keep));
            w = setUp_gaussFilt_sigma(smoothing_sigma,stepSize);

            if ismember(19,Experiment_Information.Segments(segment).Flags)
                epoch_laser_state = zeros(size(epoch_laser_state));
            elseif ismember(18,Experiment_Information.Segments(segment).Flags)
           
                 epoch_laser_state = ones(size(epoch_laser_state));
                 tf_on = tf_off;
                 tf_on_zscored_within_trial = tf_off_zscored_within_trial;
                 tf_on_zscored_within_session = tf_off_zscored_within_session;
            end

            speed_off = speed_during_stopping(epoch_laser_state==0,:);
            speed_on = speed_during_stopping(epoch_laser_state==1,:);
            if smooth_speed_plots == 1
                for row = 1:size(speed_off,1)
                    speed_off(row,:) = conv(speed_off(row,:),w,'same');
                end
                for row = 1:size(speed_on)
                    speed_on(row,:) = conv(speed_on(row,:),w,'same');
                end
            end
            all_sessions_speed_off = [all_sessions_speed_off; nanmean(speed_off,1)];
            all_sessions_speed_on = [all_sessions_speed_on; nanmean(speed_on,1)];
%      
%             if length(find(isnan(tf_off))) > 0
%                 keyboard
%             end

            %             all_sessions_tf = cat(3,all_sessions_tf,tf);
            all_sessions_tf_off = cat(3,all_sessions_tf_off,tf_off);
            all_sessions_tf_on = cat(3,all_sessions_tf_on,tf_on);

            %             all_sessions_tf_zscored_within_trial = cat(3,all_sessions_tf_zscored_within_trial,tf_zscored_within_trial);
            all_sessions_tf_off_zscored_within_trial = cat(3,all_sessions_tf_off_zscored_within_trial,tf_off_zscored_within_trial);
            all_sessions_tf_on_zscored_within_trial = cat(3,all_sessions_tf_on_zscored_within_trial,tf_on_zscored_within_trial);

            %             all_sessions_tf_zscored_within_session = cat(3,all_sessions_tf_zscored_within_session,tf_zscored_within_session);
            all_sessions_tf_off_zscored_within_session = cat(3,all_sessions_tf_off_zscored_within_session,tf_off_zscored_within_session);
            all_sessions_tf_on_zscored_within_session = cat(3,all_sessions_tf_on_zscored_within_session,tf_on_zscored_within_session);

            time_freq = struct();
            time_freq.tf = tf;
            time_freq.tf_off = tf_off;
            time_freq.tf_on = tf_on;
            %             time_freq.tf_zscored_within_trial = tf_zscored_within_trial;
            time_freq.tf_off_zscored_within_trial = tf_off_zscored_within_trial;
            time_freq.tf_on_zscored_within_trial = tf_on_zscored_within_trial;
            %             time_freq.tf_zscored_within_session = tf_zscored_within_session;
            time_freq.tf_off_zscored_within_session = tf_off_zscored_within_session;
            time_freq.tf_on_zscored_within_session = tf_on_zscored_within_session;
            time_freq_groups = fieldnames(time_freq);


            for group = 1:length(time_freq_groups)
                for frequency_band = 1:length(frequency_bands)
                    session_power_in_band = nanmean(time_freq.(time_freq_groups{group})(freq_inds.(frequency_bands{frequency_band}),:));;
                    if smooth_frequency_band_plots == 1
                        session_power_in_band = conv(nanmean(session_power_in_band,1),w,'same');
                    end
                    time_freq_band.(time_freq_groups{group}).(frequency_bands{frequency_band}) = [time_freq_band.(time_freq_groups{group}).(frequency_bands{frequency_band}); session_power_in_band];
                end
            end
        end
    end
end
%%
time_freq_band_mean = struct();
time_freq_band_sem = struct();
for group = 1:length(time_freq_groups)
    for frequency_band = 1:length(frequency_bands)
        time_freq_band_mean.(time_freq_groups{group}).(frequency_bands{frequency_band}) = nanmean(time_freq_band.(time_freq_groups{group}).(frequency_bands{frequency_band}),1);
        time_freq_band_sem.(time_freq_groups{group}).(frequency_bands{frequency_band}) = nanstd(time_freq_band.(time_freq_groups{group}).(frequency_bands{frequency_band}),[],1)./...
            sqrt(height(time_freq_band.(time_freq_groups{group}).(frequency_bands{frequency_band})));
    end
end
signal_length = length(times2keep);

time_freq_band_mean = struct();
time_freq_band_sem = struct();
for group = 1:length(time_freq_groups)
    for frequency_band = 1:length(frequency_bands)
        time_freq_band_mean.(time_freq_groups{group}).(frequency_bands{frequency_band}) = nanmean(time_freq_band.(time_freq_groups{group}).(frequency_bands{frequency_band}),1);
        time_freq_band_sem.(time_freq_groups{group}).(frequency_bands{frequency_band}) = nanstd(time_freq_band.(time_freq_groups{group}).(frequency_bands{frequency_band}),[],1)./...
            sqrt(height(time_freq_band.(time_freq_groups{group}).(frequency_bands{frequency_band})));
    end
end
signal_length = length(times2keep);

%%
regions_to_outline = [];
regions_to_outline(1).freq_low = dsearchn(frex',delta_band(1));
regions_to_outline(1).freq_high = dsearchn(frex',delta_band(2));
regions_to_outline(1).time_start = dsearchn(times2keep',time_of_peak_reverse_rate-1.5);
regions_to_outline(1).time_end = dsearchn(times2keep',time_of_peak_reverse_rate+1.5);
regions_to_outline(2).freq_low = dsearchn(frex',low_gamma_band(1));
regions_to_outline(2).freq_high = dsearchn(frex',low_gamma_band(2));
regions_to_outline(2).time_start = dsearchn(times2keep',time_of_peak_forward_rate-1);
regions_to_outline(2).time_end = dsearchn(times2keep',time_of_peak_forward_rate+1);
regions_to_outline(3).freq_low = dsearchn(frex',medium_gamma_band(1));
regions_to_outline(3).freq_high = dsearchn(frex',medium_gamma_band(2));
regions_to_outline(3).time_start = dsearchn(times2keep',time_of_gamma_burst-0.5);
regions_to_outline(3).time_end = dsearchn(times2keep',time_of_gamma_burst + 0.5);



% plot results
figure('Position',[590 46 500 350])
% subplot(2,2,1)
% imagesc(log(nanmean(all_sessions_tf_off,3))); set(gca,'YDir','normal'); colormap jet;
% xtick_positions = dsearchn(times2keep',[-2, 0, 2, 4, 6, 8, 10]');
% set(gca,'xtick',xtick_positions,'xticklabel',round(-2:2:10))
% if linear==1
%     set(gca,'ytick',linspace(1,num_frex,num_yticks),'yticklabel',round(linspace(min_freq,max_freq,num_yticks)))
% elseif linear == 0
%     set(gca,'ytick',linspace(1,num_frex,num_yticks),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),num_yticks)*10)/10)
% end
% xlim([xtick_positions(1) xtick_positions(end)])
% colorbar
% xlabel('Time since drink onset (s)')
% ylabel('f[Hz]')
% title('Laser-OFF trials')
% drawnow()

% subplot(2,2,2)
% imagesc(log(nanmean(all_sessions_tf_on,3))); set(gca,'YDir','normal'); colormap jet;
% xtick_positions = dsearchn(times2keep',[-2, 0, 2, 4, 6, 8, 10]');
% set(gca,'xtick',xtick_positions,'xticklabel',round(-2:2:10))
% if linear==1
%     set(gca,'ytick',linspace(1,num_frex,num_yticks),'yticklabel',round(linspace(min_freq,max_freq,num_yticks)))
% elseif linear == 0
%     set(gca,'ytick',linspace(1,num_frex,num_yticks),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),num_yticks)*10)/10)
% end
% xlim([xtick_positions(1) xtick_positions(end)])
% colorbar
% xlabel('Time since drink onset (s)')
% ylabel('f[Hz]')
% title('Laser-ON trials')
% drawnow ()

subplot(2,2,3)
imagesc(nanmean(all_sessions_tf_off_zscored_within_trial,3)); set(gca,'YDir','normal','clim',[-1 1]); colormap jet;
xtick_positions = dsearchn(times2keep',[-2, 0, 2, 4, 6, 8, 10]');
set(gca,'xtick',xtick_positions,'xticklabel',round(-2:2:10))
if linear==1
    set(gca,'ytick',linspace(1,num_frex,num_yticks),'yticklabel',round(linspace(min_freq,max_freq,num_yticks)))
elseif linear == 0
    set(gca,'ytick',linspace(1,num_frex,num_yticks),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),num_yticks)*10)/10)
end
xlim([xtick_positions(1) xtick_positions(end)])
% clim([-1 1])
clim([-0.5 0.5])
hcb = colorbar; hcb.Label.String = 'Power (z-score)'; hcb.FontSize = 11;
xlabel('Time since drink onset (s)')
ylabel('f[Hz]')
title('Laser-OFF trials')

% Draw boxes around frequencies/regions of interest if requested
if plot_regions_of_interest==1
    for region = 1:size(regions_to_outline,2)
        rectangle('Position',[regions_to_outline(region).time_start, regions_to_outline(region).freq_low,...
            regions_to_outline(region).time_end - regions_to_outline(region).time_start,...
            regions_to_outline(region).freq_high - regions_to_outline(region).freq_low],...
            'EdgeColor','k','LineWidth',2,'LineStyle',':');
    end
end



subplot(2,2,4)
imagesc(nanmean(all_sessions_tf_on_zscored_within_trial,3)); set(gca,'YDir','normal','clim',[-1 1]); colormap jet;
xtick_positions = dsearchn(times2keep',[-2, 0, 2, 4, 6, 8, 10]');
set(gca,'xtick',xtick_positions,'xticklabel',round(-2:2:10))
if linear==1
    set(gca,'ytick',linspace(1,num_frex,num_yticks),'yticklabel',round(linspace(min_freq,max_freq,num_yticks)))
elseif linear == 0
    set(gca,'ytick',linspace(1,num_frex,num_yticks),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),num_yticks)*10)/10)
end
clim([-1 1])
xlim([xtick_positions(1) xtick_positions(end)])
hcb = colorbar; hcb.Label.String = 'Power (z-score)'; hcb.FontSize = 11;
xlabel('Time since drink onset (s)')
ylabel('f[Hz]')
title('Laser-ON trials')
% Draw boxes around frequencies/regions of interest if requested
if plot_regions_of_interest==1
    for region = 1:size(regions_to_outline,2)
        rectangle('Position',[regions_to_outline(region).time_start, regions_to_outline(region).freq_low,...
            regions_to_outline(region).time_end - regions_to_outline(region).time_start,...
            regions_to_outline(region).freq_high - regions_to_outline(region).freq_low],...
            'EdgeColor','k','LineWidth',2,'LineStyle',':');
    end
end

set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
saveas(gcf,fullfile(fig_path,['tf_spectrum','_reward_' num2str(reward_size(1),2) '_to_' num2str(reward_size(2),2)]),'pdf')

%%
empty_sessions = 0;
for i = 1:size(all_sessions_tf_off,3)
    if sum(sum(all_sessions_tf_off(:,:,i)))==0
        empty_session = 1;
    else
        empty_session = 0;
    end
    empty_sessions = empty_sessions+empty_session;
end

%% Look at speed during stopping period

figure('Position',[1963 907 900 125])
tiledlayout(1,7,'TileSpacing','tight')
mean_speed_off = nanmean(all_sessions_speed_off,1);
sem_speed_off = nanstd(all_sessions_speed_off,1)./sqrt(sum(~isnan(all_sessions_speed_off),1));
mean_speed_on = nanmean(all_sessions_speed_on,1);
sem_speed_on = nanstd(all_sessions_speed_on,1)./sqrt(sum(~isnan(all_sessions_speed_on),1));
nexttile
shadedErrorBar(times2keep,mean_speed_off,sem_speed_off)
if plot_laser_on_data==1
    shadedErrorBar(times2keep,mean_speed_on,sem_speed_on,'lineprops','r');
end
xtick_positions = [-2, 0, 2, 4, 6, 8, 10];
set(gca,'xtick',xtick_positions,'xticklabel',round(-2:2:10))
xtickangle(0);
xlim([-2 10])
ylim([0 60])

ylabel('Speed (cm/s)')
title('speed')
set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
saveas(gcf,fullfile(fig_path,['movement_speed_reward_' num2str(reward_size(1),2) '_to_' num2str(reward_size(2),2)]),'pdf')



%% zscored within session

frequency_bands_to_plot = {'ripple','high_gamma','medium_gamma','low_gamma','beta','theta','delta'};
%frequency_bands_to_plot = {'beta','low_gamma','delta','theta','high_gamma','medium_gamma'};
figure('Position',[1963 907 900 125])
tiledlayout(1,length(frequency_bands_to_plot))
for frequency = 1:length(frequency_bands_to_plot)
    nexttile
    shadedErrorBar(times2keep,time_freq_band_mean.tf_off_zscored_within_session.(frequency_bands_to_plot{frequency}),time_freq_band_sem.tf_off_zscored_within_session.(frequency_bands_to_plot{frequency}))
    xtick_positions = [-2, 0, 2, 4, 6, 8, 10];
    set(gca,'xtick',xtick_positions,'xticklabel',round(-2:2:10))
    xlabel('Time since arrival (s)');
    ylabel('z-scored within session')
    hold on
    if plot_laser_on_data == 1
        shadedErrorBar(times2keep,time_freq_band_mean.tf_on_zscored_within_session.(frequency_bands_to_plot{frequency}),time_freq_band_sem.tf_on_zscored_within_session.(frequency_bands_to_plot{frequency}),'lineprops','r');
        if frequency == length(frequency_bands_to_plot)
            h = legend({'laser-OFF','laser-ON','',''},'box','off');
        end
    end
    xlim([-2 10])
    hold on
    if mark_forward_reverse_times==1
    v1 = xline(time_of_peak_forward_rate); v1.Color = [.4660 0.6740 0.1880]; v1.LineWidth = 2;
    v2 = xline(time_of_peak_reverse_rate); v2.Color = [0.4940 0.1840 0.5560]; v2.LineWidth = 2;
    end
    title(frequency_bands_to_plot{frequency},'Interpreter','none')
end

set(gcf, 'Color', 'white','Renderer','painters');
%export_fig(fullfile(fig_path,['power_in_bands_over_time_zscored_session','_reward_' num2str(reward_size(1),2) '_to_' num2str(reward_size(2),2)]),'-jpeg')

%% zscored within trial


figure('Position',[1963 907 900 125])
tiledlayout(1,length(frequency_bands_to_plot),'TileSpacing','tight')
for frequency = 1:length(frequency_bands_to_plot)
    nexttile
    shadedErrorBar(times2keep,time_freq_band_mean.tf_off_zscored_within_trial.(frequency_bands_to_plot{frequency}),time_freq_band_sem.tf_off_zscored_within_trial.(frequency_bands_to_plot{frequency}))
    xtick_positions = [-2, 0, 2, 4, 6, 8, 10];
    set(gca,'xtick',xtick_positions,'xticklabel',round(-2:2:10))
    xtickangle(0)
%     xlabel('Time since arrival (s)');
    if frequency==1
    ylabel('Power (Zscore)')
    end
    hold on
    if plot_laser_on_data == 1
        shadedErrorBar(times2keep,time_freq_band_mean.tf_on_zscored_within_trial.(frequency_bands_to_plot{frequency}),time_freq_band_sem.tf_on_zscored_within_trial.(frequency_bands_to_plot{frequency}),'lineprops','r');
        if frequency == length(frequency_bands_to_plot)
            h = legend({'laser-OFF','laser-ON','',''},'box','off');
        end
    end
    xlim([-2 10])
%     ylim([-0.5 0.5])
    hold on
    if mark_forward_reverse_times==1
    v1 = xline(time_of_peak_forward_rate); v1.Color = [.4660 0.6740 0.1880]; v1.LineWidth = 2;
    v2 = xline(time_of_peak_reverse_rate); v2.Color = [0.4940 0.1840 0.5560]; v2.LineWidth = 2;
    end

     title(frequency_bands_to_plot{frequency},'Interpreter','none')
end

set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
saveas(gcf,fullfile(fig_path,['power_in_bands_over_time_zcored_within_trial','_reward_' num2str(reward_size(1),2) '_to_' num2str(reward_size(2),2)]),'pdf')


%%
figure()
y = time_freq_band_mean.tf_off_zscored_within_trial.('theta');
plot(times2keep,(y-min(y))./(max(y)-min(y))); hold on;
y = time_freq_band_mean.tf_off_zscored_within_trial.('ripple');
plot(times2keep,(y-min(y))./(max(y)-min(y))); hold on;
y = time_freq_band_mean.tf_off_zscored_within_trial.('beta');
plot(times2keep,(y-min(y))./(max(y)-min(y))); hold on;
y = time_freq_band_mean.tf_off_zscored_within_trial.('delta');
plot(times2keep,(y-min(y))./(max(y)-min(y))); hold on;
xtick_positions = [-2, 0, 2, 4, 6, 8, 10];
legend({'theta','ripple','beta','delta'})
 
%% Raw


figure('Position',[1963 907 1500 250])
tiledlayout(1,length(frequency_bands_to_plot))
for frequency = 1:length(frequency_bands_to_plot)
    nexttile
    shadedErrorBar(times2keep,time_freq_band_mean.tf_off.(frequency_bands_to_plot{frequency}),time_freq_band_sem.tf_off.(frequency_bands_to_plot{frequency}))
    xtick_positions = [-2, 0, 2, 4, 6, 8, 10];
    set(gca,'xtick',xtick_positions,'xticklabel',round(-2:2:10))
    xlabel('time since reward onset (s)');
    ylabel('power, (uV^2)/hz')
    hold on
    if plot_laser_on_data == 1
        shadedErrorBar(times2keep,time_freq_band_mean.tf_on.(frequency_bands_to_plot{frequency}),time_freq_band_sem.tf_on.(frequency_bands_to_plot{frequency}),'lineprops','r');
        if frequency == length(frequency_bands_to_plot)
            h = legend({'laser-OFF','laser-ON','',''},'box','off');
        end
    end
    xlim([-2 10])
    hold on
    if mark_forward_reverse_times==1
    v1 = xline(time_of_peak_forward_rate); v1.Color = [.4660 0.6740 0.1880]; v1.LineWidth = 2;
    v2 = xline(time_of_peak_reverse_rate); v2.Color = [0.4940 0.1840 0.5560]; v2.LineWidth = 2;
    end

    title(frequency_bands_to_plot{frequency},'Interpreter','none')
end

set(gcf, 'Color', 'white','Renderer','painters');
export_fig(fullfile(fig_path,['power_in_bands_over_time_raw','_reward_' num2str(reward_size(1),2) '_to_' num2str(reward_size(2),2)]),'-jpeg')

%%
% Ratios

% figure()
% a = nanmean(time_freq_band.tf_off.low_gamma./time_freq_band.tf_off.high_gamma); hold on;
% b = nanmean(time_freq_band.tf_on.low_gamma./time_freq_band.tf_on.high_gamma); hold on;
% a_sem = nanstd(time_freq_band.tf_off.low_gamma./time_freq_band.tf_off.high_gamma)./sqrt(height(time_freq_band.tf_off.high_gamma));
% b_sem = nanstd(time_freq_band.tf_on.low_gamma./time_freq_band.tf_on.high_gamma)./sqrt(height(time_freq_band.tf_on.high_gamma));
% shadedErrorBar(x,a,a_sem,'lineprops','k'); hold on
% shadedErrorBar(x,b,b_sem,'lineprops','r')
%
% figure()
% a = nanmean(time_freq_band.tf_off.beta./time_freq_band.tf_off.delta); hold on;
% b = nanmean(time_freq_band.tf_on.beta./time_freq_band.tf_on.delta); hold on;
% a_sem = nanstd(time_freq_band.tf_off.beta./time_freq_band.tf_off.delta)./sqrt(height(time_freq_band.tf_off.delta));
% b_sem = nanstd(time_freq_band.tf_on.beta./time_freq_band.tf_on.delta)./sqrt(height(time_freq_band.tf_on.delta));
% shadedErrorBar(x,a,a_sem,'lineprops','k'); hold on
% shadedErrorBar(x,b,b_sem,'lineprops','r')
%
% figure()
% a = nanmean(time_freq_band.tf_off.low_gamma./time_freq_band.tf_off.delta); hold on;
% b = nanmean(time_freq_band.tf_on.low_gamma./time_freq_band.tf_on.delta); hold on;
% a_sem = nanstd(time_freq_band.tf_off.low_gamma./time_freq_band.tf_off.delta)./sqrt(height(time_freq_band.tf_off.delta));
% b_sem = nanstd(time_freq_band.tf_on.low_gamma./time_freq_band.tf_on.delta)./sqrt(height(time_freq_band.tf_on.delta));
% shadedErrorBar(times2keep,a,a_sem,'lineprops','k'); hold on
% shadedErrorBar(times2keep,b,b_sem,'lineprops','r')
% xlim([-2 10])
% title('low gamma:delta ratio')

% figure()
% a = nanmean(time_freq_band.tf_off.low_gamma./time_freq_band.tf_off.medium_gamma); hold on;
% b = nanmean(time_freq_band.tf_on.low_gamma./time_freq_band.tf_on.medium_gamma); hold on;
% a_sem = nanstd(time_freq_band.tf_off.low_gamma./time_freq_band.tf_off.medium_gamma)./sqrt(height(time_freq_band.tf_off.medium_gamma));
% b_sem = nanstd(time_freq_band.tf_on.low_gamma./time_freq_band.tf_on.medium_gamma)./sqrt(height(time_freq_band.tf_on.medium_gamma));
% shadedErrorBar(times2keep,a,a_sem,'lineprops','k'); hold on
% shadedErrorBar(times2keep,b,b_sem,'lineprops','r')
% xlim([-2 10])
% title('low gamma:medium gamma ratio')

% figure()
% a = nanmean(time_freq_band.tf_off.low_gamma./time_freq_band.tf_off.theta); hold on;
% b = nanmean(time_freq_band.tf_on.low_gamma./time_freq_band.tf_on.theta); hold on;
% a_sem = nanstd(time_freq_band.tf_off.low_gamma./time_freq_band.tf_off.theta)./sqrt(height(time_freq_band.tf_off.theta));
% b_sem = nanstd(time_freq_band.tf_on.low_gamma./time_freq_band.tf_on.theta)./sqrt(height(time_freq_band.tf_on.theta));
% shadedErrorBar(times2keep,a,a_sem,'lineprops','k'); hold on
% shadedErrorBar(times2keep,b,b_sem,'lineprops','r')
% xlim([-2 10])
% title('low gamma:theta ratio')

%%
inds = find(times2keep>= 0 & times2keep<= 3);
time_freq_band.tf_off_zscored_within_trial.theta;
forward_window_theta = nanmean(time_freq_band.tf_off_zscored_within_trial.theta(:,inds),2);
forward_window_theta(isnan(forward_window_theta)) = [];

inds = find(times2keep> 3 & times2keep<= 10);
time_freq_band.tf_off_zscored_within_trial.theta;
reverse_window_theta = nanmean(time_freq_band.tf_off_zscored_within_trial.theta(:,inds),2);
reverse_window_theta(isnan(reverse_window_theta)) = [];

inds = find(times2keep>= -2 & times2keep<= -1);
time_freq_band.tf_off_zscored_within_trial.theta;
run_window_theta = nanmean(time_freq_band.tf_off_zscored_within_trial.theta(:,inds),2);
run_window_theta(isnan(run_window_theta)) = [];

data = [run_window_theta forward_window_theta reverse_window_theta]
[p, tbl, stats] = friedman(data)
posthocs = multcompare(stats)

[p,h,z] = signrank(run_window_theta,reverse_window_theta)
% [h,p] = ttest(run_window_theta,reverse_window_theta)

figure()
data = [mean(run_window_theta) mean(forward_window_theta) mean(reverse_window_theta)];
err = [std(run_window_theta)/sqrt(length(run_window_theta)) std(forward_window_theta)/sqrt(length(forward_window_theta)) std(reverse_window_theta)/sqrt(length(reverse_window_theta))]

figure('Position',[1986 1051 100 100])
colors = [.4660 0.6740 0.1880; 0.4940 0.1840 0.5560;];



b1 = bar(1,data(1)); hold on;
e1 = errorbar(1,data(1),err(1),'k','linestyle','none');
e1.CapSize = 4;
hold on
b2 = bar(2,data(2)); hold on;
e2 = errorbar(2,data(2),err(2),'k','linestyle','none');
e2.CapSize = 4;
b3 = bar(3,data(3)); hold on;
e3 = errorbar(3,data(3),err(3),'k','linestyle','none');
e3.CapSize = 4;

b1.FaceColor = [0 0 0];
b1.EdgeColor = 'none';
b2.FaceColor = colors(1,:);
b2.EdgeColor = 'none';
b3.FaceColor = colors(2,:);
b3.EdgeColor = 'none';

figure('Position',[100 100 100 100])
data = [run_window_theta; forward_window_theta; reverse_window_theta];
data_labels = [ones(size(run_window_theta)); 2*ones(size(run_window_theta)); 3*ones(size(run_window_theta))]

boxplot(data,data_labels)


set(gcf, 'Color', 'white','Renderer','painters', 'PaperPositionMode', 'auto');
fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures';
saveas(gcf,fullfile(fig_path,'theta_power_in_different_epochs'),'pdf')

