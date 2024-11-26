

clearvars -except dayFiles day directory rat windows hand_clustered_only
plot_figure = 0;
zscored_ripple_power = [];
zscored_sd = [];


% Look at posterior quality in a time window surrounding an event (ripple
% peak or laser onset)
load Analysis_Information;
load Experiment_Information;
if exist('candidateEvents.mat')==2
    load('candidateEvents.mat');
end

ripple_tetrode = Experiment_Information.ripple_refs;
lfp_file_denoised = ['LFP_Data' num2str(ripple_tetrode) '_denoised.mat'];
lfp_file = ['LFP_Data' num2str(ripple_tetrode) '.mat'];


load(lfp_file_denoised);
load(lfp_file);
load spikeDensity_pyr;
load Position_Data;

save_laser_events = 1;
time_window_post_laser_onset = 0.400; % time after the laser onset that you wish to examine (s)
time_window_pre_laser_onset = 0.400; % time before the laser onset that you wish to examine (s)

Run_Times = Experiment_Information.Run_Times;
Sleep_Times = Experiment_Information.Sleep_Times;


unique_Session_Times = vertcat(Experiment_Information.Segments.Times);

Times_day = [min(unique_Session_Times) max(unique_Session_Times)];

%positions
times = load_timeBins_cm(Times_day,0.005*spikeSampRate,0.02*spikeSampRate);
Position_Data = load_positions_full(Run_Times,Sleep_Times,times,Position_Data,4);
Position_Data(isnan(Position_Data(:,5)),5) = 0;


% filter LFP for ripples and smooth
LFP = LFP_Data;
clear LFP_data
[~,LFP_filt_SWRAmp,LFP_filt_SWR] = compute_filteredLFP_cm(SWRFreqRange,LFP(:,2),lfpSampRate);
LFP = [LFP,LFP_filt_SWR,LFP_filt_SWRAmp];
% Column 1 = time
% Column 2 = raw lfp
% Column 3 = SWR filtered lfp
% Column 4 = SWR amplitude
% Column 5 = SWR amplitude, smoothed 
% Column 6 = SWR zscored



%smooth ripple power:
smoothing_sigma = 0.0125; % desired standard deviation of the Gaussian used for smoothing
stepSize = mean(diff(LFP(:,1)))/spikeSampRate;
w = setUp_gaussFilt_sigma(smoothing_sigma,stepSize);
ripplePower_smoothed = conv(LFP(:,4),w,'same');
LFP(:,5) = ripplePower_smoothed;

% Smooth spike density
smoothing_sigma = 0.0125; % desired standard deviation of the Gaussian used for smoothing
stepSize = mean(diff(spikeDensity(:,1)))/spikeSampRate;
w = setUp_gaussFilt_sigma(smoothing_sigma,stepSize);
spikeDensity_smoothed = conv(spikeDensity(:,2),w,'same');
spikeDensity = [spikeDensity(:,1),spikeDensity_smoothed];

ripple_events = [];
spike_events = [];

if exist('LaserTimestampData.mat')
    load LaserTimestampData
 
    laser_on_timestamps = laser_timestamps_cat(laser_timestamps_cat(:,2)==1,1);
    laser_off_timestamps = laser_timestamps_cat(laser_timestamps_cat(:,2)==0,1);
    laser_on_duration = (laser_off_timestamps-laser_on_timestamps)./spikeSampRate;
    laser_events = [laser_on_timestamps - time_window_pre_laser_onset*spikeSampRate laser_on_timestamps + time_window_post_laser_onset*spikeSampRate];
else
    laser_events = [];
end


for i = 1:size(unique_Session_Times,1)

    LFP_sub = compute_dataTemporalConcatenation(LFP,unique_Session_Times(i,:));
    LFP_Data_denoised_sub = compute_dataTemporalConcatenation(LFP_Data_denoised,unique_Session_Times(i,:));
    
    % Do the following separately for each session:
    %load LFP ripple amp within session and zscore
    
    % position:
    Position_Data_sub = compute_dataTemporalConcatenation(Position_Data,unique_Session_Times(i,:));
    Position_Data_sub = compute_dataInterpolation(Position_Data_sub, LFP_sub(:,1), []);
    
    %limit to stationary periods and mask out artifacts prior to z-scoring:
    artifact_inds = find(isnan(LFP_Data_denoised_sub(:,2)));
    moving_inds = find(Position_Data_sub(:,5) > speedThr);
    bad_inds = unique([artifact_inds; moving_inds]);

    % save the good inds (stationary) that will be used for computing mean
    % and variance. 
    good_inds = setdiff(1:length(LFP_sub),bad_inds)'; 

    % compute mean and variance of the smoothed ripple amplitude, stationary periods only:
    mean_ripple_power = nanmean(LFP_sub(good_inds,5));
    std_ripple_power = nanstd(LFP_sub(good_inds,5));

    LFP_sub(:,6) = (LFP_sub(:,5)-mean_ripple_power)./std_ripple_power;
    LFP_sub(:,7) = LFP_sub(:,6);
    LFP_sub(bad_inds,7) = nan;

    %save a copy of the z-scored ripple power, without removing any data
    %points.

    zscored_ripple_power = [zscored_ripple_power; LFP_sub];
    
    [ripple_events_sub, ripple_event_amplitude_sub] = find_candidate_events_2(LFP_sub(:,7),LFP_sub(:,1));

    % To plot:
    % plot_candidate_events(ripple_events_sub, LFP_sub, Position_Data_sub, spikeSampRate)
    
    % Spike Density Events:
    
    spikeDensity_sub = compute_dataTemporalConcatenation(spikeDensity,unique_Session_Times(i,:));
    spikeDensity_sub = compute_dataInterpolation(spikeDensity_sub,Position_Data_sub(:,1),[]);
    
    
    good_inds = setdiff(1:length(spikeDensity_sub),moving_inds)';
    spike_density_mean = nanmean(spikeDensity_sub(good_inds,2));
    spike_density_std = nanstd(spikeDensity_sub(good_inds,2));

    spikeDensity_sub(:,2) = (spikeDensity_sub(:,2)-spike_density_mean)./spike_density_std;
    spikeDensity_sub(:,3) = spikeDensity_sub(:,2);
    spikeDensity_sub(moving_inds,3) = nan;

    
    zscored_sd = [zscored_sd; spikeDensity_sub];

    [spike_density_events_sub, spike_density_event_amplitude_sub] = find_candidate_events_2(spikeDensity_sub(:,3),spikeDensity_sub(:,1));

    % To plot:
    %plot_candidate_events(spike_density_events_sub, spikeDensity_sub, Position_Data_sub, spikeSampRate)
    
    ripple_events = [ripple_events; ripple_events_sub];
    spike_events = [spike_events; spike_density_events_sub];

    
end

if exist('candidateEvents.mat')==2
    load candidateEvents.mat
else
    candidateEvents = struct();
end
candidateEvents.ripple_events = ripple_events;
candidateEvents.spike_events = spike_events;
if save_laser_events
    candidateEvents.laser_events = laser_events;
end
save('candidateEvents','candidateEvents')

save('zscored_ripple_power','zscored_ripple_power');
save('zscored_sd','zscored_sd');


%  figure()
% plot(Position_Data(:,1)./30000-Position_Data(1,1)/30000,Position_Data(:,5)./10,'linewidth',2);
% hold on
% plot(zscored_ripple_power(:,1)./30000-Position_Data(1,1)/30000,zscored_ripple_power(:,6))
% legend({'0.1*speed'; 'z-scored ripple power'})
% xlabel('time, s')