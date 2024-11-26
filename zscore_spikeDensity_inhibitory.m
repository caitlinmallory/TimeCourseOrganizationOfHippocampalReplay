function zscore_spikeDensity_inhibitory(smoothing_sigma)

%clearvars -except dayFiles day directory rat windows hand_clustered_only
plot_figure = 0;
zscored_sd = [];

% Look at posterior quality in a time window surrounding an event (ripple
% peak or laser onset)
load Analysis_Information;
load Experiment_Information;

ripple_tetrode = Experiment_Information.ripple_refs;
lfp_file_denoised = ['LFP_Data' num2str(ripple_tetrode) '_denoised.mat'];
lfp_file = ['LFP_Data' num2str(ripple_tetrode) '.mat'];

load(lfp_file);
load spikeDensity_inhibitory.mat;
load Position_Data;


Run_Times = Experiment_Information.Run_Times;
Sleep_Times = Experiment_Information.Sleep_Times;
unique_Session_Times = vertcat(Experiment_Information.Segments.Times);
Times_day = [min(unique_Session_Times) max(unique_Session_Times)];

%positions
times = load_timeBins_cm(Times_day,0.005*spikeSampRate,0.03*spikeSampRate);
Position_Data = load_positions_full(Run_Times,Sleep_Times,times,Position_Data,4);
Position_Data(isnan(Position_Data(:,5)),5) = 0;

% Smooth spike density
%smoothing_sigma = 0.0125; % desired standard deviation of the Gaussian used for smoothing
%smoothing_sigma = 0.0015; % desired standard deviation of the Gaussian used for smoothing

stepSize = mean(diff(spikeDensity(:,1)))/spikeSampRate;
w = setUp_gaussFilt_sigma(smoothing_sigma,stepSize);
spikeDensity_smoothed = conv(spikeDensity(:,2),w,'same');
spikeDensity = [spikeDensity(:,1),spikeDensity_smoothed];

spike_events = [];


for i = 1:size(unique_Session_Times,1)

    % Do the following separately for each session:
    %load spike density within session and zscore
    LFP_sub = compute_dataTemporalConcatenation(LFP_Data,unique_Session_Times(i,:));
    % Spike Density Events:
    
    spikeDensity_sub = compute_dataTemporalConcatenation(spikeDensity,unique_Session_Times(i,:));
    spikeDensity_sub = compute_dataInterpolation(spikeDensity_sub,LFP_sub(:,1),[]);

    % position:

    Position_Data_sub = compute_dataTemporalConcatenation(Position_Data,unique_Session_Times(i,:));
    Position_Data_sub = compute_dataInterpolation(Position_Data_sub, spikeDensity_sub(:,1), []);    
    moving_inds = find(Position_Data_sub(:,5) > speedThr);
    good_inds = setdiff([1:length(spikeDensity_sub)],moving_inds)'; 
    
    spike_density_mean = nanmean(spikeDensity_sub(good_inds,2));
    spike_density_std = nanstd(spikeDensity_sub(good_inds,2));

    spikeDensity_sub(:,2) = (spikeDensity_sub(:,2)-spike_density_mean)./spike_density_std;
    spikeDensity_sub(:,3) = spikeDensity_sub(:,2);
    spikeDensity_sub(moving_inds,3) = nan;

    
    zscored_sd = [zscored_sd; spikeDensity_sub];

    [spike_density_events_sub, spike_density_event_amplitude_sub] = find_candidate_events_2(spikeDensity_sub(:,3),spikeDensity_sub(:,1));

    % To plot:
    %plot_candidate_events(spike_density_events_sub, spikeDensity_sub, Position_Data_sub, spikeSampRate)
    
 
    spike_events = [spike_events; spike_density_events_sub];

    
end



% if exist('candidateEvents.mat')==2
%     load candidateEvents.mat
% else
%     candidateEvents = struct();
% end
% candidateEvents.ripple_events = ripple_events;
% candidateEvents.spike_events = spike_events;
% if save_laser_events
%     candidateEvents.laser_events = laser_events;
% end
% save('candidateEvents','candidateEvents')

if smoothing_sigma == 0.0125 % defuault
 
  save('zscored_sd_inhibitory','zscored_sd');
else
 save(['zscored_sd_inhibitory_smoothing_' num2str(smoothing_sigma*10000)],'zscored_sd');
end

%  figure()
% plot(Position_Data(:,1)./30000-Position_Data(1,1)/30000,Position_Data(:,5)./10,'linewidth',2);
% hold on
% plot(zscored_ripple_power(:,1)./30000-Position_Data(1,1)/30000,zscored_ripple_power(:,6))
% legend({'0.1*speed'; 'z-scored ripple power'})
% xlabel('time, s')