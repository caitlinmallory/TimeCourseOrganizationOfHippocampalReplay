clearvars -except dayFiles day directory rat windows hand_clustered_only

plot_fig=0;
% Logic of the single cell analysis:

% Find times when the rat stopped on the left hand of the track.
% Look at cells that only have a place field in the left running direction
% compared to cells that only have a place field in the right running
% direction.

% Prediction is that cells with left-run fields will be suppressed for ~3
% seconds, and then fire after 3 seconds. Cells with right-fields will
% firing right away and then ramp down. Cells with BOTH right and left
% fields will initially be suppressed, and then firing strongly.

% Repeat this for stops on the right hand of the track (and combine).

% First, load stops. 

load Experiment_Information.mat
load clusters

if Experiment_Information.spatialDim == 2
load true_drink_periods.mat
else
load Behavior_Data.mat
end

load Position_Data;

bin_width =0.025;
spikeDensityStepSize = 0.025;
smoothing_sigma = 0.25;
% compute the instantaneous firing rate of each cell over the session

spikeSampRate = Experiment_Information.spikeSampRate;
Run_Times = Experiment_Information.Run_Times;
Sleep_Times = Experiment_Information.Sleep_Times;

times_list = cat(2, (cat(2, Run_Times{:})), Sleep_Times{:});
Times_day = [min(times_list) max(times_list)];

%times to estimate spike density
timeEdges = (Times_day(1):spikeDensityStepSize*spikeSampRate:Times_day(2))';

%get speed at the time edges
speed = compute_dataTemporalConcatenation(Position_Data,[Times_day(1) Times_day(2)]);
speed = compute_dataInterpolation(speed,timeEdges,[]);
speed = speed(:,5);

spikeDensity = zeros(length(clusters),length(timeEdges));


left_field = zeros(length(clusters),1);
right_field = zeros(length(clusters),1);
left_right_max_spatial_firing_rate_ratio = nan(length(clusters),1);
left_right_mean_spatial_firing_rate_ratio = nan(length(clusters),1);
left_mean_infield_firing_rate = nan(length(clusters),1);
right_mean_infield_firing_rate = nan(length(clusters),1);
left_max_spatial_firing_rate = nan(length(clusters),1);
right_max_spatial_firing_rate = nan(length(clusters),1);
left_com_spatial_firing_bin_location = nan(length(clusters),1);
right_com_spatial_firing_bin_location = nan(length(clusters),1);
left_com_spatial_firing_bin_location_pcnt = nan(length(clusters),1);
right_com_spatial_firing_bin_location_pcnt = nan(length(clusters),1);

for i = 1:length(clusters)
    spkTimes = clusters(i).spkTime;
    spkTimes(find(diff(spkTimes)==0),:) = [];
    if ~isempty(spkTimes)
        spikeDensity(i,:) = histc(spkTimes,timeEdges);
    end
    if clusters(i).runs.directions(1).num_subFields>0
    left_field(i)=1;
    end
    if clusters(i).runs.directions(2).num_subFields>0
    right_field(i)=1;
    end
    left_right_max_spatial_firing_rate_ratio(i) = clusters(i).runs.directions(1).maxSpatialFiringRate/clusters(i).runs.directions(2).maxSpatialFiringRate;
    left_right_mean_spatial_firing_rate_ratio(i) = clusters(i).runs.directions(1).meanSpatialFiringRate/clusters(i).runs.directions(2).meanSpatialFiringRate;

    left_mean_infield_firing_rate(i) = clusters(i).runs.directions(1).Mean_InField_Firing_Rate;
    right_mean_infield_firing_rate(i) = clusters(i).runs.directions(2).Mean_InField_Firing_Rate;
    left_max_spatial_firing_rate(i) = clusters(i).runs.directions(1).maxSpatialFiringRate;
    right_max_spatial_firing_rate(i) = clusters(i).runs.directions(2).maxSpatialFiringRate;
    left_com_spatial_firing_bin_location(i) = clusters(i).runs.directions(1).COMSpatialFiringBinLoc;
    right_com_spatial_firing_bin_location(i) = clusters(i).runs.directions(2).COMSpatialFiringBinLoc;
    left_com_spatial_firing_bin_location_pcnt(i) = left_com_spatial_firing_bin_location(i)/(Experiment_Information.maze_size/2);
    right_com_spatial_firing_bin_location_pcnt(i) = right_com_spatial_firing_bin_location(i)/(Experiment_Information.maze_size/2);
end




clusters_left_field_only = (left_field==1 & right_field==0);
clusters_right_field_only = (left_field==0 & right_field==1);
clusters_right_and_left_field = (left_field==1 & right_field==1);

clusters_left_max_rate_higher = left_right_max_spatial_firing_rate_ratio>1;
clusters_left_mean_rate_higher = left_right_mean_spatial_firing_rate_ratio>1;
clusters_right_max_rate_higher = left_right_max_spatial_firing_rate_ratio<1;
clusters_right_mean_rate_higher = left_right_mean_spatial_firing_rate_ratio<1;
max_spatial_firing_rate_ratio = left_right_max_spatial_firing_rate_ratio;
max_spatial_firing_rate_ratio(max_spatial_firing_rate_ratio<1)=1./(max_spatial_firing_rate_ratio(max_spatial_firing_rate_ratio<1));
mean_spatial_firing_rate_ratio = left_right_mean_spatial_firing_rate_ratio;
mean_spatial_firing_rate_ratio(mean_spatial_firing_rate_ratio<1)=1./(mean_spatial_firing_rate_ratio(mean_spatial_firing_rate_ratio<1));

spikeDensity = spikeDensity./spikeDensityStepSize;
w = setUp_gaussFilt_sigma(smoothing_sigma,spikeDensityStepSize);
spikeDensity_smoothed = zeros(size(spikeDensity));
for i = 1:size(spikeDensity,1)
spikeDensity_smoothed(i,:) = conv(spikeDensity(i,:),w,'same');
end

spikeDensity = spikeDensity_smoothed;

spikeDensity_stationary = spikeDensity;
spikeDensity_stationary(:,speed>5)=nan; % remove any periods from spike density where rat was moving.

%% Limit spike density to either replay or SDEs if desired
%% past or future representations, using SDEs as candidate events:
% Start with SDEs. Then look for representations of the past or future map
replay_weighted_r_thr=0.6;
replay_coverage_thr=0.2;
replay_posterior_diff_thr = 0.33;

events = 'spike_filtered';
sde_thr = 3;
ripple_thr = -inf;
use_duration_og = 0;
posterior_diff_thr = replay_posterior_diff_thr;
coverage_thr = replay_coverage_thr;
weighted_r_thr = replay_weighted_r_thr;
min_time_into_stopping_period=0;
max_time_into_stopping_period = inf;

sub_events = load_replays_from_individual_session(events,use_duration_og,coverage_thr,...
    weighted_r_thr,posterior_diff_thr,sde_thr,ripple_thr,min_time_into_stopping_period,max_time_into_stopping_period);
replays = sub_events(sub_events.replay==1,:);
replay_times=replays.timePoints_og;

spikeDensity_during_replay = nan(size(spikeDensity));
for event=1:height(replays)
    [~,ind_start] = min(abs(timeEdges-replay_times(event,1)));
    [~,ind_stop] = min(abs(timeEdges-replay_times(event,2)));
    spikeDensity_during_replay(:,ind_start:ind_stop) = spikeDensity(:,ind_start:ind_stop);
end

events = 'spike_filtered';
sde_thr = 3;
ripple_thr = -inf;
use_duration_og = 1;
posterior_diff_thr = 0;
coverage_thr = -inf;
weighted_r_thr = -inf;
min_time_into_stopping_period = 0;
max_time_into_stopping_period = 10;

sdes = load_replays_from_individual_session(events,use_duration_og,coverage_thr,...
    weighted_r_thr,posterior_diff_thr,sde_thr,ripple_thr,min_time_into_stopping_period,max_time_into_stopping_period);
sde_times=sdes.timePoints_og;

spikeDensity_during_sdes = nan(size(spikeDensity));
for event=1:height(sdes)
    [~,ind_start] = min(abs(timeEdges-sde_times(event,1)));
    [~,ind_stop] = min(abs(timeEdges-sde_times(event,2)));
    spikeDensity_during_sdes(:,ind_start:ind_stop) = spikeDensity(:,ind_start:ind_stop);
end
%%

stopping_period_times = Reward_Epoch_Time_Boundaries_speed_thresholded;
end_zone = Reward_Epoch_Time_Boundaries_endzone;


%%
length_to_plot = 10;
bins_forward = length_to_plot/spikeDensityStepSize;

firing_rates = nan(length(clusters),2*bins_forward,length(stopping_period_times));
firing_rates_replay = nan(length(clusters),2*bins_forward,length(stopping_period_times));
firing_rates_sde = nan(length(clusters),2*bins_forward,length(stopping_period_times));
cluster_peak_firing_rates_pre=nan(length(clusters),length(stopping_period_times));
cluster_mean_firing_rates_pre =nan(length(clusters),length(stopping_period_times));
cluster_mean_firing_rates_replay = nan(length(clusters),length(stopping_period_times));
cluster_mean_firing_rates_after_stopping = nan(length(clusters),length(stopping_period_times));

for i = 2:length(stopping_period_times)
    run_times = Run_Epoch_Time_Boundaries;
    run_times(run_times>stopping_period_times(i,1))= nan;
    [~,ind] = min(abs(run_times(:,1)-stopping_period_times(i,1)));
    start_time_of_previous_run = run_times(ind,1);

    stopping_period_start_time = stopping_period_times(i,1);
    if isnan(stopping_period_start_time)
        continue
    end
    [~,stopping_period_start_ind] = min(abs(timeEdges-stopping_period_start_time));
    [~,stopping_period_start_3s_ind] = min(abs(timeEdges-(stopping_period_start_time+3.*30000)));
    [~,ind_of_start_of_previous_run] = min(abs(timeEdges-start_time_of_previous_run));

    % look back either 10 sec or one lap
    firing_rate_pre_stop = spikeDensity(:,ind_of_start_of_previous_run:stopping_period_start_ind);
    peak_firing_rate_in_pre_lap = max(firing_rate_pre_stop,[],2);
    cluster_peak_firing_rates_pre(:,i) = peak_firing_rate_in_pre_lap;
    cluster_mean_firing_rates_pre(:,i) = nanmean(firing_rate_pre_stop,2);
    pre_period_start_ind = max(ind_of_start_of_previous_run,stopping_period_start_ind-bins_forward);
    num_bins_backward = stopping_period_start_ind-pre_period_start_ind;
  
    % look at the first 3 sec of replay, post-stop
    length_of_replay_post_stop = 3;
    firing_rate_post_stop = spikeDensity_during_replay(:,stopping_period_start_ind:stopping_period_start_3s_ind);
    mean_firing_rate_in_first_three_sec_replay = nanmean(firing_rate_post_stop,2);
    cluster_mean_firing_rates_replay(:,i) = mean_firing_rate_in_first_three_sec_replay;
    
    %look at the first 3 of data after stopping (not just replay)
    firing_rate_post_stop = spikeDensity_stationary(:,stopping_period_start_ind:stopping_period_start_3s_ind);
    mean_firing_rate_in_first_three_sec = nanmean(firing_rate_post_stop,2);
    cluster_mean_firing_rates_after_stopping(:,i) = mean_firing_rate_in_first_three_sec;
    

[stopping_period_end_time]= stopping_period_times(i,2);
[~,stopping_period_end_ind] = min(abs(timeEdges-stopping_period_end_time));
inds_to_plot_forward = stopping_period_start_ind:min((stopping_period_start_ind + (bins_forward-1)),stopping_period_end_ind); 
num_bins_forward = length(inds_to_plot_forward);
inds_to_plot_backward = pre_period_start_ind:(stopping_period_start_ind-1);

if any(inds_to_plot_forward<0)
    continue 
end
if any(inds_to_plot_backward<0)
    continue
end
if any(inds_to_plot_forward>length(spikeDensity))
    continue
end

sub_rates = nan(length(clusters),2*bins_forward);
sub_rates(:,(bins_forward+1):(bins_forward+1)+length(inds_to_plot_forward)-1) = spikeDensity(:,inds_to_plot_forward);
sub_rates(:,(bins_forward+1)-num_bins_backward:bins_forward) = spikeDensity(:,inds_to_plot_backward);

sub_rates_replay = nan(length(clusters),2*bins_forward);
sub_rates_replay(:,(bins_forward+1):(bins_forward+1)+length(inds_to_plot_forward)-1) = spikeDensity_during_replay(:,inds_to_plot_forward);
sub_rates_replay(:,(bins_forward+1)-num_bins_backward:bins_forward) = spikeDensity(:,inds_to_plot_backward);

sub_rates_sdes = nan(length(clusters),2*bins_forward);
sub_rates_sdes(:,(bins_forward+1):(bins_forward+1)+length(inds_to_plot_forward)-1) = spikeDensity_during_sdes(:,inds_to_plot_forward);
sub_rates_sdes(:,(bins_forward+1)-num_bins_backward:bins_forward) = spikeDensity(:,inds_to_plot_backward);

firing_rates(:,:,i) =sub_rates ;
firing_rates_replay(:,:,i) = sub_rates_replay;
firing_rates_sde(:,:,i) = sub_rates_sdes;
end

firing_rates([clusters.Excitatory]==0,:,:)=nan;
firing_rates_replay([clusters.Excitatory]==0,:,:)=nan;
firing_rates_sde([clusters.Excitatory]==0,:,:)=nan;

mean_frs = nan(length(clusters),1);
std_frs = nan(length(clusters),1);
mean_frs_replay = nan(length(clusters),1);
std_frs_replay = nan(length(clusters),1);
mean_frs_sde = nan(length(clusters),1);
std_frs_sde = nan(length(clusters),1);


for i=1:length(clusters)
    sub = squeeze(firing_rates(i,:,:));
    mean_frs(i) = nanmean(sub(:));
    std_frs(i) = nanstd(sub(:));

    sub = squeeze(firing_rates_replay(i,:,:));
    mean_frs_replay(i) = nanmean(sub(:));
    std_frs_replay(i) = nanstd(sub(:));

    sub = squeeze(firing_rates_sde(i,:,:));
    mean_frs_sde(i) = nanmean(sub(:));
    std_frs_sde(i) = nanstd(sub(:));
end

peak_fr_left_laps_pre= cluster_peak_firing_rates_pre(:,end_zone==1);
peak_fr_left_laps_pre_mean = nanmean(peak_fr_left_laps_pre,2);
peak_fr_left_laps_pre_std = nanstd(peak_fr_left_laps_pre,[],2);
peak_fr_right_laps_pre = cluster_peak_firing_rates_pre(:,end_zone==2);
peak_fr_right_laps_pre_mean = nanmean(peak_fr_right_laps_pre,2);
peak_fr_right_laps_pre_std = nanstd(peak_fr_right_laps_pre,[],2);

ave_fr_left_laps_pre= cluster_mean_firing_rates_pre(:,end_zone==1);
ave_fr_left_laps_pre_mean = nanmean(ave_fr_left_laps_pre,2);
ave_fr_left_laps_pre_std = nanstd(ave_fr_left_laps_pre,[],2);
ave_fr_right_laps_pre = cluster_mean_firing_rates_pre(:,end_zone==2);
ave_fr_right_laps_pre_mean = nanmean(ave_fr_right_laps_pre,2);
ave_fr_right_laps_pre_std = nanstd(ave_fr_right_laps_pre,[],2);


ave_fr_in_replay_left_laps_post= cluster_mean_firing_rates_replay(:,end_zone==1);
ave_fr_in_replay_left_laps_post_mean = nanmean(cluster_mean_firing_rates_replay,2);
ave_fr_in_replay_left_laps_post_std = nanstd(cluster_mean_firing_rates_replay,[],2);
ave_fr_in_replay_right_laps_post = cluster_mean_firing_rates_replay(:,end_zone==2);
ave_fr_in_replay_right_laps_post_mean = nanmean(cluster_mean_firing_rates_replay,2);
ave_fr_in_replay_right_laps_post_std = nanstd(cluster_mean_firing_rates_replay,[],2);

ave_fr_left_laps_post= cluster_mean_firing_rates_after_stopping(:,end_zone==1);
ave_fr_left_laps_post_mean = nanmean(cluster_mean_firing_rates_after_stopping,2);
ave_fr_left_laps_post_std = nanstd(cluster_mean_firing_rates_after_stopping,[],2);
ave_fr_right_laps_post = cluster_mean_firing_rates_after_stopping(:,end_zone==2);
ave_fr_right_laps_post_mean = nanmean(cluster_mean_firing_rates_after_stopping,2);
ave_fr_right_laps_post_std = nanstd(cluster_mean_firing_rates_after_stopping,[],2);

zscored_peak_pre_rates = nan(length(clusters),length(stopping_period_times));
zscored_peak_pre_rates(:,end_zone==1) = (peak_fr_left_laps_pre-(peak_fr_left_laps_pre_mean))./peak_fr_left_laps_pre_std;
zscored_peak_pre_rates(:,end_zone==2) = (peak_fr_right_laps_pre-(peak_fr_right_laps_pre_mean))./peak_fr_right_laps_pre_std;
zscored_peak_pre_rates_left=zscored_peak_pre_rates(:,end_zone==1);
zscored_peak_pre_rates_right=zscored_peak_pre_rates(:,end_zone==2);

zscored_ave_pre_rates = nan(length(clusters),length(stopping_period_times));
zscored_ave_pre_rates(:,end_zone==1) = (ave_fr_left_laps_pre-(ave_fr_left_laps_pre_mean))./ave_fr_left_laps_pre_std;
zscored_ave_pre_rates(:,end_zone==2) = (ave_fr_right_laps_pre-(ave_fr_right_laps_pre_mean))./ave_fr_right_laps_pre_std;
zscored_ave_pre_rates_left=zscored_ave_pre_rates(:,end_zone==1);
zscored_ave_pre_rates_right=zscored_ave_pre_rates(:,end_zone==2);

zscored_ave_post_rates_ = nan(length(clusters),length(stopping_period_times));
zscored_ave_post_rates(:,end_zone==1) = (ave_fr_left_laps_post-(ave_fr_left_laps_post_mean))./ave_fr_left_laps_post_std;
zscored_ave_post_rates(:,end_zone==2) = (ave_fr_right_laps_post-(ave_fr_right_laps_post_mean))./ave_fr_right_laps_post_std;
zscored_ave_post_rates_left=zscored_ave_post_rates(:,end_zone==1);
zscored_ave_post_rates_right=zscored_ave_post_rates(:,end_zone==2);

zscored_ave_post_rates_in_replay = nan(length(clusters),length(stopping_period_times));
zscored_ave_post_rates_in_replay(:,end_zone==1) = (ave_fr_in_replay_left_laps_post-(ave_fr_in_replay_left_laps_post_mean))./ave_fr_in_replay_left_laps_post_std;
zscored_ave_post_rates_in_replay(:,end_zone==2) = (ave_fr_in_replay_right_laps_post-(ave_fr_in_replay_right_laps_post_mean))./ave_fr_in_replay_right_laps_post_std;
zscored_ave_post_rates_left_in_replay=zscored_ave_post_rates_in_replay(:,end_zone==1);
zscored_ave_post_rates_right_in_replay=zscored_ave_post_rates_in_replay(:,end_zone==2);



zscored_firing_rates = (firing_rates-mean_frs)./std_frs;
zscored_firing_rates_replay = (firing_rates_replay-mean_frs_replay)./std_frs_replay;
zscored_firing_rates_sde = (firing_rates_sde-mean_frs_sde)./std_frs_sde;

firing_rates_at_left = squeeze(nanmean(firing_rates(:,:,end_zone==1),3));
firing_rates_at_right = squeeze(nanmean(firing_rates(:,:,end_zone==2),3));
zscored_firing_rates_at_left = squeeze(nanmean(zscored_firing_rates(:,:,end_zone==1),3));
zscored_firing_rates_at_right = squeeze(nanmean(zscored_firing_rates(:,:,end_zone==2),3)); 

firing_rates_at_left_replay = squeeze(nanmean(firing_rates_replay(:,:,end_zone==1),3));
firing_rates_at_right_replay = squeeze(nanmean(firing_rates_replay(:,:,end_zone==2),3));
zscored_firing_rates_at_left_replay = squeeze(nanmean(zscored_firing_rates_replay(:,:,end_zone==1),3));
zscored_firing_rates_at_right_replay = squeeze(nanmean(zscored_firing_rates_replay(:,:,end_zone==2),3)); 

firing_rates_at_left_sde = squeeze(nanmean(firing_rates_sde(:,:,end_zone==1),3));
firing_rates_at_right_sde = squeeze(nanmean(firing_rates_sde(:,:,end_zone==2),3));
zscored_firing_rates_at_left_sde = squeeze(nanmean(zscored_firing_rates_sde(:,:,end_zone==1),3));
zscored_firing_rates_at_right_sde = squeeze(nanmean(zscored_firing_rates_sde(:,:,end_zone==2),3)); 

firing_rates_at_left_low = nan(length(clusters),size(firing_rates,2));
firing_rates_at_left_high = nan(length(clusters),size(firing_rates,2));
firing_rates_at_right_low = nan(length(clusters),size(firing_rates,2));
firing_rates_at_right_high = nan(length(clusters),size(firing_rates,2));
firing_rates_at_left_replay_low = nan(length(clusters),size(firing_rates,2));
firing_rates_at_left_replay_high = nan(length(clusters),size(firing_rates,2));
firing_rates_at_right_replay_low = nan(length(clusters),size(firing_rates,2));
firing_rates_at_right_replay_high = nan(length(clusters),size(firing_rates,2));
firing_rates_at_left_sde_low = nan(length(clusters),size(firing_rates,2));
firing_rates_at_left_sde_high = nan(length(clusters),size(firing_rates,2));
firing_rates_at_right_sde_low = nan(length(clusters),size(firing_rates,2));
firing_rates_at_right_sde_high = nan(length(clusters),size(firing_rates,2));


for i = 1:length(clusters)
    low = (zscored_peak_pre_rates(i,:)<0)';
    high = (zscored_peak_pre_rates(i,:)>0)';

    firing_rates_at_left_low(i,:) = nanmean(squeeze(firing_rates(i,:,end_zone==1&low)),2);
    firing_rates_at_left_high(i,:) = nanmean(squeeze(firing_rates(i,:,end_zone==1&high)),2);
    firing_rates_at_right_low(i,:) = nanmean(squeeze(firing_rates(i,:,end_zone==2&low)),2);
    firing_rates_at_right_high(i,:) = nanmean(squeeze(firing_rates(i,:,end_zone==2&high)),2);    

    firing_rates_at_left_replay_low(i,:) = nanmean(squeeze(firing_rates_replay(i,:,end_zone==1&low)),2);
    firing_rates_at_left_replay_high(i,:) = nanmean(squeeze(firing_rates_replay(i,:,end_zone==1&high)),2);
    firing_rates_at_right_replay_low(i,:) = nanmean(squeeze(firing_rates_replay(i,:,end_zone==2&low)),2);
    firing_rates_at_right_replay_high(i,:) = nanmean(squeeze(firing_rates_replay(i,:,end_zone==2&high)),2);  
    
    firing_rates_at_left_sde_low(i,:) = nanmean(squeeze(firing_rates_sde(i,:,end_zone==1&low)),2);
    firing_rates_at_left_sde_high(i,:) = nanmean(squeeze(firing_rates_sde(i,:,end_zone==1&high)),2);
    firing_rates_at_right_sde_low(i,:) = nanmean(squeeze(firing_rates_sde(i,:,end_zone==2&low)),2);
    firing_rates_at_right_sde_high(i,:) = nanmean(squeeze(firing_rates_sde(i,:,end_zone==2&high)),2);  
end



if plot_fig==1
subplot(1,2,1)
imagesc(firing_rates_at_left_replay./max(max(firing_rates_at_right_replay,[],2),max(firing_rates_at_left_replay,[],2)))
xticks(0.5:5:50.5)
xticklabels(0:1:10)
subplot(1,2,2)
imagesc(firing_rates_at_right./max(max(firing_rates_at_right,[],2),max(firing_rates_at_left,[],2)))
xticks(0.5:5:50.5)
xticklabels(0:1:10)

subplot(2,2,1)
norm_left_rates = firing_rates_at_left./max(max(firing_rates_at_right,[],2),max(firing_rates_at_left,[],2));
imagesc(norm_left_rates(left_field==1,:))
subplot(2,2,2)
imagesc(norm_left_rates(right_field==1,:))


subplot(1,3,1)
norm_left_rates = firing_rates_at_left./max(max(firing_rates_at_right,[],2),max(firing_rates_at_left,[],2));
imagesc(norm_left_rates(clusters_left_field_only==1,:))
subplot(1,3,2)
imagesc(norm_left_rates(clusters_right_field_only==1,:))
subplot(1,3,3)
imagesc(norm_left_rates(clusters_right_and_left_field==1,:))

subplot(1,3,1)
imagesc(firing_rates_at_left(clusters_left_field_only==1,:))
subplot(1,3,2)
imagesc(firing_rates_at_left(clusters_right_field_only==1,:))
subplot(1,3,3)
imagesc(firing_rates_at_left(clusters_right_and_left_field==1,:))




% a = clusters_left_field_only;
% b = clusters_right_field_only;
% c = clusters_right_and_left_field;
a = max_spatial_firing_rate_ratio>5 & clusters_left_max_rate_higher;
b = max_spatial_firing_rate_ratio>5 & clusters_right_max_rate_higher;
c = setdiff(1:length(clusters),unique([find([clusters.Excitatory]==0)';find(a==1);find(b==1)]));
figure()
subplot(2,2,1)
% plot cells that were just active
plot(nanmean(firing_rates_at_left(a,:)))
hold on
% plot cells that are usually active in other direction
plot(nanmean(firing_rates_at_left(b,:)))
hold on
% plot cells that are usually active in both
plot(nanmean(firing_rates_at_left(c,:)))
subplot(2,2,2)
norm_left_rates = firing_rates_at_left./max(firing_rates_at_left,[],2);
%norm_left_rates = firing_rates_at_left./max(max(firing_rates_at_left,[],2),max(firing_rates_at_right,[],2));
plot(nanmean(norm_left_rates(a,:)))
hold on
plot(nanmean(norm_left_rates(b,:)))
hold on
plot(nanmean(norm_left_rates(c,:)))

subplot(2,2,3)
plot(nanmean(firing_rates_at_right(b,:)))
hold on
plot(nanmean(firing_rates_at_right(a,:)))
hold on
plot(nanmean(firing_rates_at_right(c,:)))

subplot(2,2,4)
norm_right_rates = firing_rates_at_right./max(firing_rates_at_right,[],2);
%norm_right_rates = firing_rates_at_right./max(max(firing_rates_at_right,[],2),max(firing_rates_at_left,[],2));
plot(nanmean(norm_right_rates(b,:)))
hold on
plot(nanmean(norm_right_rates(a,:)))
hold on
plot(nanmean(norm_right_rates(c,:)))
end

cluster_rate_table = table();
cluster_rate_table.Excitatory=[clusters.Excitatory]';
cluster_rate_table.Well_Isolated = [clusters.well_isolated]';
cluster_rate_table.lr_peak_firing_rate_difference = [clusters.lr_peak_firing_rate_difference]';
cluster_rate_table.lr_firing_rate_modulation = [clusters.lr_firing_rate_modulation]';
cluster_rate_table.lr_map_directional_correlation = [clusters.lr_map_directional_correlation]';
cluster_rate_table.left_mean_infield_firing_rate = left_mean_infield_firing_rate;
cluster_rate_table.right_mean_infield_firing_rate = right_mean_infield_firing_rate;
cluster_rate_table.left_max_spatial_firing_rate = left_max_spatial_firing_rate;
cluster_rate_table.right_max_spatial_firing_rate = right_max_spatial_firing_rate;
cluster_rate_table.left_com_spatial_firing_bin_location = left_com_spatial_firing_bin_location;
cluster_rate_table.right_com_spatial_firing_bin_location = right_com_spatial_firing_bin_location;
cluster_rate_table.left_com_spatial_firing_bin_location = left_com_spatial_firing_bin_location_pcnt;
cluster_rate_table.right_com_spatial_firing_bin_location = right_com_spatial_firing_bin_location_pcnt;


cluster_rate_table.clusters_left_max_rate_higher = clusters_left_max_rate_higher;
cluster_rate_table.clusters_right_max_rate_higher = clusters_right_max_rate_higher;
cluster_rate_table.max_spatial_firing_rate_ratio = max_spatial_firing_rate_ratio;
cluster_rate_table.mean_spatial_firing_rate_ratio = mean_spatial_firing_rate_ratio;
cluster_rate_table.left_field = left_field;
cluster_rate_table.right_field = right_field;
cluster_rate_table.right_field = right_field;
cluster_rate_table.max_spatial_firing_rate_ratio = max_spatial_firing_rate_ratio;
cluster_rate_table.mean_spatial_firing_rate_ratio = mean_spatial_firing_rate_ratio;

cluster_rate_table.firing_rates_at_left = firing_rates_at_left;
cluster_rate_table.firing_rates_at_right = firing_rates_at_right;
cluster_rate_table.zscored_firing_rates_at_left = zscored_firing_rates_at_left;
cluster_rate_table.zscored_firing_rates_at_right = zscored_firing_rates_at_right;

cluster_rate_table.firing_rates_at_left_replay = firing_rates_at_left_replay;
cluster_rate_table.firing_rates_at_right_replay = firing_rates_at_right_replay;
cluster_rate_table.zscored_firing_rates_at_left_replay = zscored_firing_rates_at_left_replay;
cluster_rate_table.zscored_firing_rates_at_right_replay = zscored_firing_rates_at_right_replay;

cluster_rate_table.firing_rates_at_left_sde = firing_rates_at_left_sde;
cluster_rate_table.firing_rates_at_right_sde = firing_rates_at_right_sde;
cluster_rate_table.zscored_firing_rates_at_left_sde = zscored_firing_rates_at_left_sde;
cluster_rate_table.zscored_firing_rates_at_right_sde = zscored_firing_rates_at_right_sde;

cluster_rate_table.clusters_left_field_only = clusters_left_field_only ;
cluster_rate_table.clusters_right_field_only = clusters_right_field_only;
cluster_rate_table.clusters_right_and_left_field = clusters_right_and_left_field;
cluster_rate_table.clusters_left_max_rate_higher  = clusters_left_max_rate_higher;
cluster_rate_table.clusters_left_mean_rate_higher = clusters_left_mean_rate_higher;
cluster_rate_table.clusters_right_max_rate_higher = clusters_right_max_rate_higher;
cluster_rate_table.clusters_right_mean_rate_higher = clusters_right_mean_rate_higher;
cluster_rate_table.max_spatial_firing_rate_ratio= max_spatial_firing_rate_ratio;
cluster_rate_table.mean_spatial_firing_rate_ratio= mean_spatial_firing_rate_ratio;

cluster_rate_table.firing_rates_at_left_low = firing_rates_at_left_low;
cluster_rate_table.firing_rates_at_left_high = firing_rates_at_left_high;
cluster_rate_table.firing_rates_at_right_low = firing_rates_at_right_low;
cluster_rate_table.firing_rates_at_right_high = firing_rates_at_right_high;

cluster_rate_table.firing_rates_at_left_replay_low = firing_rates_at_left_replay_low;
cluster_rate_table.firing_rates_at_left_replay_high = firing_rates_at_left_replay_high;
cluster_rate_table.firing_rates_at_right_replay_low = firing_rates_at_right_replay_low;
cluster_rate_table.firing_rates_at_right_replay_high = firing_rates_at_right_replay_high;

cluster_rate_table.firing_rates_at_left_sde_low = firing_rates_at_left_sde_low;
cluster_rate_table.firing_rates_at_left_sde_high = firing_rates_at_left_sde_high;
cluster_rate_table.firing_rates_at_right_sde_low = firing_rates_at_right_sde_low;
cluster_rate_table.firing_rates_at_right_sde_high = firing_rates_at_right_sde_high;




cluster_rate_table.zscored_peak_pre_rates_left = zscored_peak_pre_rates_left;
cluster_rate_table.zscored_peak_pre_rates_right  = zscored_peak_pre_rates_right ;
cluster_rate_table.zscored_ave_pre_rates_left = zscored_ave_pre_rates_left;
cluster_rate_table.zscored_ave_pre_rates_right = zscored_ave_pre_rates_right;
cluster_rate_table.zscored_ave_post_rates_left_in_replay = zscored_ave_post_rates_left_in_replay;
cluster_rate_table.zscored_ave_post_rates_right_in_replay =zscored_ave_post_rates_right_in_replay;
cluster_rate_table.zscored_ave_post_rates_left = zscored_ave_post_rates_left;
cluster_rate_table.zscored_ave_post_rates_right =zscored_ave_post_rates_right;









save('cluster_rates_during_stops_v2.mat','cluster_rate_table');