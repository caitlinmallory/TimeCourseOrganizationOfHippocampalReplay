clearvars -except dayFiles day directory rat windows hand_clustered_only

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

load Experiment_Information.mat
load clusters

if Experiment_Information.spatialDim == 2
load true_drink_periods.mat
else
load Behavior_Data.mat
end

load Position_Data;

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
spikeDensity(:,speed>5)=nan; % remove any periods from spike density where rat was moving.

%% Find replays, using SDEs as candidate events:
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

stopping_period_times = Reward_Epoch_Time_Boundaries_speed_thresholded;
end_zone = Reward_Epoch_Time_Boundaries_endzone;

%%
length_to_plot = 10;
bins_forward = length_to_plot/spikeDensityStepSize;

firing_rates = nan(length(clusters),bins_forward,length(stopping_period_times));
firing_rates_replay = nan(length(clusters),bins_forward,length(stopping_period_times));
firing_rates_sde = nan(length(clusters),bins_forward,length(stopping_period_times));

for i = 1:length(stopping_period_times)
stopping_period_start_time = stopping_period_times(i,1);
[~,stopping_period_start_ind] = min(abs(timeEdges-stopping_period_start_time));
[stopping_period_end_time]= stopping_period_times(i,2);
[~,stopping_period_end_ind] = min(abs(timeEdges-stopping_period_end_time));
inds_to_plot = stopping_period_start_ind:min((stopping_period_start_ind + (bins_forward-1)),stopping_period_end_ind); 
if any(inds_to_plot<0)
    continue 
end
if any(inds_to_plot>length(spikeDensity))
    continue
end

sub_rates = nan(length(clusters),bins_forward);
sub_rates(:,1:length(inds_to_plot)) = spikeDensity(:,inds_to_plot);

sub_rates_replay = nan(length(clusters),bins_forward);
sub_rates_replay(:,1:length(inds_to_plot)) = spikeDensity_during_replay(:,inds_to_plot);

sub_rates_sde = nan(length(clusters),bins_forward);
sub_rates_sde(:,1:length(inds_to_plot)) = spikeDensity_during_sdes(:,inds_to_plot);


firing_rates(:,:,i) =sub_rates ;
firing_rates_replay(:,:,i) = sub_rates_replay;
firing_rates_sde(:,:,i) = sub_rates_sde;
end

firing_rates_participating_replay = nan(size(firing_rates_replay));
for cell = 1:size(firing_rates_replay,1)
    for column = 1:size(firing_rates_replay,3)
        stopping_period = squeeze(firing_rates_replay(cell,:,column));
        stopping_period(isnan(stopping_period)) = 5000;
        diffs = [0 diff(stopping_period)];
        replay_starts = find(diffs<-4000);
        replay_ends = find(diffs>4000)-1;
        if length(replay_ends)<length(replay_starts)
            replay_ends = [replay_ends size(firing_rates_replay,2)];
        end
         
        num_replays = length(replay_starts);
        for replay = 1:num_replays
            if sum(stopping_period(replay_starts(replay):replay_ends(replay))) == 0
                stopping_period(replay_starts(replay):replay_ends(replay)) = nan;
            end
        end
        stopping_period(stopping_period==5000)=nan;
        firing_rates_participating_replay(cell,:,column) = stopping_period;
    end
end

firing_rates([clusters.Excitatory]==0,:,:)=nan;
firing_rates_replay([clusters.Excitatory]==0,:,:)=nan;
firing_rates_participating_replay([clusters.Excitatory]==0,:,:)=nan;
firing_rates_sde([clusters.Excitatory]==0,:,:)=nan;

mean_frs = nan(length(clusters),1);
std_frs = nan(length(clusters),1);
mean_frs_replay = nan(length(clusters),1);
std_frs_replay = nan(length(clusters),1);
mean_frs_participating_replay = nan(length(clusters),1);
std_frs_participating_replay = nan(length(clusters),1);
mean_frs_sde = nan(length(clusters),1);
std_frs_sde = nan(length(clusters),1);

for i=1:length(clusters)
    sub = squeeze(firing_rates(i,:,:));
    mean_frs(i) = nanmean(sub(:));
    std_frs(i) = nanstd(sub(:));

    sub = squeeze(firing_rates_replay(i,:,:));
    mean_frs_replay(i) = nanmean(sub(:));
    std_frs_replay(i) = nanstd(sub(:));

    sub = squeeze(firing_rates_participating_replay(i,:,:));
    mean_frs_participating_replay(i) = nanmean(sub(:));
    std_frs_participating_replay(i) = nanstd(sub(:));


    sub = squeeze(firing_rates_sde(i,:,:));
    mean_frs_sde(i) = nanmean(sub(:));
    std_frs_sde(i) = nanstd(sub(:));
end
zscored_firing_rates = (firing_rates-mean_frs)./std_frs;
zscored_firing_rates_replay = (firing_rates_replay-mean_frs_replay)./std_frs_replay;
zscored_firing_rates_participating_replay = (firing_rates_participating_replay-mean_frs_replay)./std_frs_participating_replay;
zscored_firing_rates_sde = (firing_rates_sde-mean_frs_sde)./std_frs_sde;

firing_rates_at_left = squeeze(nanmean(firing_rates(:,:,end_zone==1),3));
firing_rates_at_right = squeeze(nanmean(firing_rates(:,:,end_zone==2),3));
zscored_firing_rates_at_left = squeeze(nanmean(zscored_firing_rates(:,:,end_zone==1),3));
zscored_firing_rates_at_right = squeeze(nanmean(zscored_firing_rates(:,:,end_zone==2),3)); 

firing_rates_at_left_participating_replay = squeeze(nanmean(firing_rates_participating_replay(:,:,end_zone==1),3));
firing_rates_at_right_participating_replay = squeeze(nanmean(firing_rates_participating_replay(:,:,end_zone==2),3));
zscored_firing_rates_at_left_participating_replay = squeeze(nanmean(zscored_firing_rates_participating_replay(:,:,end_zone==1),3));
zscored_firing_rates_at_right_participating_replay = squeeze(nanmean(zscored_firing_rates_participating_replay(:,:,end_zone==2),3)); 

firing_rates_at_left_replay = squeeze(nanmean(firing_rates_replay(:,:,end_zone==1),3));
firing_rates_at_right_replay = squeeze(nanmean(firing_rates_replay(:,:,end_zone==2),3));
zscored_firing_rates_at_left_replay = squeeze(nanmean(zscored_firing_rates_replay(:,:,end_zone==1),3));
zscored_firing_rates_at_right_replay = squeeze(nanmean(zscored_firing_rates_replay(:,:,end_zone==2),3)); 

firing_rates_at_left_sde = squeeze(nanmean(firing_rates_sde(:,:,end_zone==1),3));
firing_rates_at_right_sde = squeeze(nanmean(firing_rates_sde(:,:,end_zone==2),3));
zscored_firing_rates_at_left_sde = squeeze(nanmean(zscored_firing_rates_sde(:,:,end_zone==1),3));
zscored_firing_rates_at_right_sde = squeeze(nanmean(zscored_firing_rates_sde(:,:,end_zone==2),3)); 

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
cluster_rate_table.firing_rates_at_left_participating_replay = firing_rates_at_left_participating_replay;
cluster_rate_table.firing_rates_at_right_participating_replay = firing_rates_at_right_participating_replay;
cluster_rate_table.zscored_firing_rates_at_left_participating_replay = zscored_firing_rates_at_left_participating_replay;
cluster_rate_table.zscored_firing_rates_at_right_participating_replay = zscored_firing_rates_at_right_participating_replay;
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


save('cluster_rates_during_stops_v1.mat','cluster_rate_table');