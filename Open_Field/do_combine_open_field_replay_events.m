rats = [1 4 6 12 13 14]; %all rats
% rats = [1 6]; % Opto rats
% rats = 4; % Dreads rat
load data_tbl11.mat;
fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures';
downsample_to_match_angle_bins = 0;
look_at_speed_thresholded_events = 1;
convert_radians_to_degrees = 0;
align_to_drink_onset = 1;
align_to_drink_offset = 0;
min_stopping_period_duration = 0;
max_stopping_period_duration = inf;
remove_trials_with_anticipatory_licking = 0;
align_to_start_of_anticipatory_licking = 0;
session_decoding_accuracy_thr = 5;
novel_range = [10.1 10.2];
thrForCategorization_include = 20;
categories_must_be_exclusive = 1;
future_wins = 0;
include_antiPast = 0;
flags_to_include = [11];
% 18 = dreadds; 19 = saline; 13 = laser off; 14 = laser on; what is 20?
% flags_to_exclude = [100 20 19 13]; % MEC inactive
flags_to_exclude = [100 20 18 14]; % MEC active
% flags_to_exclude = [100 20]; % flag 20 = rat sleeping on maze. flag 100 =
% decoding worse than 5 cm.
angle_between_past_future_trajectory_Thr = [0 180];
use_mean_replay_path_angle = 1;
use_distribution_replay_path_angle = 0;
restrict_to_particular_ring = 1;
ring_range = [7 90]; % set to [1 1] to include just the first ring, etc.
restrict_past_future_paths_to_particular_ring = 1;
past_future_ring_range = [7 30];
% first ring is 2 cm away, second ring is 4 cm away, etc...
flags.flags_to_include = flags_to_include;
flags.flags_to_exclude = flags_to_exclude;
must_include_all_flags = 1; % 1: requires that all flags in Flags_to_include are present. 0: requires that at least one flag in Flags_to_include are present.
replay_dispersionThr = 0; % in bins
replay_distanceThr = 0; % in bins
replay_durationThr = 0.05;
replay_startEndDistThr = 0; %20 bins = 40cm;
sd_thr = -inf;
ripple_thr = -inf;
home_trials_only = 0;
away_trials_only = 0;
load_shuffles = 0;
plot_shuffles = 0;
num_trials_away_low = -inf;
num_trials_away_high = inf;
%shuffle_type = 'path_swap_future_past_preserved';
shuffle_type = 'path_swap';
speedThr = inf;
trial_range = [1 300];
windowSize = 0.5;
windowShift = 0.5;
start_time = 0;
end_time = 14;
ylim_deg = [50 120]; %[60 120][
ylim_deg_difference_plot = [-40 10]; %[-10 40]

% Need to run this if data_tbl is not already loaded
% combine_replayEvents_behavioralSimilarity
%%
data_tbl.laser_state_binary(cellfun(@sum,cellfun(@(x) ismember(x,18),data_tbl.session_flags,'UniformOutput',false))==1) = 1; % Dreadds session-flagging it 'laser on'
data_tbl.laser_state_binary(cellfun(@sum,cellfun(@(x) ismember(x,19),data_tbl.session_flags,'UniformOutput',false))==1) = 0; % Saline session-flagging it 'laser on'
data_tbl.ripple_power(isnan(data_tbl.ripple_power)) = -inf;
% If the laser was on at any point in the session, remove any replays detected without the laser on. 
bad_events = intersect(find(cellfun(@sum,cellfun(@(x) ismember(x,14),data_tbl.session_flags,'UniformOutput',false))==1), find(data_tbl.laser_state_binary==0));
data_tbl(bad_events,:) = [];

number_of_flags_met = cellfun(@sum,cellfun(@(x) ismember(x,flags_to_include),data_tbl.session_flags,'UniformOutput',false));
contains_exclusion_flags = cellfun(@sum,cellfun(@(x) ismember(x,flags_to_exclude),data_tbl.session_flags,'UniformOutput',false));
replay = data_tbl(ismember(data_tbl.rat_label, rats) & number_of_flags_met==length(flags.flags_to_include) & contains_exclusion_flags == 0 &...
    data_tbl.dispersion > replay_dispersionThr & ...
    data_tbl.distance > replay_distanceThr  & ...
    data_tbl.start_to_end_distance > replay_startEndDistThr & ...
    data_tbl.duration >= replay_durationThr & ...
    data_tbl.drink_period_time >= min_stopping_period_duration & data_tbl.drink_period_time <= max_stopping_period_duration & ...
    data_tbl.spike_density_power>=sd_thr & ...
    data_tbl.ripple_power>=ripple_thr & ...
    data_tbl.ratSpeed <= speedThr & ...
    data_tbl.drink_period_number >= trial_range(1) & data_tbl.drink_period_number <= trial_range(2) & ... 
    data_tbl.time_since_real_drink_onset >= start_time & data_tbl.time_since_real_drink_onset <= end_time,:);

% convert from radians to degrees
if convert_radians_to_degrees == 1
    replay.mean_abs_angular_displacement_ratLoc = rad2deg(replay.mean_abs_angular_displacement_ratLoc);
    replay.meanAngDisplacement_futPath = rad2deg(replay.meanAngDisplacement_futPath);
    replay.meanAngDisplacement_pastPath = rad2deg(replay.meanAngDisplacement_pastPath);
    replay.angle_between_past_future_trajectory = rad2deg(replay.angle_between_past_future_trajectory);

    if iscell(replay.all_angles_between_past_future_trajectory)
        replay.all_angles_between_past_future_trajectory = cellfun(@(x) rad2deg(x),replay.all_angles_between_past_future_trajectory, 'UniformOutput', false);
        replay.all_angles_between_past_future_trajectory_inner_90 = cellfun(@(x) rad2deg(x),replay.all_angles_between_past_future_trajectory_inner_90, 'UniformOutput', false);
        replay.all_angles_between_past_future_trajectory_inner_80 = cellfun(@(x) rad2deg(x),replay.all_angles_between_past_future_trajectory_inner_80, 'UniformOutput', false);
    else
        replay.all_angles_between_past_future_trajectory = rad2deg(replay.all_angles_between_past_future_trajectory);
        replay.all_angles_between_past_future_trajectory_inner_90 = rad2deg(replay.all_angles_between_past_future_trajectory_inner_90);
        replay.all_angles_between_past_future_trajectory_inner_80 = rad2deg(replay.all_angles_between_past_future_trajectory_inner_80);
    end

    replay.angDisplacement_futPath = cellfun(@(x) rad2deg(x),replay.angDisplacement_futPath, 'UniformOutput', false);
    replay.angDisplacement_pastPath = cellfun(@(x) rad2deg(x),replay.angDisplacement_pastPath, 'UniformOutput', false);

    if load_shuffles==1
        replay.replay_shuffle_past = cellfun(@(x) rad2deg(x),replay.replay_shuffle_past, 'UniformOutput', false);
        replay.replay_shuffle_future = cellfun(@(x) rad2deg(x),replay.replay_shuffle_future, 'UniformOutput', false);
        replay.all_angularDisplacement_past_shuffled = cellfun(@(x) rad2deg(x),replay.all_angularDisplacement_past_shuffled, 'UniformOutput', false);
        replay.all_angularDisplacement_future_shuffled = cellfun(@(x) rad2deg(x),replay.all_angularDisplacement_future_shuffled, 'UniformOutput', false);
    end
end
replay.dispersion = (replay.dispersion).*2; % convert from bins to cm
replay.distance = (replay.distance).*2; % convert from bins to cm
replay.startDistFromRat = replay.startDistFromRat.*2; % convert from bins to cm
replay.start_to_end_distance = (replay.start_to_end_distance).*2; % convert from bins to cm
replay.percent_participation_excitatory = (replay.percent_participation_excitatory).*100; % convert from proportion to percent
replay.percent_participation_inhibitory = (replay.percent_participation_inhibitory).*100; % convert from proportion to percent
replay.slope_dispersion = (replay.dispersion)./replay.duration;
replay.slope_distance = (replay.distance)./replay.duration;
replay.slope_start_to_end_distance = (replay.start_to_end_distance)./replay.duration;

if restrict_to_particular_ring == 1
    replay.meanAngDisplacement_futPath = cellfun(@(x) nanmean(abs(x(ring_range(1):ring_range(2)))), replay.angDisplacement_futPath);
    replay.meanAngDisplacement_pastPath = cellfun(@(x) nanmean(abs(x(ring_range(1):ring_range(2)))), replay.angDisplacement_pastPath);
%     replay.meanAngDisplacement_HD = cellfun(@(x) nanmean(abs(x(ring_range(1):ring_range(2)))), replay.angDisplacement_HD);
%     replay.meanAngDisplacement_HD_and_past_path = cellfun(@(x) nanmean(abs(x(ring_range(1):ring_range(2)))), replay.angDisplacement_HD_and_past_path);
%     replay.meanAngDisplacement_HD_and_future_path = cellfun(@(x) nanmean(abs(x(ring_range(1):ring_range(2)))), replay.angDisplacement_HD_and_future_path);
end
replay.meanAngDisplacement_antiPastPath = 180-replay.meanAngDisplacement_pastPath;
replay.meanAngDisplacement_antiFutPath = 180-replay.meanAngDisplacement_futPath;

% remove shuffles that came BEFORE or AFTER the current trial
% if load_shuffles==1
%     modifyCell = @(shuffle,num_trials) shuffle(num_trials>= num_trials_away_low & num_trials <= num_trials_away_high,:);
%     replay.replay_shuffle_past = cellfun(modifyCell, replay.replay_shuffle_past, replay.num_trials_between_current_stopping_period_and_shuffle, 'UniformOutput', false);
%     replay.replay_shuffle_future = cellfun(modifyCell, replay.replay_shuffle_future, replay.num_trials_between_current_stopping_period_and_shuffle, 'UniformOutput', false);
% end
if load_shuffles==1
    modifyCell = @(shuffle,num_trials) shuffle(num_trials>= num_trials_away_low & num_trials <= num_trials_away_high,:);
    replay.all_angularDisplacement_past_shuffled = cellfun(modifyCell, replay.all_angularDisplacement_past_shuffled, replay.num_trials_between_current_stopping_period_and_shuffle, 'UniformOutput', false);
    replay.all_angularDisplacement_future_shuffled = cellfun(modifyCell, replay.all_angularDisplacement_future_shuffled, replay.num_trials_between_current_stopping_period_and_shuffle, 'UniformOutput', false);
end
if load_shuffles==1
    for i = 1:height(replay)
        if isempty(replay.all_angularDisplacement_past_shuffled{i})
            replay.past_shuffled{i} = nan;
        else
        replay.past_shuffled{i} = nanmean(abs(replay.all_angularDisplacement_past_shuffled{i}(:,ring_range(1):ring_range(2))),2);
        end
        if isempty(replay.all_angularDisplacement_future_shuffled{i})
            replay.future_shuffled{i} = nan;
        else
        replay.future_shuffled{i} = nanmean(abs(replay.all_angularDisplacement_future_shuffled{i}(:,ring_range(1):ring_range(2))),2);
        end
    end
end
if plot_shuffles==1    
inds_with_shuffles = cellfun(@(x) size(x,2), replay.all_angularDisplacement_past_shuffled) == 90;
replay.past_shuffled_selected_rings = cell(height(replay),1);
replay.past_shuffled_selected_rings(inds_with_shuffles) = cellfun(@(x)  nanmean(abs(x(:,ring_range(1):ring_range(2))),2), replay.all_angularDisplacement_past_shuffled(inds_with_shuffles), 'UniformOutput', false);

inds_with_shuffles = cellfun(@(x) size(x,2), replay.all_angularDisplacement_future_shuffled) == 90;
replay.future_shuffled_selected_rings = cell(height(replay),1);
replay.future_shuffled_selected_rings(inds_with_shuffles) = cellfun(@(x)  nanmean(abs(x(:,ring_range(1):ring_range(2))),2), replay.all_angularDisplacement_future_shuffled(inds_with_shuffles), 'UniformOutput', false);
  
replay.all_paths_shuffled_selected_rings = cellfun(@(x,y) [x; y], replay.future_shuffled_selected_rings, replay.past_shuffled_selected_rings, 'UniformOutput', false);

% Do we want all shuffles, or the average of all for each replay?
replay.mean_shuffle_future = cellfun(@nanmean,replay.future_shuffled_selected_rings);
replay.mean_shuffle_past = cellfun(@nanmean,replay.past_shuffled_selected_rings);
replay.mean_shuffle_all = cellfun(@nanmean,replay.all_paths_shuffled_selected_rings);

all_shuffles = nan(height(replay)*100,1); %upper boundary on size
all_shuffles_times = nan(height(replay)*100,1); %upper boundary on size
all_shuffles_home = nan(height(replay)*100,1); %upper boundary on size
start_index = 1;
for i = 1:height(replay)
    n = 0;
    if ~isempty(replay.all_paths_shuffled_selected_rings{i})
    n = size(replay.all_paths_shuffled_selected_rings{i},1);
    all_shuffles(start_index:start_index+n-1) = replay.all_paths_shuffled_selected_rings{i};
    all_shuffles_times(start_index:start_index+n-1) = repmat(replay.time_since_real_drink_onset(i),[n,1]);
    all_shuffles_home(start_index:start_index+n-1) = repmat(replay.home_event(i),[n,1]);
    end 
    start_index = start_index+n;
end
inds_to_keep = intersect(find(~isnan(all_shuffles)), find(~isnan(all_shuffles_times)));
shuffled_replay = table();
shuffled_replay.all_shuffles = all_shuffles(inds_to_keep);
shuffled_replay.time_since_real_drink_onset = all_shuffles_times(inds_to_keep);
shuffled_replay.home_events = all_shuffles_home(inds_to_keep);
end

if restrict_past_future_paths_to_particular_ring == 1
replay.angle_between_past_future_trajectory(:) = nanmean(abs(replay.all_angles_between_past_future_trajectory(:,past_future_ring_range(1):past_future_ring_range(2))),2);
end

replay.angle_between_past_future_trajectory(isnan(replay.angle_between_past_future_trajectory))= -1; % NaN's are messing up my sorting
% Select for replays where the angle between the past and future paths is
% within the desired limits
if iscell(replay.all_angles_between_past_future_trajectory_inner_90)
    replay.all_angles_between_past_future_trajectory_inner_90 = cell2mat(replay.all_angles_between_past_future_trajectory_inner_90);
end
if iscell(replay.all_angles_between_past_future_trajectory_inner_80)
    replay.all_angles_between_past_future_trajectory_inner_80 = cell2mat(replay.all_angles_between_past_future_trajectory_inner_80);
end
if use_mean_replay_path_angle == 1
    replay(replay.angle_between_past_future_trajectory < angle_between_past_future_trajectory_Thr(1) | ...
        replay.angle_between_past_future_trajectory > angle_between_past_future_trajectory_Thr(2),:) = [];
elseif use_distribution_replay_path_angle == 1
    replay = replay(replay.all_angles_between_past_future_trajectory_inner_80(:,1) >= angle_between_past_future_trajectory_Thr(1) & replay.all_angles_between_past_future_trajectory_inner_80(:,2) <= angle_between_past_future_trajectory_Thr(2),:);
end

replay(isnan(replay.meanAngDisplacement_futPath) & isnan(replay.meanAngDisplacement_pastPath),:) = [];
replay.replay = ones(height(replay),1);
replay.future = replay.meanAngDisplacement_futPath <= thrForCategorization_include;
replay.past = replay.meanAngDisplacement_pastPath <= thrForCategorization_include;
replay.home_ending = replay.distance_end_replay_to_home < 7.5;
replay.random_ending = replay.distance_end_replay_to_random_well < 7.5;

replay.more_future = replay.future==1 & (replay.meanAngDisplacement_futPath < replay.meanAngDisplacement_pastPath);
replay.more_past = replay.past==1 & (replay.meanAngDisplacement_pastPath < replay.meanAngDisplacement_futPath);
replay.either_future_or_past = replay.more_future==1 | replay.more_past==1;

if include_antiPast == 1
    replay.antiPast = replay.meanAngDisplacement_antiPastPath <= thrForCategorization_include;
else
    replay.antiPast = zeros(height(replay),1);
end

if categories_must_be_exclusive == 1
    mixed = find(replay.future + replay.past + replay.antiPast > 1);
    replay.future(mixed) = 0;
    replay.past(mixed) = 0;
    replay.antiPast(mixed) = 0;
elseif future_wins == 1
    replay.past(replay.future==1)=0;
end
replay.categorized = replay.future+replay.past+replay.antiPast>0;
replay.uncategorized = ~replay.categorized;

replay.home_ending_future = replay.home_ending & replay.future;
replay.random_ending_future = replay.random_ending & replay.future;

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

%replay(isnan(replay.meanAngDisplacement_futPath) | isnan(replay.meanAngDisplacement_pastPath),:) = [];

%properties = {'meanAngDisplacement_futPath'; 'meanAngDisplacement_pastPath'; 'mean_shuffle_all'; 'mean_shuffle_future'; 'mean_shuffle_past'};
properties = {'meanAngDisplacement_futPath'; 'meanAngDisplacement_pastPath'};

windowSize = 0.5;
windowShift = 0.5;
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
        binned_data.median(time_bin) = nanmedian(binned_data.data{time_bin});
        binned_data.sem(time_bin) = nanstd(binned_data.data{time_bin})/sqrt(sum(~isnan(binned_data.data{time_bin})));
        binned_data.n(time_bin) = sum(~isnan(binned_data.data{time_bin}));        
        binned_data.low(time_bin) = quantile(binned_data.data{time_bin},0.025);
        binned_data.high(time_bin) = quantile(binned_data.data{time_bin},0.975);
        binned_data.sum(time_bin) = nansum(binned_data.data{time_bin});
        binned_data.pcnt(time_bin) = nansum(binned_data.data{time_bin})./length(binned_data.data{time_bin});
    end
    binned_properties.(properties{property}) = binned_data;
end

binned_data = table();
for time_bin = 1:length(bin_edges)
    a = binned_properties.meanAngDisplacement_futPath(time_bin,:).data{:};
    b = binned_properties.meanAngDisplacement_pastPath(time_bin,:).data{:};
    c = a-b;
    binned_data.data{time_bin} = c;
    binned_data.mean(time_bin) = nanmean(c);
    binned_data.sem(time_bin) = nanstd(c)/sqrt(sum(~isnan(c)));
    binned_data.n(time_bin) = sum(~isnan(c));    
    binned_data.low(time_bin) = quantile(c,0.025);
    binned_data.high(time_bin) = quantile(c,0.975);
end
binned_properties.future_minus_past_angle = binned_data;
%% 
keyboard
% Plot angle between replay and future path, angle between replay and past path over time into stopping period.
plot_replay_angles_relative_to_past_future_path % Fig 1,j-k
plot_open_field_replay_rate_by_stopping_period %Fig 2k-m. make sure bins are not overlapping before running this!

plot_time_since_arrival_versus_replay_angle %Fig 2g
plot_figure2_angle_means %Fig 2i; 
plot_fig2i_angles_in_early_or_late_window_2_by_2.m % Fig 2i
combine_replayEvents_behavioralSimilarity_by_session_or_rat_v2 % Figure 2h
plot_open_field_replay_properties_over_stopping_period % Figure 2L, 4D, 4K, supp figure 1
plot_timing_of_pro_retro_replays_early_late_novel_familiar


plot_replay_angles_relative_to_past_future_path % in Fig 3, Fig 4 (collapses across time into stopping period)
plot_crossing_angles_relative_to_past_future_path % (collapses across time into stopping period; Fig S6B)
plot_angle_replay_future_past_diff_two_groups % Panel d from Fig 3, also used in Fig 4.
plot_crossings_future_past_diff_two_groups 

plot_distribution_of_past_versus_future_path_angles % in Fig 3b
plot_past_angle_versus_all_shuffle_figure_3 % In Fig 3

combine_sd_ripple_rates_across_stopping_periods % Supp fig 4

