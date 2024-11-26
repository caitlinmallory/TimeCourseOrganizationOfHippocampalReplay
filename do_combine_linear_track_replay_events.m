reload_data = 0;
save_data_table = 0;
plot_proportions = 1;
plot_features_over_stopping_period = 1;
% Select which replays you want to look at:
candidate_events_to_plot = 'spike_filtered';
final_fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures';
posterior_difference_to_use = [{'percent_non_local_posterior_in_left_map'},{'percent_non_local_posterior_in_right_map'}];%

downsample_to_match_time_into_stopping_period_and_session = 0;
downsample_to_match_time_into_session = 0;
downsample_to_match_time_into_stopping_period = 0;

% MANUSCRIPT:
shuffle_labels = 0;
override_coverage_thr = 0.2;
override_weighted_r_thr = 0.6;
override_posterior_diff_thr = 0.33;
num_cells_participating_thr = 10;
fraction_cells_participating_thr = 0;
override_sde_thr = -inf;
override_ripple_thr = -inf;
combine_control_epoch_data = 0;
rats = [1:11];
% rats = [4 5]; % Laser control (GFP) rats.
% rats = [1 2 3 6]; % Jaws-GFP rats


fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures';
binSize = 2;
novel_range = [10.1 inf]; %inclusive range

flags_to_include = [11];
flags_to_exclude = [0 1002 1003];
flags.flags_to_include = flags_to_include;
flags.flags_to_exclude = flags_to_exclude;
must_include_all_flags = 1; % 1: requires that all flags in Flags_to_include are present. 0: requires that at least one flag in Flags_to_include are present.
long_replays_only = 1;

engagment_time = 5;
time_into_session_thr = inf; % set to inf to including everything
endzone_thr = 1; % set to 0 to include everthing;
distance_from_end_tight_threshold = 30; % set to inf to include everthing
engagement_classifications_to_use = [-1 1 0]; % set to  [-1,1,0] to include everything

upper_limit_time_since_reward_zone_entry = 10; % inf is most permissive
lower_limit_time_since_reward_zone_entry = 0; % -inf is most permissive

% Session quality criterion:
min_num_replays_per_session = 0; % this will remove sessions with fewer than min_num_replays_per_session in each epoch (Laser on and Laser off) from analysis
pnct_directional_decoding_correct_thr = 0.70;
median_decoding_error_thr = 5;

time_to_run_threshold = inf;
num_passes_criterion = 0; % 10 laps
if strcmp(candidate_events_to_plot,'ripple_filtered')
    ripple_power_thr = override_ripple_thr;
    sde_amplitude_thr = -inf;
elseif strcmp(candidate_events_to_plot,'spike_filtered')
    ripple_power_thr = override_ripple_thr;
    sde_amplitude_thr = override_sde_thr;
end

drinking_periods_to_examine = [0 1];

Rat_Names = {'Clover','Bo','MEC1','Bolt','Dash','CM1','Janni','Harpy','Imp','W18','W19'};
experimental_rats = [1 2 3 6]; % opto; with Jaws
control_rats1 = [4 5 ]; % opto, gfp only
control_rats2 = [7 8 9 10 11]; % no laser, no injections

if strcmp(candidate_events_to_plot,'ripple_filtered')
    event_choice = 5;
elseif strcmp(candidate_events_to_plot,'spike_filtered')
    event_choice = 6;
end
%%
if reload_data==1
    combine_replayPropertiesAcrossSessions
else
    load(['/home/caitlin/Data/Processed_Data/data_tbl_' candidate_events_to_plot]);
    load(['/home/caitlin/Data/Processed_Data/t_session_' candidate_events_to_plot]);
end

if save_data_table == 1
    save(['/home/caitlin/Data/Processed_Data/data_tbl_' candidate_events_to_plot],'data_tbl','-v7.3')
    save(['/home/caitlin/Data/Processed_Data/t_session_' candidate_events_to_plot],'t_session','-v7.3')
end
%%
[criterion_str, duration_thr, max_jump_distance_thr, coverage_thr_percent_of_track, coverage_thr, weighted_r_thr, posterior_in_map_thr] = load_replay_criterion(long_replays_only,override_coverage_thr, override_weighted_r_thr, override_posterior_diff_thr);
coverage_thr = coverage_thr/binSize;
unique_session_id = unique(data_tbl.unique_session_id);

% make a copy of the data and filter according to your desired criterion
t = data_tbl;
t.weighted_r = abs(t.weighted_r);
t.end_distance_from_rat(t.incongruent_with_rat_location==1) = nan;
t.range = t.range*2; % conver to bins
t.mean_CI95_posterior_cropped_normalized = 1-t.mean_CI95_posterior_cropped_normalized;
t.mean_HPD95_posterior_cropped_normalized = t.mean_HPD95_posterior_cropped_normalized.*100;

pull_ripple_sd_properties_from_decoder_events;

% select the data from the rats/sessions you want to analayze
t(~ismember(t.rat_label,rats),:) = [];
t_session(~ismember(t_session.rat_num,rats),:) = [];

for flag = 1:length(flags_to_include)
    t(cellfun(@(x) ismember(flags_to_include(flag),x), cellfun(@round, t.session_flags, 'UniformOutput', false))==0,:) = [];
    t_session(cellfun(@(x) ismember(flags_to_include(flag),x), cellfun(@round, t_session.session_flags, 'UniformOutput', false))==0,:) = [];

    if flags_to_include(flag) == 10 % novelty flag
        % step 1: find the flag containing 10
        % step 2: check whether this flag is between the lower and upper
        % bounds desired
        [inds2keep] = [];
        for ind = 1:height(t)
            session_flags = t.session_flags{ind};
            [~,index_of_novelty_flag] = find(floor(session_flags) == 10);
            if session_flags(index_of_novelty_flag) >= novel_range(1) && session_flags(index_of_novelty_flag) <= novel_range(2)
                inds2keep = [inds2keep; ind];
            end
        end

        t = t(inds2keep,:);

        [inds2keep] = [];
        for ind = 1:height(t_session)
            session_flags = t_session.session_flags{ind};
            [~,index_of_novelty_flag] = find(floor(session_flags) == 10);
            if session_flags(index_of_novelty_flag) >= novel_range(1) && session_flags(index_of_novelty_flag) <= novel_range(2)
                inds2keep = [inds2keep; ind];
            end
        end

        t_session = t_session(inds2keep,:);
    end

end
for flag = 1:length(flags_to_exclude)
    t(cellfun(@(x) ismember(flags_to_exclude(flag),x), cellfun(@round, t.session_flags, 'UniformOutput', false))==1,:) = [];
    t_session(cellfun(@(x) ismember(flags_to_exclude(flag),x), cellfun(@round, t_session.session_flags, 'UniformOutput', false))==1,:) = [];
end

% If combine_control_epoch_data is set to one, the laser OFF epochs of control rats can be combined with the laser OFF epochs of exp. rats
if combine_control_epoch_data ==1
    % keep the laser off epochs from the control animals, but remove
    % the laser-on epochs.
    t(t.laser_state == 1 & ismember(t.rat_label,control_rats1),:) = [];

end
% Add flag for experimental versus control rats
t.control_experimental_flag(ismember(t.rat_label,experimental_rats)) = 1;
t.control_experimental_flag(ismember(t.rat_label,control_rats2)) = 2;
% Add properties
t.disengaged = zeros(height(t),1);
t.disengaged(t.time_since_reward_zone_entry > engagment_time & t.time_till_reward_zone_exit > engagment_time) = 1;
t.engaged_entry = zeros(height(t),1);
t.engaged_entry(t.time_since_reward_zone_entry > 0 & t.time_since_reward_zone_entry < engagment_time) = 1;
t.engaged_exit = zeros(height(t),1);
t.engaged_exit(t.time_till_reward_zone_exit < engagment_time) = 1;
t.engagement_classification = nan(height(t),1);
t.engagement_classification(t.disengaged==1) = 0;
t.engagement_classification(t.engaged_entry==1) = 1;
t.engagement_classification(t.engaged_exit==1) = -1;

t.meets_engagement_requirement(ismember(t.engagement_classification,engagement_classifications_to_use),:) = 1;
t.meets_time_since_reward_zone_entry_requirement = zeros(height(t),1);
t.meets_time_since_reward_zone_entry_requirement(t.time_since_reward_zone_entry < upper_limit_time_since_reward_zone_entry & ...
    t.time_since_reward_zone_entry > lower_limit_time_since_reward_zone_entry) = 1;
t.meets_in_endzone_requirement = ones(height(t),1);
t.meets_in_endzone_requirement(t.in_reward_zone == 0) = 0;
t.distance_from_end(t.reward_zone == 1) = t.distance_from_end_1(t.reward_zone == 1);
t.distance_from_end(t.reward_zone == 2) = t.distance_from_end_2(t.reward_zone == 2);
t.meets_tight_in_endzone_requirement = ones(height(t),1);
t.meets_tight_in_endzone_requirement(t.distance_from_end > distance_from_end_tight_threshold) = 0;
t.meets_drinking_requirement = zeros(height(t),1);
t.meets_drinking_requirement(ismember(t.drinking,drinking_periods_to_examine)) = 1;
t.posterior_difference = abs(t.(posterior_difference_to_use{1}) - t.(posterior_difference_to_use{2}));

% good_candidate_event = all events that meet position/timing
% criterion (but nothing about weighted correlation, posterior
% etc). replay = meets all criterion; replay = meets all
% criterion and is long
t.good_candidate_event = ...
    t.replay_spikeDensity_power >= sde_amplitude_thr & ...
    t.replay_ripple_power >= ripple_power_thr & ...
    t.duration >= duration_thr  & ...
    t.time_in_session < time_into_session_thr & ...
    t.in_reward_zone >= endzone_thr & ...
    t.meets_tight_in_endzone_requirement == 1 & ...
    t.meets_engagement_requirement == 1 & ...
    t.meets_time_since_reward_zone_entry_requirement == 1 & ...
    t.meets_drinking_requirement == 1 & ...
    t.num_of_cells_participating >= num_cells_participating_thr & ...
    t.fraction_of_cells_participating >= fraction_cells_participating_thr & ...
    t.session_percent_correct_directional_assignment >= pnct_directional_decoding_correct_thr & ...
    t.session_median_decoding_error <= median_decoding_error_thr & ...
    t.laser_state_matches_intended_laser_state_for_stopping_period == 1; 

t.replay = ...
    t.good_candidate_event == 1 & ...
    abs(t.weighted_r)>= weighted_r_thr & ...
    t.posterior_difference>= posterior_in_map_thr & ...
    t.max_jump_distance <= max_jump_distance_thr & ...
    t.range_normalized >= coverage_thr_percent_of_track & ...
    t.range >= coverage_thr;

t.forward_replay_laser_on = t.replay == 1 & t.direction == 1 & t.laser_state == 1;
t.forward_replay_laser_off = t.replay== 1 & t.direction == 1 & t.laser_state == 0;
t.reverse_replay_laser_on = t.replay == 1 & t.direction == 2 & t.laser_state == 1;
t.reverse_replay_laser_off = t.replay == 1 & t.direction == 2 & t.laser_state == 0;
t.past_map_replay = (t.reward_zone==1 & t.best_map == 1) | (t.reward_zone == 2 & t.best_map==2);
t.future_map_replay = (t.reward_zone==1 & t.best_map ==2) | (t.reward_zone == 2 & t.best_map == 1);

for i = 1:height(t_session)
    t_session.num_replays(i) = length(find(t.replay == 1 & t.unique_session_id == t_session.unique_id(i)));
    t_session.num_replays_laser_on(i) = length(find(t.replay == 1 & t.laser_state == 1 &  t.unique_session_id == t_session.unique_id(i)));
    t_session.num_replays_laser_off(i) = length(find(t.replay == 1 & t.laser_state == 0 &  t.unique_session_id == t_session.unique_id(i)));
end

bad_session_inds = unique([find(t_session.num_passes < num_passes_criterion); find(t_session.num_replays_laser_on < min_num_replays_per_session); find(t_session.num_replays_laser_off < min_num_replays_per_session)]);
bad_session_unique_ids = t_session.unique_id(bad_session_inds);

t_session(bad_session_inds,:) = [];
t(ismember(t.unique_session_id,bad_session_unique_ids),:) = [];

t_session.ratio_num_replays_session_off_on = t_session.num_replays_laser_off./t_session.num_replays_laser_on;

unique_session_ids = setdiff(t_session.unique_id,bad_session_unique_ids);

%%
classifiable_replays_laser_off = length(find(t.replay==1 & t.laser_state == 0));
classifiable_replays_laser_on = length(find(t.replay==1 & t.laser_state == 1));

forward_congruent_laser_off = length(find(t.replay==1 & t.laser_state==0 & t.direction==1 & t.congruent_with_rat_location==1));
reverse_congruent_laser_off = length(find(t.replay==1 & t.laser_state==0 & t.direction==2 & t.congruent_with_rat_location==1));

forward_congruent_laser_on = length(find(t.replay==1 & t.laser_state==1  & t.direction==1 & t.congruent_with_rat_location==1));
reverse_congruent_laser_on = length(find(t.replay==1 & t.laser_state==1 & t.direction==2 & t.congruent_with_rat_location==1));

forward_incongruent_laser_off = length(find(t.replay==1 & t.laser_state==0 & t.direction==1 & t.incongruent_with_rat_location==1));
reverse_incongruent_laser_off = length(find(t.replay==1 & t.laser_state==0 & t.direction==2 & t.incongruent_with_rat_location==1));

forward_incongruent_laser_on = length(find(t.replay==1 & t.laser_state==1 & t.direction==1 & t.incongruent_with_rat_location==1));
reverse_incongruent_laser_on = length(find(t.replay==1 & t.laser_state==1 & t.direction==2 & t.incongruent_with_rat_location==1));

%Classify all replays into: congruent-for, incongruent-for, congruent-rev,
%or incongruent-rev
rc_0 = reverse_congruent_laser_off;
ri_0 = reverse_incongruent_laser_off;
fc_0 = forward_congruent_laser_off;
fi_0 = forward_incongruent_laser_off;
cr_0 = classifiable_replays_laser_off;

rc_1 = reverse_congruent_laser_on;
ri_1 = reverse_incongruent_laser_on;
fc_1 = forward_congruent_laser_on;
fi_1 = forward_incongruent_laser_on; 
cr_1 = classifiable_replays_laser_on;


if plot_proportions == 1
    figure('Position',[600 600 125 125])
    colors = [.4660 0.6740 0.1880; 0.4940 0.1840 0.5560];
    transparency_pcnt = 0.3;
    colors_light = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
        [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];
    bar(1,(fc_0/cr_0).*100,'facecolor',colors(1,:),'EdgeColor','none'); hold on;
    bar(2,(rc_0/cr_0).*100,'facecolor',colors(2,:),'EdgeColor','none'); hold on;
    bar(3.5,(fi_0/cr_0).*100,'facecolor',colors_light(1,:),'EdgeColor','none'); hold on;
    bar(4.5,(ri_0/cr_0).*100,'facecolor',colors_light(2,:),'EdgeColor','none'); hold on;
    xticks([1.5,4])
    box off
    ylim([0 50])
    yticks([0:10:50])
    xticklabels({})
    set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
    saveas(gcf,fullfile(final_fig_path,['ForRevConIncongruent' criterion_str num2str(rats)]),'jpeg')
    saveas(gcf,fullfile(final_fig_path,['ForRevConIncongruent' criterion_str num2str(rats)]),'pdf')    
end

%% Optionally, downsample to match the time into the stopping period and session of laser off versus laser on replays
t_replay_og = t(t.replay==1,:);
t_congruent_replay_og = t(t.replay==1 & t.congruent_with_rat_location==1,:);
t_incongruent_replay_og = t(t.replay==1 & t.incongruent_with_rat_location==1,:);
t_candidate_event_og = t(t.good_candidate_event==1,:);
t_replay = [];
t_candidate_events = [];

unique_sessions = unique(t.unique_session_id);
control_experimental_groups = unique(t.control_experimental_flag);

for groups = 1:length(control_experimental_groups)
    group = control_experimental_groups(groups);

    t_replay_og_session = t(t.replay==1 & t.control_experimental_flag==group,:);
    t_candidate_events_og_session = t(t.good_candidate_event==1 & t.control_experimental_flag==group,:);

    names_of_events_to_downsample = {'replay','candidate_event'};
    events_to_downsample = [{t_replay_og_session};{t_candidate_events_og_session}];
    events_ds = cell(2,1);

    for i = 1:length(events_to_downsample)
        events = events_to_downsample{i};

        x0 = events.time_in_session(events.laser_state==0);
        x1 = events.time_in_session(events.laser_state==1);
        y0 = events.time_since_reward_zone_entry(events.laser_state==0);
        y1 = events.time_since_reward_zone_entry(events.laser_state==1);

        if downsample_to_match_time_into_stopping_period_and_session == 1
           [downsampleindices0, downsampleindices1] = downsample_match_time_in_session_and_stopping_period(x0,y0,x1,y1);
            saveas(gcf,fullfile(fig_path,[names_of_events_to_downsample{i} 'downsampled to match time in session and stopping period']),'jpeg')
            events_ds_0 = events(events.laser_state==0,:); events_ds_0 = events_ds_0(downsampleindices0,:);
            events_ds_1 = events(events.laser_state==1,:); events_ds_1 = events_ds_1(downsampleindices1,:);
            events_ds{i} = [events_ds_0; events_ds_1];
        elseif downsample_to_match_time_into_session == 1
            [downsampleindices0, downsampleindices1] = downsample_match_time_in_session(x0,y0,x1,y1);
            saveas(gcf,fullfile(fig_path,[names_of_events_to_downsample{i} 'downsampled to match time in session']),'jpeg')
            events_ds_0 = events(events.laser_state==0,:); events_ds_0 = events_ds_0(downsampleindices0,:);
            events_ds_1 = events(events.laser_state==1,:); events_ds_1 = events_ds_1(downsampleindices1,:);
            events_ds{i} = [events_ds_0; events_ds_1];
        elseif downsample_to_match_time_into_stopping_period == 1
            [downsampleindices0, downsampleindices1] = downsample_match_time_in_stopping_period(x0,y0,x1,y1);
            saveas(gcf,fullfile(fig_path,[names_of_events_to_downsample{i} 'downsampled to match time in stopping period']),'jpeg')
            events_ds_0 = events(events.laser_state==0,:); events_ds_0 = events_ds_0(downsampleindices0,:);
            events_ds_1 = events(events.laser_state==1,:); events_ds_1 = events_ds_1(downsampleindices1,:);
            events_ds{i} = [events_ds_0; events_ds_1];
        else
            events_ds{i} = events;
        end

    end
    t_replay_session = events_ds{1};
    t_candidate_events_session = events_ds{2};

    t_replay = [t_replay; t_replay_session];
    t_candidate_events = [t_candidate_events; t_candidate_events_session];
end

%%
if plot_replay_properties == 1
    plot_replay_properties_laser_on_versus_off %Fig s8D
end

if plot_features_over_stopping_period == 1
    plot_replay_properties_over_time %Fig 2E
end

plot_timing_of_forward_reverse_early_late_novel_familiar % Fig 2F

plot_forward_reverse_timing_individual_animals_or_sessions % Fig 2A

plot_pcnt_replays_w_opposing_content_versus_time_between_them % Fig 3J




% Overall timing:
% forward_replays = t_congruent_replay_og(t_congruent_replay_og.direction==1 & t_congruent_replay_og.laser_state==0,:)
% reverse_replays = t_congruent_replay_og(t_congruent_replay_og.direction==2 & t_congruent_replay_og.laser_state==0,:)
% mean(forward_replays.time_since_reward_zone_entry)
% std(forward_replays.time_since_reward_zone_entry)./(sqrt(height(forward_replays)))
% 
% quantile(forward_replays.time_since_reward_zone_entry,0.25)
% quantile(forward_replays.time_since_reward_zone_entry,0.50)
% quantile(forward_replays.time_since_reward_zone_entry,0.75)
% 
% mean(reverse_replays.time_since_reward_zone_entry)
% std(reverse_replays.time_since_reward_zone_entry)./(sqrt(height(reverse_replays)))
% 
% [p,h,z] = ranksum(forward_replays.time_since_reward_zone_entry, reverse_replays.time_since_reward_zone_entry)
% 
% quantile(reverse_replays.time_since_reward_zone_entry,0.25)
% quantile(reverse_replays.time_since_reward_zone_entry,0.50)
% quantile(reverse_replays.time_since_reward_zone_entry,0.75)