% Analyzes the rate of replays over the stopping period.
% The number of replays in 0.5 time windows has been previously computed
% for each animal/session. This code turns hist counts into rates over time
% and plots these curves.

load_data=1;
plot_individual_graphs = 1;
align_to_stopping_period_start = 0;
align_to_stopping_period_departure = 0;
align_to_reward_onset = 0;
align_to_reward_offset = 1;

min_stopping_period_duration = 0;
max_stopping_period_duration = inf;

rats = 1:11;
trials_to_include = [1 300]; % Can restrict trials, but set to a large number to include everything.
session_halves_to_include = [1 2]; % Can set to 1 to just look at the first half of the session, or 2 to just look at the secound half.
fig_ylim = [0 0.08];
fig_ylim_diff = [-0.04 0.06];

must_include_all_flags = 1; % 1: requires that all flags in Flags_to_include are present. 0: requires that at least one flag in Flags_to_include are present.
novel_range = [10.2 10.9];
flags.flags_to_include = [11];
% 18 = dreadds; 19 = saline; 13 = laser off; 14 = laser on
flags.flags_to_exclude = [1003 1002 0]; %1003 and 1002 are sessions with increased or omitted reward from Ambrose et al. Do not include these. 0 is a couple sessions where laser light was delivered for many laps in a row. Don't include these
hist_bin_width = 0.5; %seconds
reward_size = [0 inf];
% Session quality criterion:
min_num_replays_per_session = 0; % this will remove sessions with fewer than min_num_replays_per_session in each epoch (Laser on and Laser off) from analysis
pnct_directional_decoding_correct_thr = 0.70;
median_decoding_error_thr = 5;

smooth_rate_plots = 0;
smoothing_sigma = 1; % bins

Rat_Names = {'Clover','Bo','MEC1','Bolt','Dash','CM1','Janni','Harpy','Imp','W18','W19'};
experimental_rats = [1 2 3 6]; % opto; with Jaws
control_rats1 = [4 5 ]; % opto, gfp only
control_rats2 = [7 8 9 10 11]; % no laser, no injections
session_id = 0;

% Where to save the data:
fig_path_prefix = '/home/caitlin/Data/Processed_Data/Replay_summary';

fig_path = fullfile(fig_path_prefix, 'replay_rates_over_stopping_period');
final_fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures';
% make directory to save figs if it doesn't exist
if ~exist(fig_path,'dir')
    mkdir(fig_path)
end

filter_length = smoothing_sigma*6;
if(mod(filter_length,2)==0) % iseven
    filter_length = filter_length+1;
end

combined_table_first_event_times = table();
combined_table_last_event_times = table();
combined_table = table();
combined_hists = table();

if load_data==1
    for rat = rats
        load_linear_track_session_list
        for day = 1:length(dayFiles)
            display([Rat_Names{rat},' Day ', dayFiles{day}])
            cd(fullfile(directory,dayFiles{day}))
            load Experiment_Information
            load session_wide_properties
            if ~(Experiment_Information.reward_size >= reward_size(1) && Experiment_Information.reward_size < reward_size(2))
                continue
            end
            if percent_correct_directional_assigment < pnct_directional_decoding_correct_thr
                continue
            end
            if median_decoding_error > median_decoding_error_thr
                continue
            end
            sessions_that_meet_criterion_day = check_for_relevant_sessions(flags,1,Experiment_Information,must_include_all_flags,1,novel_range);
            for session_count = 1:length(sessions_that_meet_criterion_day)
                session_id = session_id + 1;
                sub_session_num = sessions_that_meet_criterion_day(session_count);
                session_flags = Experiment_Information.Segments(sub_session_num).Flags;

                if sum(ismember(session_flags,18)) == 1
                    saline = 1;
                else
                    saline = 0;
                end

                if sum(ismember(session_flags,19)) == 1
                    cno = 1;
                else
                    cno = 0;
                end

                if sum(ismember(session_flags,1)) == 1
                    reward_stim = 1;
                else
                    reward_stim = 0;
                end

                if sum(ismember(session_flags,2)) == 1
                    run_stim = 1;
                else
                    run_stim = 0;
                end

                if sum(ismember(session_flags,10)) == 1
                    novel = 1;
                else
                    novel = 0;
                end

                if align_to_stopping_period_start==1
                    load('replay_by_stopping_period.mat');
                end
                if align_to_stopping_period_departure==1
                    load('replay_by_stopping_period_aligned_to_departure.mat')
                end
                if align_to_reward_onset==1
                    load('replay_by_stopping_period_aligned_to_drink_onset.mat')
                end
                if align_to_reward_offset==1
                    load('replay_by_stopping_period_aligned_to_drink_offset.mat')
                end

                % add some session properties to the table
                % if one of the ends had increased reward (Ambrose et al
                % experiment), indicated which one.
                if isfield(Experiment_Information,'increased_reward_end')
                    t2.increased_reward_end(:) = repmat(Experiment_Information.increased_reward_end,size(t2,1),1);
                else
                    t2.increased_reward_end(:) = nan(size(t2,1),1);
                end
                t2.session_label(:) = repmat({[Rat_Names{rat},' Day ', dayFiles{day}]},size(t2,1),1);
                t2.rat_label(:) = repmat(rat,size(t2,1),1);
                t2.day_label(:) = repmat(day,size(t2,1),1);
                t2.unique_session_id(:) = repmat(session_id,size(t2,1),1);
                t2.segment_label(:) = repmat(sub_session_num,size(t2,1),1);
                t2.novel_label(:) = repmat(novel,size(t2,1),1);
                t2.run_stim_label(:) = repmat(run_stim,size(t2,1),1);
                t2.reward_stim_label(:) = repmat(reward_stim,size(t2,1),1);
                t2.saline_label(:) = repmat(saline,size(t2,1),1);
                t2.cno_label(:) = repmat(cno,size(t2,1),1);

                if align_to_stopping_period_start == 1
                    combined_table_first_event_times = [combined_table_first_event_times; first_event_time_table];
                    combined_table_last_event_times = [combined_table_last_event_times; last_event_time_table];
                end
                combined_table = [combined_table; t2];
                combined_hists = [combined_hists; struct2table(event_hists)];
            end

        end
    end
    sz = combined_hists.sde;
    sz=size(sz,2);
    time_to_plot = 10; % sec
    hist_edges = 0:hist_bin_width:sz;
    %% saline sessions will be labeled as 'laser OFF' and cno sessions 'laser
    % on'
    combined_table.laser_state(combined_table.saline_label == 1) = 0;
    combined_table.laser_state(combined_table.cno_label == 1) = 1;

else
    load('combined_table.mat');
    load('combined_hists.mat');
    load('combined_table_first_event_times');
    load('combined_table_last_event_times');
end
%% Remove stopping periods that were too long or too short
rows_to_remove = find(combined_table.duration < min_stopping_period_duration | combined_table.duration > max_stopping_period_duration);
combined_table(rows_to_remove,:) = [];
combined_hists(rows_to_remove,:) = [];
if align_to_stopping_period_start==1
    combined_table_first_event_times(rows_to_remove,:) = [];
    combined_table_last_event_times(rows_to_remove,:) = [];
end

unique_sessions = unique(combined_table.unique_session_id);
for i = 1:height(unique_sessions)
    inds = find(combined_table.unique_session_id == unique_sessions(i));
    num_trials = length(inds);
    first_half_inds=inds(1:round(length(inds)/2));
    second_half_inds = inds(round(length(inds)/2+1:end));
    combined_table.session_half(first_half_inds)=1;
    combined_table.session_half(second_half_inds)=2;
end

%% Limit to the desired trial range
rows_to_remove = find(~ismember(combined_table.session_half,session_halves_to_include));
    combined_table(rows_to_remove,:) = [];
    combined_hists(rows_to_remove,:) = [];
    combined_table_first_event_times(rows_to_remove,:) = [];
    combined_table_last_event_times(rows_to_remove,:) = [];
%% Limit to the desired trial range
rows_to_remove = find(combined_table.pass_number < trials_to_include(1) | combined_table.pass_number > trials_to_include(2));
    combined_table(rows_to_remove,:) = [];
    combined_hists(rows_to_remove,:) = [];
    combined_table_first_event_times(rows_to_remove,:) = [];
    combined_table_last_event_times(rows_to_remove,:) = [];
%% Adding a few more event types to the table:
combined_table.replays_spike = combined_table.reverse_replays_spike + combined_table.forward_replays_spike;
combined_table.replays_ripple = combined_table.reverse_replays_ripple + combined_table.forward_replays_ripple;

event_types = [{'replays_spike'},{'replays_ripple'},...
    combined_hists.Properties.VariableNames];

combined_table_rate = table();
for event = 1:length(event_types)
    combined_table_rate.([event_types{event}]) = combined_table.(event_types{event})./combined_table.duration;
end

%% Adding a few more event types to the hists:
combined_hists.replays_spike = combined_hists.reverse_replays_spike + combined_hists.forward_replays_spike;
combined_hists.replays_ripple = combined_hists.reverse_replays_ripple + combined_hists.forward_replays_ripple;
combined_hists.forward_minus_reverse_spike = combined_hists.forward_replays_spike - combined_hists.reverse_replays_spike;
combined_hists.congruent_forward_minus_reverse_spike = combined_hists.forward_congruent_replays_spike - combined_hists.reverse_congruent_replays_spike;
combined_hists.incongruent_forward_minus_reverse_spike = combined_hists.forward_incongruent_replays_spike - combined_hists.reverse_incongruent_replays_spike;
combined_hists.forward_minus_reverse_ripple = combined_hists.forward_replays_ripple - combined_hists.reverse_replays_ripple;
combined_hists.congruent_forward_minus_reverse_ripple = combined_hists.forward_congruent_replays_ripple - combined_hists.reverse_congruent_replays_ripple;
combined_hists.incongruent_forward_minus_reverse_ripple = combined_hists.forward_incongruent_replays_ripple - combined_hists.reverse_incongruent_replays_ripple;
combined_hists.replays_spike = combined_hists.reverse_replays_spike + combined_hists.forward_replays_spike;
combined_hists.replays_ripple = combined_hists.reverse_replays_ripple + combined_hists.forward_replays_ripple;

event_types = combined_hists.Properties.VariableNames;
% Convert hists from counts to rates

combined_hists_rates = table();
for i = 1:length(event_types)
    combined_hists_rates.(event_types{i}) = combined_hists.(event_types{i})./hist_bin_width;
end

w = setUp_gaussFilt([1 filter_length ],smoothing_sigma);
if smooth_rate_plots==1
    for i = 1:length(event_types)
        for j = 1:size(combined_hists_rates,1)
            rate_sub = combined_hists_rates.(event_types{i})(j,:);
            no_nan_inds = find(~isnan(rate_sub));
            rate_sub(isnan(rate_sub)) = [];
            combined_hists_rates.(event_types{i})(j,no_nan_inds) = conv(rate_sub,w,'same');
        end
    end
end
mean_combined_hists_rates = table();
sem_combined_hists_rates = table();
mean_combined_hists_rates_off = table();
mean_combined_hists_rates_on = table();
sem_combined_hists_rates_off = table();
sem_combined_hists_rates_on = table();
for i = 1:length(event_types)
    laser_off_data = combined_hists_rates.(event_types{i})(combined_table.laser_state==0,:);
    laser_on_data = combined_hists_rates.(event_types{i})(combined_table.laser_state==1,:);
    data = combined_hists_rates.(event_types{i});
    laser_off_n = sum(~isnan(laser_off_data),1);
    laser_on_n = sum(~isnan(laser_on_data),1);
    n = sum(~isnan(data),1);
    mean_combined_hists_rates_off.(event_types{i}) = nanmean(laser_off_data,1);
    mean_combined_hists_rates_on.(event_types{i}) = nanmean(laser_on_data,1);
    mean_combined_hists_rates.(event_types{i}) = nanmean(data,1);
    sem_combined_hists_rates_off.(event_types{i}) = nanstd(laser_off_data,1)./(sqrt(laser_off_n));
    sem_combined_hists_rates_on.(event_types{i}) = nanstd(laser_on_data,1)./sqrt(laser_on_n);
    sem_combined_hists_rates.(event_types{i}) = nanstd(data,1)./sqrt(n);    
end

num_trials_off = sum(~isnan(combined_hists_rates.sde(combined_table.laser_state==0,:)));
[min(num_trials_off(1:20)) max(num_trials_off(1:20))]
%%
mean_combined_table_rate_off = table();
mean_combined_table_rate_on = table();
sem_combined_table_rate_off = table();
sem_combined_table_rate_on = table();
sem_combined_table_rate = table();

event_types = combined_table_rate.Properties.VariableNames;
for i = 1:length(event_types)
    laser_off_data = combined_table_rate.(event_types{i})(combined_table.laser_state==0,:);
    laser_on_data = combined_table_rate.(event_types{i})(combined_table.laser_state==1,:);
    data = combined_table_rate.(event_types{i});  
    laser_off_n = sum(~isnan(laser_off_data),1);
    laser_on_n = sum(~isnan(laser_on_data),1);
    n=sum(~isnan(data),1);
    mean_combined_table_rate_off.(event_types{i}) = nanmean(laser_off_data,1);
    mean_combined_table_rate_on.(event_types{i}) = nanmean(laser_on_data,1);
    mean_combined_table_rate.(event_types{i}) = nanmean(data,1);
    sem_combined_table_rate_off.(event_types{i}) = nanstd(laser_off_data,1)./(sqrt(laser_off_n));
    sem_combined_table_rate_on.(event_types{i}) = nanstd(laser_on_data,1)./(sqrt(laser_on_n));
    sem_combined_table_rate.(event_types{i}) = nanstd(data,1)./sqrt(n);
end

%% For laser-OFF stopping periods that have both a forward and reverse event, which came first?

if align_to_stopping_period_start==1
    events_to_plot = {'forward_congruent_replays_spike','reverse_congruent_replays_spike'};
   %  events_to_plot = {'forward_congruent_replays_ripple','reverse_congruent_replays_ripple'};

    laser_off_forward = combined_table_first_event_times.(events_to_plot{1})(combined_table.laser_state==0,:);
    laser_off_reverse = combined_table_first_event_times.(events_to_plot{2})(combined_table.laser_state==0,:);
    good_stopping_periods = find(laser_off_forward < 10 & laser_off_reverse < 10);
    reverse_minus_forward_time = laser_off_reverse(good_stopping_periods)-laser_off_forward(good_stopping_periods);
    figure('Position',[600 600 100 100])
    [a,b] = hist(reverse_minus_forward_time,linspace(-8,8,16));
    h = bar(b,a);
    box off
    h.BarWidth = 1; h.FaceColor = [0 0 0];
    xline(0,'--k','LineWidth',2)
    axis square
    xlabel('Time (s)')
    ylabel('Stopping periods')
    set(gca,'FontSize',8)
    [p,~,z] = signrank(reverse_minus_forward_time)
    title(['p=' num2str(p,2)],'FontSize',8)
    set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
    saveas(gcf,[fullfile(final_fig_path,'time_of_first_reverse_minus_forward')],'jpeg')
    saveas(gcf,[fullfile(final_fig_path,'time_of_first_reverse_minus_forward')],'pdf')    

    median(reverse_minus_forward_time)
    quantile(reverse_minus_forward_time,0.25)
    quantile(reverse_minus_forward_time,0.5)
    quantile(reverse_minus_forward_time,0.75)
end

%% plot distribution of rates: (Fig 1E)
figure('Position',[803 408 100 100])
colors = [.4660 0.6740 0.1880; 0.4940 0.1840 0.5560];
a=combined_table_rate.forward_congruent_replays_spike(combined_table.laser_state==0);
b=combined_table_rate.reverse_congruent_replays_spike(combined_table.laser_state==0);

data = [mean(a); mean(b)];
err = [std(a)./sqrt(length(a)); std(b)./sqrt(length(b))];

b1 = bar(1,data(1)); hold on;
e1 = errorbar(1,data(1),err(1),'k','linestyle','none');
e1.CapSize = 4;
hold on
b2 = bar(2,data(2)); hold on;
e2 = errorbar(2,data(2),err(2),'k','linestyle','none');
e2.CapSize = 4;
b1.FaceColor = colors(1,:);
b1.EdgeColor = 'none';
b2.FaceColor = colors(2,:);
b2.EdgeColor = 'none';
set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
box off
saveas(gcf,fullfile(final_fig_path,'mean forward reverse rates'),'pdf')

%% Stats:
% for_minus_rev_spike_off_sig = signrank(combined_hists_rates.("forward_minus_reverse_spike")(combined_table.laser_state==0,:)); for_minus_rev_spike_off_sig(isnan(for_minus_rev_spike_off_sig)) = 0; 
% for_minus_rev_spike_on_sig = signrank(combined_hists_rates.("forward_minus_reverse_spike")(combined_table.laser_state==1,:)); for_minus_rev_spike_on_sig(isnan(for_minus_rev_spike_on_sig)) = 0; 
% congruent_for_minus_rev_spike_off_sig = signrank(combined_hists_rates.("congruent_forward_minus_reverse_spike")(combined_table.laser_state==0,:)); congruent_for_minus_rev_spike_off_sig(isnan(congruent_for_minus_rev_spike_off_sig)) = 0; 
% congruent_for_minus_rev_spike_on_sig = signrank(combined_hists_rates.("congruent_forward_minus_reverse_spike")(combined_table.laser_state==1,:)); congruent_for_minus_rev_spike_on_sig(isnan(congruent_for_minus_rev_spike_on_sig)) = 0; 
% for_minus_rev_ripple_off_sig = signrank(combined_hists_rates.("forward_minus_reverse_ripple")(combined_table.laser_state==0,:)); for_minus_rev_ripple_off_sig(isnan(for_minus_rev_ripple_off_sig)) = 0; 
% for_minus_rev_ripple_on_sig = signrank(combined_hists_rates.("forward_minus_reverse_ripple")(combined_table.laser_state==1,:)); for_minus_rev_ripple_on_sig(isnan(for_minus_rev_ripple_on_sig)) = 0; 
% congruent_for_minus_rev_ripple_off_sig = signrank(combined_hists_rates.("congruent_forward_minus_reverse_ripple")(combined_table.laser_state==0,:)); congruent_for_minus_rev_ripple_off_sig(isnan(congruent_for_minus_rev_ripple_off_sig)) = 0; 
% congruent_for_minus_rev_ripple_on_sig = signrank(combined_hists_rates.("congruent_forward_minus_reverse_ripple")(combined_table.laser_state==1,:)); congruent_for_minus_rev_ripple_on_sig(isnan(congruent_for_minus_rev_ripple_on_sig)) = 0; 

for_minus_rev_spike_off_sig = nan(size(combined_hists_rates.('forward_minus_reverse_spike'),2),1);
for bin = 1:size(combined_hists_rates.('forward_minus_reverse_spike'),2)
    [~,for_minus_rev_spike_off_sig(bin)] = signrank(combined_hists_rates.("forward_minus_reverse_spike")(combined_table.laser_state==0,bin)); 
end

for_minus_rev_spike_on_sig = nan(size(combined_hists_rates.('forward_minus_reverse_spike'),2),1);
for bin = 1:size(combined_hists_rates.('forward_minus_reverse_spike'),2)
    [~,for_minus_rev_spike_on_sig(bin)] = signrank(combined_hists_rates.("forward_minus_reverse_spike")(combined_table.laser_state==1,bin)); 
end

for_minus_rev_spike_sig = nan(size(combined_hists_rates.('forward_minus_reverse_spike'),2),1);
for bin = 1:size(combined_hists_rates.('forward_minus_reverse_spike'),2)
    [~,for_minus_rev_spike_sig(bin)] = signrank(combined_hists_rates.("forward_minus_reverse_spike")(:,bin)); 
end

congruent_for_minus_rev_spike_off_sig = nan(size(combined_hists_rates.('congruent_forward_minus_reverse_spike'),2),1);
for bin = 1:size(combined_hists_rates.('congruent_forward_minus_reverse_spike'),2)
    [~,congruent_for_minus_rev_spike_off_sig(bin)] = signrank(combined_hists_rates.("congruent_forward_minus_reverse_spike")(combined_table.laser_state==0,bin)); 
end

congruent_for_minus_rev_spike_on_sig = nan(size(combined_hists_rates.('congruent_forward_minus_reverse_spike'),2),1);
for bin = 1:size(combined_hists_rates.('congruent_forward_minus_reverse_spike'),2)
    [~,congruent_for_minus_rev_spike_on_sig(bin)] = signrank(combined_hists_rates.("congruent_forward_minus_reverse_spike")(combined_table.laser_state==1,bin)); 
end

congruent_for_minus_rev_spike_sig = nan(size(combined_hists_rates.('congruent_forward_minus_reverse_spike'),2),1);
for bin = 1:size(combined_hists_rates.('congruent_forward_minus_reverse_spike'),2)
    [~,congruent_for_minus_rev_spike_sig(bin)] = signrank(combined_hists_rates.("congruent_forward_minus_reverse_spike")(:,bin)); 
end


for_minus_rev_ripple_off_sig = nan(size(combined_hists_rates.('forward_minus_reverse_ripple'),2),1);
for bin = 1:size(combined_hists_rates.('forward_minus_reverse_ripple'),2)
    [~,for_minus_rev_ripple_off_sig(bin)] = signrank(combined_hists_rates.("forward_minus_reverse_ripple")(combined_table.laser_state==0,bin)); 
end

for_minus_rev_ripple_on_sig = nan(size(combined_hists_rates.('forward_minus_reverse_ripple'),2),1);
for bin = 1:size(combined_hists_rates.('forward_minus_reverse_ripple'),2)
    [~,for_minus_rev_ripple_on_sig(bin)] = signrank(combined_hists_rates.("forward_minus_reverse_ripple")(combined_table.laser_state==1,bin)); 
end

for_minus_rev_ripple_sig = nan(size(combined_hists_rates.('forward_minus_reverse_ripple'),2),1);
for bin = 1:size(combined_hists_rates.('forward_minus_reverse_ripple'),2)
    [~,for_minus_rev_ripple_sig(bin)] = signrank(combined_hists_rates.("forward_minus_reverse_ripple")(:,bin)); 
end

congruent_for_minus_rev_ripple_off_sig = nan(size(combined_hists_rates.('congruent_forward_minus_reverse_ripple'),2),1);
for bin = 1:size(combined_hists_rates.('congruent_forward_minus_reverse_ripple'),2)
    [~,congruent_for_minus_rev_ripple_off_sig(bin)] = signrank(combined_hists_rates.("congruent_forward_minus_reverse_ripple")(combined_table.laser_state==0,bin)); 
end

congruent_for_minus_rev_ripple_on_sig = nan(size(combined_hists_rates.('congruent_forward_minus_reverse_ripple'),2),1);
for bin = 1:size(combined_hists_rates.('congruent_forward_minus_reverse_ripple'),2)
    [~,congruent_for_minus_rev_ripple_on_sig(bin)] = signrank(combined_hists_rates.("congruent_forward_minus_reverse_ripple")(combined_table.laser_state==1,bin)); 
end

congruent_for_minus_rev_ripple_sig = nan(size(combined_hists_rates.('congruent_forward_minus_reverse_ripple'),2),1);
for bin = 1:size(combined_hists_rates.('congruent_forward_minus_reverse_ripple'),2)
    [~,congruent_for_minus_rev_ripple_sig(bin)] = signrank(combined_hists_rates.("congruent_forward_minus_reverse_ripple")(:,bin)); 
end

plot_fig_2_CD; % Generates Figure 2C-D
plot_supp_fig_S8_FG % Generates Figures S8 C-D


