


% 1) Add in stopping period duration threshold option so that you can
% ensure that none of the period you are looking at here corresponds to the
% rat's drinking time.

% 2) Change the definition of exit time to be the last time that the rat
% was in the reward zone moving under the speed threshold.



%close all;
clear
windows=0;
plot_by_session = 0;
plot_by_rat = 0;
hand_clustered_only = 0;
rats = [1 2 3 6];
run_stats = 0;

must_include_all_flags = 1; % 1: requires that all flags in Flags_to_include are present. 0: requires that at least one flag in Flags_to_include are present.
novel_cutoff = 1; % of days that will be considered novel

novel_range = [10.2 10.5]

flags.flags_to_include = [11];
flags.flags_to_exclude = [0];

events_to_compare = [0 1]; % 0 = laser was off, 1 = laser was on
candidate_events_file_name = 'decoder_candidateEvents.mat';
candidate_events_to_plot = 'spike_filtered';
min_num_replays_per_condition =0; % If this is set to 1, the number below must be met for both laser ON and laser OFF epochs. otherwise, the number is for all replays in the session.
min_num_replays_per_session = 0; % this will remove sessions with fewer than min_num_replays_per_session in both laser on and laser off epochs from analysis
min_num_passes_per_session = 0;
min_number_replays_in_limited_stopping_period =  0;
Rat_Names = {'Clover','Bo','MEC1','Bolt','Dash','CM1'};

bin_width = 0.5;
smoothing_sigma = 3 %bins;
filter_length = smoothing_sigma*6;
if(mod(filter_length,2)==0) % iseven
    filter_length = filter_length+1;
end

min_stopping_period_duration = 0; % s
max_stopping_period_duration = inf;

time_to_plot = 10; % sec
num_bins_to_plot = time_to_plot*(1/bin_width)
smooth_rate_plots = 1;

hist_edges = 0-time_to_plot:bin_width:0;


plot_high_reverse_sessions_only = 0;
plot_low_reverse_sessions_only = 0;

% if strcmp(candidate_events_to_plot,'ripple')
% event_choice = 1;
% elseif strcmp(candidate_events_to_plot,'spike')
% event_choice = 2;
% elseif strcmp(candidate_events_to_plot,'laser')
% event_choice = 3;
% end

% each entry will be one candidate event
data = [];
rat_label = [];
day_label = [];
saline_label = [];
cno_label = [];
unique_session_id = [];
segment_label = [];
group_type_label = [];
novel_label = [];
reward_stim_label = [];
run_stim_label = [];

combined_congruent_replay_hist = [];
combined_incongruent_replay_hist = [];

combined_reverse_replay_hist = [];
combined_forward_replay_hist = [];

combined_reverse_congruent_replay_hist = [];
combined_forward_congruent_replay_hist = [];

combined_reverse_incongruent_replay_hist = [];
combined_forward_incongruent_replay_hist = [];

combined_sde_hist = [];
combined_replay_hist = [];

session_id = 0;

% Where to save the data:
if windows==1
    fig_path_prefix = 'D:/Dropbox/Foster Lab/Data/Replay summary/Replay_comparisons';
else
    fig_path_prefix = '/home/caitlin/Insync/caitlinmallory@berkeley.edu/Google Drive/Processed_Data/Replay_summary';
end

if hand_clustered_only == 1
    fig_path = fullfile(fig_path_prefix,'hand_clustered');
else
    fig_path = fullfile(fig_path_prefix, 'best_decoding');
end

if strcmp(candidate_events_to_plot,'spike')
    fig_path = fullfile(fig_path,'spike');
elseif strcmp(candidate_events_to_plot,'filtered')
    fig_path = fullfile(fig_path,'filtered');
elseif strcmp(candidate_events_to_plot,'ripple')
    fig_path = fullfile(fig_path,'ripple');
elseif strcmp(candidate_events_to_plot,'spike_filtered')
    fig_path = fullfile(fig_path,'spike_filtered');
else
    keyboard
end

% make directory to save figs if it doesn't exist
if ~exist(fig_path,'dir')
    mkdir(fig_path)
end

for rat = rats
    if rat == 1 %Clover

        if windows == 1
            directory = 'G:/My Drive/Processed_Data/Clover/linear_track';
        else
            directory = '/home/caitlin/Insync/caitlinmallory@berkeley.edu/Google Drive/Processed_Data/Clover/linear_track';
        end

        if hand_clustered_only == 1
            % only use hand clustered sessions:
            dayFiles = {'20210427/hand_clustered','20210428/hand_clustered','20210429/hand_clustered','20210430/hand_clustered','20210501/hand_clustered','20210507/hand_clustered','20210508/hand_clustered','20210509/hand_clustered','20210512/hand_clustered','20210517/hand_clustered','20210519/hand_clustered'};
        else
            % use the sessions with the best decoding:
            dayFiles = {'20210427/hand_clustered','20210428/hand_clustered','20210429/kilosort_clustered','20210430/kilosort_clustered','20210501/hand_clustered','20210507/kilosort_clustered','20210508/hand_clustered','20210509/kilosort_clustered','20210512/kilosort_clustered','20210517/hand_clustered','20210519/kilosort_clustered'};
        end
    end
    if rat == 2 %Bo
        if windows == 1
            directory = 'G:/My Drive/Processed_Data/Bo/linear_track';
        else
            directory = '/home/caitlin/Insync/caitlinmallory@berkeley.edu/Google Drive/Processed_Data/Bo/linear_track';
        end
        if hand_clustered_only == 1
            % hand clustered sessions
            dayFiles =  {'20210509/hand_clustered','20210511/hand_clustered','20210512/hand_clustered','20210513/hand_clustered','20210517/hand_clustered'};
        else
            %best decoding
            dayFiles = {'20210511/kilosort_clustered','20210512/kilosort_clustered','20210513/kilosort_clustered','20210517/kilosort_clustered'};

        end

    end
    if rat == 3 % MEC1
        if windows == 1
            directory = 'G:/My Drive/Processed_Data/MEC1';
        else
            %             directory = '/media/caitlin/Caitlin_Drive_4/MEC1_data/';
            directory = '/home/caitlin/Insync/caitlinmallory@berkeley.edu/Google Drive/Processed_Data/MEC1';
        end

        if hand_clustered_only == 1
            % hand clustered sessions
            dayFiles = {'20200926/original_clustering','20200927/hand_clustered','20200928/hand_clustered_clusters'};
        else
            %best decoding:
            % dayFiles = {'20200926/kilosort_clustered','20200927/hand_clustered'};
            %'20200927/kilosort_clustered','20200928/kilosort_clustered;
            dayFiles = {'20200926/kilosort_clustered','20200927/hand_clustered2'};
        end
    end
    if rat == 4
        if windows == 1
            % directory = 'L:/';
            directory = 'G:/My Drive/Processed_Data/Bolt/linear_track';
        else
            directory = '/home/caitlin/Insync/caitlinmallory@berkeley.edu/Google Drive/Processed_Data/Bolt/linear_track/';
        end
              dayFiles = {'20220714/msort_clustered','20220715/msort_clustered','20220718_1/msort_clustered','20220718_2/msort_clustered','20220722/msort_clustered_full','20220729/msort_clustered'};
        %
        %     dayFiles = {
        %
        %     '20220531',...
        %      '20220607_linear_track/mountainsort_clustered_full',...
        %
        %             }

        %         '20220608_162649_novel_lin_track_1/processed/msort_clustered_full', ...
        %             '20220608_172158_bolt_novel_lin_track2_cno/processed/msort_clustered_full',...
        %             ...
        %             '20220531',...
        %             '20220601/mountainsort_clustered_full',...
        %             '20220603/mountainsort_clustered_full',...
        %             '20220607_linear_track/mountainsort_clustered_full',...
        %             '20220608_162649_novel_lin_track_1/processed/msort_clustered_full', ...
        %             '20220608_172158_bolt_novel_lin_track2_cno/processed/msort_clustered_full',...
        %             '20220609_novel_lintrack_3/msort_clustered_full',...
        %             '20220614_102255_novel_lintrack_5_saline/processed/msort_clustered_full',...
        %              '20220614_105754_novel_lintrack_6_saline/processed/msort_clustered_full', ...
        %              '20220613_novel_lintrack_4_cno/msort_clustered_full'}
      

    elseif rat == 5

        if windows == 1
            directory = 'G:/My Drive/Processed_Data/Dash';
        else
            directory = '/home/caitlin/Insync/caitlinmallory@berkeley.edu/Google Drive/Processed_Data/Dash';

        end
        dayFiles = {'20221025_113913'; '20221025_170200'};

    elseif rat == 6
        if windows == 1
            % directory = 'L:/';
            directory = 'G:/My Drive/Processed_Data/CM1';
        else
            directory = '/home/caitlin/Insync/caitlinmallory@berkeley.edu/Google Drive/Processed_Data/CM1';
        end
            dayFiles = ...
                {'20230217','20230219_01','20230219','20230220','20230221','20230221_2_bad','20230221_3','20230222_1','20230224_2','20230228','20230301_01','20230302_01','20230302_02','20230303_01','20230303_02','20230305_01','20230305_02'};

 

             % It's unclear whether 20230224_2 should count as novel. there was
             % a 20 minutes session earlier but the data was lost. the
             % position of the track wasn't very different from the
             % previous.



    end

    for day = 1:length(dayFiles)
        display([Rat_Names{rat},' Day ', dayFiles{day}])
        cd(fullfile(directory,dayFiles{day}))

        load Experiment_Information

        sessions_that_meet_criterion_day = check_for_relevant_sessions(flags,1,Experiment_Information,must_include_all_flags,1,novel_range);

    
            flags_to_include = flags.flags_to_include;
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


               load(['replay_by_stopping_period_aligned_to_departure_' candidate_events_to_plot '_segment_' num2str(sub_session_num) '.mat']);



                properties =  t2.Properties.VariableNames;
                t2 = table2array(t2);


                rat_label = [rat_label; repmat(rat,size(t2,1),1)];
                day_label = [day_label; repmat(day,size(t2,1),1)];
                unique_session_id = [unique_session_id; repmat(session_id,size(t2,1),1)];
                segment_label = [segment_label; repmat(sub_session_num,size(t2,1),1)];
         
                novel_label = [novel_label; repmat(novel,size(t2,1),1)];
                run_stim_label = [run_stim_label; repmat(run_stim,size(t2,1),1)];
                reward_stim_label = [reward_stim_label; repmat(reward_stim,size(t2,1),1)];
                saline_label = [saline_label; repmat(saline,size(t2,1),1)];
                cno_label = [cno_label; repmat(cno,size(t2,1),1)];

                data = [data; t2];
                combined_reverse_replay_hist = [combined_reverse_replay_hist; reverse_replay_hist];
                combined_forward_replay_hist = [combined_forward_replay_hist; forward_replay_hist];


                combined_reverse_congruent_replay_hist = [combined_reverse_congruent_replay_hist; reverse_congruent_replay_hist];
                combined_forward_congruent_replay_hist = [combined_forward_congruent_replay_hist; forward_congruent_replay_hist];

                combined_reverse_incongruent_replay_hist = [combined_reverse_incongruent_replay_hist; reverse_incongruent_replay_hist];
                combined_forward_incongruent_replay_hist = [combined_forward_incongruent_replay_hist; forward_incongruent_replay_hist];

                combined_congruent_replay_hist = [combined_congruent_replay_hist; congruent_replay_hist];
                combined_incongruent_replay_hist = [combined_incongruent_replay_hist; incongruent_replay_hist];
                    
   
                combined_sde_hist = [combined_sde_hist; sde_hist];
                combined_replay_hist = [combined_replay_hist; reverse_replay_hist+forward_replay_hist];
            end

    end
    cd ..
end

%%

% put all the data into a table
% rat = rat_label
% day = day_label;
% lap - data(:,1);
% laser_state = data(:,2);
% duration = data(:,3)
tbl = table();
tbl.rat = rat_label;
tbl.day = day_label;
tbl.unique_session_id = unique_session_id;
tbl.saline = saline_label;
tbl.cno = cno_label;
tbl.lap = data(:,1);


% Real:
tbl.laser_state = data(:,2);

% saline sessions will be labeled as 'laser off' and cno sessions 'laser
% on'
tbl.laser_state(tbl.saline == 1) = 0;
tbl.laser_state(tbl.cno == 1) = 1;



tbl.laser_state_stopping_period = data(:,3);
tbl.duration = data(:,4);
tbl.reverse_replay = data(:,5);
tbl.forward_replay = data(:,6);
tbl.congruent_replay = data(:,7);
tbl.incongruent_replay = data(:,8);
tbl.reverse_congruent_replay = data(:,9);
tbl.reverse_incongruent_replay = data(:,10);
tbl.forward_congruent_replay = data(:,11);
tbl.forward_incongruent_replay = data(:,12);
tbl.sde_count = data(:,13);
tbl.novel = novel_label;
tbl.reward_stim = reward_stim_label;
tbl.run_stim = run_stim_label;



tbl(tbl.duration == 0,:) = [];

tbl_copy = tbl;

% remove sessions where there weren't enough replay events
unique_sessions = unique(tbl.unique_session_id);
% remove sessions with too few laps


num_laps = nan(length(unique_sessions),1);
num_replays_session = nan(length(unique_sessions),1);
num_good_replays_session_laser_on = nan(length(unique_sessions),1);
num_good_replays_session_laser_off = nan(length(unique_sessions),1);
for i = 1:length(unique_sessions)
    tbl_sub = tbl(tbl.unique_session_id == unique_sessions(i),:);
    num_replays_session(i) = sum(tbl_sub.forward_replay) + sum(tbl_sub.reverse_replay);
    num_good_replays_session_laser_on(i) = sum(tbl_sub.forward_replay(tbl_sub.laser_state ==1)) + sum(tbl_sub.reverse_replay(tbl_sub.laser_state==1));
    num_good_replays_session_laser_off(i) = sum(tbl_sub.forward_replay(tbl_sub.laser_state ==0)) + sum(tbl_sub.reverse_replay(tbl_sub.laser_state==0));
    num_laps(i) = height(tbl_sub);
end

if min_num_replays_per_condition == 1
    bad_session_inds_not_enough_replays = [find(num_good_replays_session_laser_on < min_num_replays_per_session); ...
        find(num_good_replays_session_laser_off < min_num_replays_per_session)]; 
else
    bad_session_inds_not_enough_replays = find(num_replays_session < min_num_replays_per_session);
end

bad_session_inds_not_enough_laps = find(num_laps < min_num_passes_per_session);

bad_session_inds_reverse_replay_out_of_range = [];
if plot_low_reverse_sessions_only == 1
    bad_session_inds_reverse_replay_out_of_range =  find(prop_reverse_laser_off >= 0.4);
end
if plot_high_reverse_sessions_only == 1
    bad_session_inds_reverse_replay_out_of_range = find(prop_reverse_laser_off <= 0.6);
end

bad_session_inds = unique([bad_session_inds_not_enough_replays; bad_session_inds_not_enough_laps; bad_session_inds_reverse_replay_out_of_range]);

%
% bad_inds = unique([find(tbl.duration ==0); find(tbl.time_since_reward_entry<min_stopping_period_duration); find(tbl.time_since_reward_entry>max_stopping_period_duration); find(ismember(tbl.unique_session_id,unique_sessions(bad_session_inds)))]);

bad_inds = unique([find(tbl.duration ==0); find(tbl.duration<min_stopping_period_duration); find(tbl.duration>max_stopping_period_duration); find(ismember(tbl.unique_session_id,unique_sessions(bad_session_inds)))]);


tbl(bad_inds,:) = [];


tbl2 = table();
tbl2.rat = [tbl.rat; tbl.rat];
tbl2.day = [tbl.day; tbl.day];
tbl2.unique_session_id = [tbl.unique_session_id; tbl.unique_session_id];
tbl2.lap = [tbl.lap; tbl.lap];
tbl2.laser_state = [tbl.laser_state; tbl.laser_state];
tbl2.laser_state_stopping_period = [tbl.laser_state_stopping_period; tbl.laser_state_stopping_period];
tbl2.duration = [tbl.duration; tbl.duration];
tbl2.replay = [tbl.reverse_replay; tbl.forward_replay];
tbl2.reverse = [ones(size(tbl.reverse_replay)); zeros(size(tbl.forward_replay))];
tbl2.replay_congruent =  [tbl.congruent_replay; tbl.incongruent_replay];
tbl2.reverse_congruent = [ones(size(tbl.congruent_replay)); zeros(size(tbl.incongruent_replay))];


duration2 = [tbl2.duration; tbl2.duration];

%%
% Look at session averages







%%
if run_stats == 1
    glme = fitglme(tbl,...
        'replay ~ 1 + laser_state + (1|rat:day)',...
        'Distribution','Poisson','Offset',log(tbl.duration))

    glme = fitglme(tbl, ...
        'sde_count ~ 1 + laser_state + (1|rat:day)', ...
        'Distribution','Poisson','Offset',log(tbl.duration))


    glme = fitglme(tbl, ...
        'congruent_replay ~ 1 + laser_state + (1|rat:day)', ...
        'Distribution','Poisson','Offset',log(tbl.duration))

    glme = fitglme(tbl, ...
        'incongruent_replay ~ 1 + laser_state +  (1|rat:day)', ...
        'Distribution','Poisson','Offset',log(tbl.duration))

    glme = fitglme(tbl, ...
        'reverse_replay ~ 1 + laser_state + (1|rat:day)', ...
        'Distribution','Poisson','Offset',log(tbl.duration))

    glme = fitglme(tbl, ...
        'reverse_congruent_replay ~ 1 + laser_state + (1|rat:day)', ...
        'Distribution','Poisson','Offset',log(tbl.duration))

    glme = fitglme(tbl, ...
        'forward_replay ~ 1 + laser_state + (1|rat:day)', ...
        'Distribution','Poisson','Offset',log(tbl.duration))

    glme = fitglme(tbl, ...
        'forward_congruent_replay ~ 1 + laser_state + (1|rat:day)', ...
        'Distribution','Poisson','Offset',log(tbl.duration))





    glme = fitglme(tbl2, ...
        'replay ~ 1 + laser_state + reverse + laser_state*reverse + (1|rat:day)',...
        'Distribution','Poisson','Offset',log(tbl2.duration))

    glme = fitglme(tbl2, ...
        'replay_congruent ~ 1 + laser_state + reverse_congruent + laser_state*reverse_congruent + (1|rat:day)',...
        'Distribution','Poisson','Offset',log(tbl2.duration))
    %
    % glme = fitglme(tbl, ...
    % %         'reverse_replay ~ 1 + laser_state + reward_stim + reward_stim*laser_state + (1|rat) + (1|day:rat)', ...
    % %     'Distribution','Poisson','Offset',log(duration))
end
%%

sde_rate = tbl.sde_count./tbl.duration;

sde_0 = sde_rate(tbl.laser_state == 0);
sde_0_mean = mean(sde_0);
sde_0_sem = std(sde_0)/(sqrt(length(sde_0)));

sde_1 = sde_rate(tbl.laser_state == 1);
sde_1_mean = mean(sde_1);
sde_1_sem = std(sde_1)/(sqrt(length(sde_1)));


y = [sde_0_mean sde_1_mean];
err = [sde_0_sem sde_1_sem];

figure(); clf;
x = [1:2];
b = bar(x,y)
b(1).FaceColor = [0.5 0.5 0.5]
hold on
er = errorbar(x,y,err,err);
er.Color = [0 0 0];
er.LineStyle = 'none';

% Set Axis properties
set(gca,'xticklabel',{'Laser off'; 'Laser on'});
set(gca,'FontSize', 16)
% ylim([0 0.8])

if strcmp(candidate_events_to_plot, 'ripple')
    ylabel('Ripple rate (Hz)')
elseif strcmp(candidate_events_to_plot, 'spike')
    ylabel('Spike density event rate (Hz)')
elseif strcmp(candidate_events_to_plot, 'filtered')
    ylabel('Replay attempt event rate (Hz)')
elseif strcmp(candidate_events_to_plot, 'spike_filtered')
    ylabel('filtered spike density event rate (Hz)')
else
    keyboard
end

set(gcf,'Position',[20 20 300 400])

saveas(gcf,fullfile(fig_path,['sde rate rats' num2str(rats)]),'jpg')


%%

reverse_replay_rate = tbl.reverse_replay./tbl.duration;
forward_replay_rate = tbl.forward_replay./tbl.duration;
tbl.reverse_replay_rate = reverse_replay_rate;
tbl.forward_replay_rate = forward_replay_rate;

r_0 = reverse_replay_rate(tbl.laser_state == 0);
r_0_mean = mean(r_0);
r_0_sem = std(r_0)/(sqrt(length(r_0)));
r_0_95_conf = 1.96*r_0_sem;

f_0 = forward_replay_rate(tbl.laser_state == 0);
f_0_mean = mean(f_0);
f_0_sem = std(f_0)/(sqrt(length(f_0)));
f_0_95_conf = 1.96*f_0_sem;


rl_1 = reverse_replay_rate(tbl.laser_state == 1);
rl_1_mean = mean(rl_1);
rl_1_sem = std(rl_1)/(sqrt(length(rl_1)));
rl_1_95_conf = 1.96*rl_1_sem;

fl_1 = forward_replay_rate(tbl.laser_state == 1);
fl_1_mean = mean(fl_1);
fl_1_sem = std(fl_1)/(sqrt(length(fl_1)));
fl_1_95_conf = 1.96*fl_1_sem;

[p,h] = ranksum(r_0,rl_1);
[p,h] = ranksum(f_0,fl_1);

y = [r_0_mean f_0_mean; rl_1_mean fl_1_mean];
err = [r_0_sem f_0_sem; rl_1_sem fl_1_sem;];

figure(); clf;
hb = bar(y); % get the bar handles
hb(1).FaceColor = [0.4940 0.1840 0.5560];
hb(2).FaceColor = [.4660 0.6740 0.1880];

hold on;
for k = 1:size(y,2)
    % get x positions per group
    xpos = hb(k).XData + hb(k).XOffset;
    % draw errorbar
    errorbar(xpos, y(:,k), err(:,k), 'LineStyle', 'none', ...
        'Color', 'k', 'LineWidth', 1);
end

% Set Axis properties
set(gca,'xticklabel',{'Laser off'; 'Laser on'});
set(gca,'FontSize', 16)

if ismember(4,rats)
    ylim([0 0.06])
    yticks(linspace(0,0.06,5))
else
    ylim([0 0.02])
    yticks(linspace(0, 0.02, 5))
end
set(gcf,'Position', [10 10 300 400])
ylabel('Long replay rate (events/s)')
legend({'Reverse','Forward'})
fig_title = ['Replay rate, n = ' num2str(length(r_0)) ' laser off, ' num2str(length(rl_1)) ' laser on stopping periods'];
title(fig_title)

% ylim([0 0.015])
% yticks(linspace(0,0.015,4))
saveas(gcf,fullfile(fig_path,['replay rate by stopping period rats' num2str(rats)]),'jpg')


dv  = [r_0 f_0 zeros(length([r_0 f_0]),1);rl_1 fl_1 ones(length([rl_1 fl_1]),1)];
dv(any(isnan(dv), 2), :) = [];
dv = array2table(dv);
% remove cells where there was no data for any of the conditions (i.e., no
% replays of that type)

% codes for column names:
dv.Properties.VariableNames = {'laserOff', 'laserON', 'direction'};
dv.direction = categorical(dv.direction);

% create the within-subjects design
withinDesign = table([0 1]','VariableNames',{'laser_state'});
withinDesign.laser_state = categorical(withinDesign.laser_state);

% create repeated measures model
rm = fitrm(dv, 'laserOff-laserON~direction', 'WithinDesign', withinDesign);

% perform anova (remove semicolon to view table generated by ranova)
AT = ranova(rm, 'WithinModel', 'laser_state');

% output an improved version of the anova table
disp(anovaTable(AT, 'Value'));

% BAR PLOT (MEAN and SEM)
% fig = two_by_two_barplot(data, group_colors, figure_xticklabels, figure_legend, figure_stats)
AT = table2array(AT);
repeated_measure_pval = AT(4,5);
group_effect_pval = AT(2,5);
interaction_pval = AT(5,5);
stats = ['rm p= ' num2str(repeated_measure_pval,2) ', group p=' num2str(group_effect_pval,2) ', interaction p=' num2str(interaction_pval,2)];
if ismember(rats,4)
    fig_ylim = [0 0.05];
else
    fig_ylim = [0 0.02];
end
set(gcf,'Position', [10 10 300 400])
two_by_two_barplot([{r_0} {f_0} {rl_1} {fl_1}], [1 0 0; 0 0 1], {'laser off'; 'laser on'}, fig_ylim, fig_title, {'reverse' 'forward'}, stats, '');
% figure_title = fullfile(figure_save_path,[fig_title num2str(rats)]);
% saveas(fig2,figure_title,'jpg')


%%

replay_rate = (tbl.reverse_replay + tbl.forward_replay)./tbl.duration;


replay_0 = replay_rate(tbl.laser_state == 0);
replay_0_mean = mean(replay_0);
replay_0_sem = std(replay_0)/(sqrt(length(replay_0)));

replay_1 = replay_rate(tbl.laser_state == 1);
replay_1_mean = mean(replay_1);
replay_1_sem = std(replay_1)/(sqrt(length(replay_1)));


[p,h] = ranksum(replay_0,replay_1);

y = [replay_0_mean replay_1_mean];
err = [replay_0_sem replay_1_sem];

figure(); clf;
x = [1:2];
b = bar(x,y)
b(1).FaceColor = [0.5 0.5 0.5]
hold on
er = errorbar(x,y,err,err);
er.Color = [0 0 0];
er.LineStyle = 'none';

% Set Axis properties
set(gca,'xticklabel',{'Laser off'; 'Laser on'});
set(gca,'FontSize', 16)

ylabel('Replay rate (Hz)')


set(gcf,'Position',[20 20 300 400])

saveas(gcf,fullfile(fig_path,['replay rate' num2str(rats)]),'jpg')






%% plot the rates of long reverse and long forward congruent replays
tbl.reverse_congruent_replay_rate = tbl.reverse_congruent_replay./tbl.duration;
tbl.forward_congruent_replay_rate = tbl.forward_congruent_replay./tbl.duration;

r_0 = tbl.reverse_congruent_replay_rate(tbl.laser_state == 0);
r_0_mean = mean(r_0);
r_0_sem = std(r_0)/(sqrt(length(r_0)));
r_0_95_conf = [1.96*r_0_sem];

r_1 = tbl.reverse_congruent_replay_rate(tbl.laser_state == 1);
r_1_mean = mean(r_1);
r_1_sem = std(r_1)/(sqrt(length(r_1)));
r_1_95_conf = [1.96*r_1_sem];

f_0 = tbl.forward_congruent_replay_rate(tbl.laser_state == 0);
f_0_mean = mean(f_0);
f_0_sem = std(f_0)/(sqrt(length(f_0)));
f_0_95_conf = [1.96*f_0_sem];


f_1 = tbl.forward_congruent_replay_rate(tbl.laser_state == 1);
f_1_mean = mean(f_1);
f_1_sem = std(f_1)/(sqrt(length(f_1)));
f_1_95_conf = [1.96*f_1_sem];

[p,h] = ranksum(r_0,r_1);
[p,h] = ranksum(f_0,f_1);

y = [r_0_mean f_0_mean; r_1_mean f_1_mean];
err = [r_0_sem f_0_sem; r_1_sem f_1_sem;];

figure(); clf;
hb = bar(y); % get the bar handles
hb(1).FaceColor = [0.4940 0.1840 0.5560];
hb(2).FaceColor = [.4660 0.6740 0.1880];
hold on;
for k = 1:size(y,2)
    % get x positions per group
    xpos = hb(k).XData + hb(k).XOffset;
    % draw errorbar
    errorbar(xpos, y(:,k), err(:,k), 'LineStyle', 'none', ...
        'Color', 'k', 'LineWidth', 1);
end

% Set Axis properties
set(gca,'xticklabel',{'Laser off'; 'Laser on'});
set(gca,'FontSize', 16)
ylabel('Replay rate (events/s)')
legend({'Reverse','Forward'})
fig_title = ['Congruent replay rate, n = ' num2str(length(r_0)) ' laser off, ' num2str(length(rl_1)) ' laser on stopping periods'];
title(fig_title)

saveas(gcf,fullfile(fig_path,['replay rate rats' num2str(rats)]),'jpg')

%%
if plot_by_session==1

% Average across session:
good_sessions = unique(tbl.unique_session_id);
r_0_session = nan(length(unique(tbl.unique_session_id)),1);
r_1_session = nan(length(unique(tbl.unique_session_id)),1);
f_0_session = nan(length(unique(tbl.unique_session_id)),1);
f_1_session = nan(length(unique(tbl.unique_session_id)),1);
r_0_session_wide = nan(length(unique(tbl.unique_session_id)),1);
r_1_session_wide = nan(length(unique(tbl.unique_session_id)),1);
f_0_session_wide = nan(length(unique(tbl.unique_session_id)),1);
f_1_session_wide = nan(length(unique(tbl.unique_session_id)),1);

r_0_session_wide_counts = nan(length(unique(tbl.unique_session_id)),1);
r_1_session_wide_counts = nan(length(unique(tbl.unique_session_id)),1);
f_0_session_wide_counts = nan(length(unique(tbl.unique_session_id)),1);
f_1_session_wide_counts = nan(length(unique(tbl.unique_session_id)),1);

duration_0_session_wide = nan(length(unique(tbl.unique_session_id)),1);
duration_1_session_wide = nan(length(unique(tbl.unique_session_id)),1);


for i = 1:length(good_sessions)
    r_0_session(i) = mean(reverse_replay_rate(tbl.laser_state == 0 & tbl.unique_session_id == good_sessions (i)));
    r_1_session(i) = mean(reverse_replay_rate(tbl.laser_state == 1 & tbl.unique_session_id == good_sessions (i)));
    f_0_session(i) = mean(forward_replay_rate(tbl.laser_state == 0 & tbl.unique_session_id == good_sessions (i)));
    f_1_session(i) = mean(forward_replay_rate(tbl.laser_state == 1 & tbl.unique_session_id == good_sessions (i)));

    duration_0_session_wide(i) = sum(tbl.duration(tbl.laser_state == 0 & tbl.unique_session_id == good_sessions (i)));
    duration_1_session_wide(i) = sum(tbl.duration(tbl.laser_state == 1 & tbl.unique_session_id == good_sessions (i)));


    r_0_session_wide_counts(i) = sum(tbl.reverse_replay(tbl.laser_state == 0 & tbl.unique_session_id == good_sessions (i)));
    r_1_session_wide_counts(i) = sum(tbl.reverse_replay(tbl.laser_state == 1 & tbl.unique_session_id == good_sessions (i)));
    f_0_session_wide_counts(i) = sum(tbl.forward_replay(tbl.laser_state == 0 & tbl.unique_session_id == good_sessions (i)));
    f_1_session_wide_counts(i) = sum(tbl.forward_replay(tbl.laser_state == 1 & tbl.unique_session_id == good_sessions (i)));

    r_0_session_wide(i) = sum(tbl.reverse_replay(tbl.laser_state == 0 & tbl.unique_session_id == good_sessions (i)))/duration_0_session_wide(i) ;
    r_1_session_wide(i) = sum(tbl.reverse_replay(tbl.laser_state == 1 & tbl.unique_session_id == good_sessions (i)))/duration_1_session_wide(i) ;
    f_0_session_wide(i) = sum(tbl.forward_replay(tbl.laser_state == 0 & tbl.unique_session_id == good_sessions (i)))/duration_0_session_wide(i) ;
    f_1_session_wide(i) = sum(tbl.forward_replay(tbl.laser_state == 1 & tbl.unique_session_id == good_sessions (i)))/duration_1_session_wide(i) ;
end

r_0_session_mean = mean(r_0_session);
r_1_session_mean = mean(r_1_session);
f_0_session_mean = mean(f_0_session);
f_1_session_mean = mean(f_1_session);
r_0_session_sem = std(r_0_session)/sqrt(length(r_0_session));
r_1_session_sem = std(r_1_session)/sqrt(length(r_1_session));
f_0_session_sem = std(f_0_session)/sqrt(length(f_0_session));
f_1_session_sem = std(f_1_session)/sqrt(length(f_1_session));


[p,h] = signrank(r_0_session,r_1_session)
[p,h] = signrank(f_0_session,f_1_session)

y = [r_0_session_mean f_0_session_mean; r_1_session_mean f_1_session_mean];
err = [r_0_session_sem f_0_session_sem; r_1_session_sem f_1_session_sem];

figure(); clf;
hb = bar(y); % get the bar handles
hb(1).FaceColor = [0.4940 0.1840 0.5560];
hb(2).FaceColor = [.4660 0.6740 0.1880];
hold on;
for k = 1:size(y,2)
    % get x positions per group
    xpos = hb(k).XData + hb(k).XOffset;
    % draw errorbar
    errorbar(xpos, y(:,k), err(:,k), 'LineStyle', 'none', ...
        'Color', 'k', 'LineWidth', 1);
end

% Set Axis properties
set(gca,'xticklabel',{'Laser off'; 'Laser on'});
set(gca,'FontSize', 16)
ylabel('Replay rate (events/s)')
legend({'Reverse','Forward'})


fig_title = ['Replay rate, n = ' num2str(length(r_0_session)) ' sessions'];
title(fig_title)
saveas(gcf,fullfile(fig_path,['Replay rate by session rats' num2str(rats)]),'jpg')

r_0_session_wide_mean = nanmean(r_0_session_wide);
r_1_session_wide_mean = nanmean(r_1_session_wide);
f_0_session_wide_mean = nanmean(f_0_session_wide);
f_1_session_wide_mean = nanmean(f_1_session_wide);
r_0_session_wide_sem = nanstd(r_0_session_wide)/sqrt(length(r_0_session_wide));
r_1_session_wide_sem = nanstd(r_1_session_wide)/sqrt(length(r_1_session_wide));
f_0_session_wide_sem = nanstd(f_0_session_wide)/sqrt(length(f_0_session_wide));
f_1_session_wide_sem = nanstd(f_1_session_wide)/sqrt(length(f_1_session_wide));


[p,h] = signrank(r_0_session_wide,r_1_session_wide)
[p,h] = signrank(f_0_session_wide,f_1_session_wide)

y = [r_0_session_wide_mean f_0_session_wide_mean; r_1_session_wide_mean f_1_session_wide_mean];
err = [r_0_session_wide_sem f_0_session_wide_sem; r_1_session_wide_sem f_1_session_wide_sem];

figure(); clf;
hb = bar(y); % get the bar handles
hb(1).FaceColor = [0.4940 0.1840 0.5560];
hb(2).FaceColor = [.4660 0.6740 0.1880];
hold on;
for k = 1:size(y,2)
    % get x positions per group
    xpos = hb(k).XData + hb(k).XOffset;
    % draw errorbar
    errorbar(xpos, y(:,k), err(:,k), 'LineStyle', 'none', ...
        'Color', 'k', 'LineWidth', 1);
end

% Set Axis properties
set(gca,'xticklabel',{'Laser off'; 'Laser on'});
set(gca,'FontSize', 16)
ylabel('replay rate session wide (events/s)')
legend({'Reverse','Forward'})
set(gcf,'Position',[20 20 300 400])

saveas(gcf,fullfile(fig_path,['replay rate session wide' num2str(rats)]),'jpg')



session_tbl = table;
session_tbl.reverse_count = [r_0_session_wide_counts; r_1_session_wide_counts];
session_tbl.forward_count = [f_0_session_wide_counts; f_1_session_wide_counts];
session_tbl.duration = [duration_0_session_wide; duration_1_session_wide];
session_tbl.session = [good_sessions; good_sessions];
rat_list = [];
for i = 1:length(good_sessions)
    rat_list(i) = unique(tbl.rat(tbl.unique_session_id == good_sessions(i)));
end
session_tbl.rat = [rat_list'; rat_list']
session_tbl.laser_state = [zeros(length(r_0_session_wide_counts),1); ones(length(r_1_session_wide_counts),1)];


glme = fitglme(session_tbl,...
    'reverse_count ~ 1 + laser_state + (1|rat:session)',...
    'Distribution','Poisson','Offset',log(session_tbl.duration))

glme = fitglme(session_tbl,...
    'forward_count ~ 1 + laser_state + (1|rat:session)',...
    'Distribution','Poisson','Offset',log(session_tbl.duration))

session_tbl2 = [session_tbl; session_tbl];
session_tbl2.replay_count(1:height(session_tbl2)/2) = session_tbl.reverse_count;
session_tbl2.replay_count(height(session_tbl2)/2+1:end) = session_tbl.forward_count;
session_tbl2.direction = [zeros(length(session_tbl.reverse_count),1);ones(length(session_tbl.forward_count),1)];
session_tbl2.rate = session_tbl2.replay_count./session_tbl2.duration;

glme = fitglme(session_tbl2, ...
    'replay_count ~ 1 + laser_state + direction + laser_state*direction + (1|rat:session)',...
    'Distribution','Poisson','Offset',log(session_tbl2.duration))



mean(session_tbl2.rate(session_tbl2.direction==0 & session_tbl2.laser_state==0))
mean(session_tbl2.rate(session_tbl2.direction==0 & session_tbl2.laser_state==1))
mean(session_tbl2.rate(session_tbl2.direction==1 & session_tbl2.laser_state==0))
mean(session_tbl2.rate(session_tbl2.direction==1 & session_tbl2.laser_state==1))

mean(session_tbl2.rate(session_tbl2.laser_state==0))
mean(session_tbl2.rate(session_tbl2.laser_state==1))
%%

dv  = [r_0_session f_0_session zeros(size([r_0_session f_0_session],1),1);r_1_session f_1_session ones(size([r_1_session f_1_session],1),1)];
dv(any(isnan(dv), 2), :) = [];
dv = array2table(dv);
% remove cells where there was no data for any of the conditions (i.e., no
% replays of that type)

% codes for column names:
dv.Properties.VariableNames = {'laserOff', 'laserON', 'direction'};
dv.direction = categorical(dv.direction);

% create the within-subjects design
withinDesign = table([0 1]','VariableNames',{'laser_state'});
withinDesign.laser_state = categorical(withinDesign.laser_state);

% create repeated measures model
rm = fitrm(dv, 'laserOff-laserON~direction', 'WithinDesign', withinDesign);

% perform anova (remove semicolon to view table generated by ranova)
AT = ranova(rm, 'WithinModel', 'laser_state');

% output an improved version of the anova table
disp(anovaTable(AT, 'Value'));

% BAR PLOT (MEAN and SEM)
% fig = two_by_two_barplot(data, group_colors, figure_xticklabels, figure_legend, figure_stats)
AT = table2array(AT);
repeated_measure_pval = AT(4,5);
group_effect_pval = AT(2,5);
interaction_pval = AT(5,5);
stats = ['rm p= ' num2str(repeated_measure_pval,2) ', group p=' num2str(group_effect_pval,2) ', interaction p=' num2str(interaction_pval,2)];
if ismember(4, rats)
    fig_ylim = [0 0.06];
else
    fig_ylim = [0 0.02];
end
set(gcf,'Position', [10 10 300 400])
fig2 = two_by_two_barplot([{r_0_session} {f_0_session} {r_1_session} {f_1_session}], [1 0 0; 0 0 1], {'laser off'; 'laser on'}, fig_ylim, fig_title, {'reverse' 'forward'}, stats);

end
%%
if plot_by_rat == 1
% Average across rats:
r_0_rat = nan(length(unique(tbl.rat)),1);
r_1_rat = nan(length(unique(tbl.rat)),1);
f_0_rat = nan(length(unique(tbl.rat)),1);
f_1_rat = nan(length(unique(tbl.rat)),1);
for i = 1:length(unique(tbl.rat))
    r_0_rat(i) = mean(reverse_replay_rate(tbl.laser_state == 0 & tbl.rat == rats(i)));
    r_1_rat(i) = mean(reverse_replay_rate(tbl.laser_state == 1 & tbl.rat == rats(i)));
    f_0_rat(i) = mean(forward_replay_rate(tbl.laser_state == 0 & tbl.rat == rats(i)));
    f_1_rat(i) = mean(forward_replay_rate(tbl.laser_state == 1 & tbl.rat == rats(i)));
end

rl_0_rat_mean = mean(r_0_rat);
rl_1_rat_mean = mean(r_1_rat);
fl_0_rat_mean = mean(f_0_rat);
fl_1_rat_mean = mean(f_1_rat);
rl_0_rat_sem = std(r_0_rat)/sqrt(length(r_0_rat));
rl_1_rat_sem = std(r_1_rat)/sqrt(length(r_1_rat));
fl_0_rat_sem = std(f_0_rat)/sqrt(length(f_0_rat));
fl_1_rat_sem = std(f_1_rat)/sqrt(length(f_1_rat));


[h,p] = ttest(r_0_rat,r_1_rat)
[h,p] = ttest(f_0_rat,f_1_rat)

y = [rl_0_rat_mean fl_0_rat_mean; rl_1_rat_mean fl_1_rat_mean];
err = [rl_0_rat_sem rl_0_rat_sem; rl_1_rat_sem fl_1_rat_sem];

figure(3); clf;
hb = bar(y); % get the bar handles
hb(1).FaceColor = [0.4940 0.1840 0.5560];
hb(2).FaceColor = [.4660 0.6740 0.1880];
hold on;
for k = 1:size(y,2)
    % get x positions per group
    xpos = hb(k).XData + hb(k).XOffset;
    % draw errorbar
    errorbar(xpos, y(:,k), err(:,k), 'LineStyle', 'none', ...
        'Color', 'k', 'LineWidth', 1);
end

% Set Axis properties
set(gca,'xticklabel',{'Laser off'; 'Laser on'});
set(gca,'FontSize', 16)
% if ismember(4,rats)
%     ylim([0 0.06])
%     yticks(linspace(0, 0.06, 5))
% else
%     ylim([0 0.02])
%     yticks(linspace(0, 0.02, 5))
% end
ylabel('Replay rate (events/s)')
legend({'Reverse','Forward'})
set(gcf,'Position', [10 10 300 400])
fig_title = ['replay rate, n = ' num2str(length(r_0_rat)) ' rats'];
title(fig_title)
saveas(gcf,fullfile(fig_path,['replay rate by rat rats' num2str(rats)]),'jpg')
end
%% Look at the rate of all replays (not just congruent replays)
reverse_replay_rate = tbl.reverse_replay./tbl.duration;
forward_replay_rate = tbl.forward_replay./tbl.duration;

r_0 = reverse_replay_rate(tbl.laser_state == 0);
r_0_mean = mean(r_0);
r_0_sem = std(r_0)/(sqrt(length(r_0)));


f_0 = forward_replay_rate(tbl.laser_state == 0);
f_0_mean = mean(f_0);
f_0_sem = std(f_0)/(sqrt(length(f_0)));


r_1 = reverse_replay_rate(tbl.laser_state == 1);
r_1_mean = mean(r_1);
r_1_sem = std(r_1)/(sqrt(length(r_1)));


f_1 = forward_replay_rate(tbl.laser_state == 1);
f_1_mean = mean(f_1);
f_1_sem = std(f_1)/(sqrt(length(f_1)));


[p,h] = ranksum(r_0,r_1)
[p,h] = ranksum(f_0,f_1)

y = [r_0_mean f_0_mean; r_1_mean f_1_mean];
err = [r_0_sem f_0_sem; f_1_sem f_1_sem];

figure(4); clf;
hb = bar(y); % get the bar handles
hb(1).FaceColor = [0.4940 0.1840 0.5560];
hb(2).FaceColor = [.4660 0.6740 0.1880];
hold on;
for k = 1:size(y,2)
    % get x positions per group
    xpos = hb(k).XData + hb(k).XOffset;
    % draw errorbar
    errorbar(xpos, y(:,k), err(:,k), 'LineStyle', 'none', ...
        'Color', 'k', 'LineWidth', 1);
end

% Set Axis properties
set(gca,'xticklabel',{'Laser off'; 'Laser on'});
set(gca,'FontSize', 16)
ylabel('Replay rate (events/s)')
legend({'Reverse','Forward'})
fig_title = ['Replay rate, n = ' num2str(length(r_0)) ' laser off, ' num2str(length(r_1)) ' laser on stopping periods'];
title(fig_title)
set(gcf,'Position', [10 10 300 400])
saveas(gcf,fullfile(fig_path,['replay rate by stopping period rats' num2str(rats)]),'jpg')

%%


combined_forward_replay_hist(bad_inds,:) = [];
combined_reverse_replay_hist(bad_inds,:) = [];

combined_forward_congruent_replay_hist(bad_inds,:) = [];
combined_reverse_congruent_replay_hist(bad_inds,:) = [];

combined_forward_incongruent_replay_hist(bad_inds,:) = [];
combined_reverse_incongruent_replay_hist(bad_inds,:) = [];


combined_congruent_replay_hist(bad_inds,:) = [];
combined_incongruent_replay_hist(bad_inds,:) = [];


combined_sde_hist(bad_inds,:) = [];
combined_replay_hist(bad_inds,:) = [];



tbl_copy(bad_inds,:) = [];

unique_sessions_in_histogram = unique(tbl_copy.unique_session_id);
% for each session, look at the # of good replays in the first 20 sec


n_laser_off_good_replays = nan(length(unique_sessions_in_histogram),1);
n_laser_on_good_replays = nan(length(unique_sessions_in_histogram),1);
for i = 1:length(unique_sessions_in_histogram)
    n_laser_off_good_reverse_replays = nansum(nansum(combined_reverse_replay_hist(tbl_copy.unique_session_id == unique_sessions_in_histogram(i) & ...
        tbl_copy.laser_state == 0,1:num_bins_to_plot)));
    n_laser_off_good_forward_replays = nansum(nansum(combined_forward_replay_hist(tbl_copy.unique_session_id == unique_sessions_in_histogram(i) & ...
        tbl_copy.laser_state == 0,1:num_bins_to_plot)));
    n_laser_on_good_reverse_replays = nansum(nansum(combined_reverse_replay_hist(tbl_copy.unique_session_id == unique_sessions_in_histogram(i) & ...
        tbl_copy.laser_state == 1,1:num_bins_to_plot)));
    n_laser_on_good_forward_replays = nansum(nansum(combined_forward_replay_hist(tbl_copy.unique_session_id == unique_sessions_in_histogram(i) & ...
        tbl_copy.laser_state == 1,1:num_bins_to_plot)));
    n_laser_off_good_replays(i) = n_laser_off_good_reverse_replays + n_laser_off_good_forward_replays;
    n_laser_on_good_replays(i) = n_laser_on_good_reverse_replays + n_laser_on_good_forward_replays;
end


good_sessions = unique_sessions_in_histogram(n_laser_off_good_replays >= min_number_replays_in_limited_stopping_period & n_laser_on_good_replays >= min_number_replays_in_limited_stopping_period);
combined_forward_replay_hist = combined_forward_replay_hist(ismember(tbl_copy.unique_session_id,good_sessions),:);
combined_reverse_replay_hist = combined_reverse_replay_hist(ismember(tbl_copy.unique_session_id,good_sessions),:);
combined_forward_congruent_replay_hist = combined_forward_congruent_replay_hist(ismember(tbl_copy.unique_session_id,good_sessions),:);
combined_reverse_congruent_replay_hist = combined_reverse_congruent_replay_hist(ismember(tbl_copy.unique_session_id,good_sessions),:);
combined_forward_incongruent_replay_hist = combined_forward_incongruent_replay_hist(ismember(tbl_copy.unique_session_id,good_sessions),:);
combined_reverse_incongruent_replay_hist = combined_reverse_incongruent_replay_hist(ismember(tbl_copy.unique_session_id,good_sessions),:);
combined_congruent_replay_hist = combined_congruent_replay_hist(ismember(tbl_copy.unique_session_id,good_sessions),:);
combined_incongruent_replay_hist = combined_incongruent_replay_hist(ismember(tbl_copy.unique_session_id,good_sessions),:);


combined_sde_hist = combined_sde_hist(ismember(tbl_copy.unique_session_id,good_sessions),:);
combined_replay_hist = combined_replay_hist(ismember(tbl_copy.unique_session_id,good_sessions),:);


tbl_copy = tbl_copy(ismember(tbl_copy.unique_session_id,good_sessions),:);

r_hist= combined_reverse_replay_hist(:,1:num_bins_to_plot);
r_hist_rate = r_hist./bin_width;
f_hist= combined_forward_replay_hist(:,1:num_bins_to_plot);
f_hist_rate = f_hist./bin_width;
r_hist_smooth = nan(size(r_hist));
f_hist_smooth = nan(size(f_hist));



r_0_hist= combined_reverse_replay_hist(tbl_copy.laser_state==0,1:num_bins_to_plot);
r_0_hist_rate = r_0_hist./bin_width;
f_0_hist= combined_forward_replay_hist(tbl_copy.laser_state==0,1:num_bins_to_plot);
f_0_hist_rate = f_0_hist./bin_width;
r_1_hist= combined_reverse_replay_hist(tbl_copy.laser_state==1,1:num_bins_to_plot);
r_1_hist_rate = r_1_hist./bin_width;
f_1_hist= combined_forward_replay_hist(tbl_copy.laser_state==1,1:num_bins_to_plot);
f_1_hist_rate = f_1_hist./bin_width;

r_0_hist_smooth = nan(size(r_0_hist));
r_1_hist_smooth = nan(size(r_1_hist));
f_0_hist_smooth = nan(size(f_0_hist));
f_1_hist_smooth = nan(size(f_1_hist));


r_cong_0_hist = combined_reverse_congruent_replay_hist(tbl_copy.laser_state==0,1:num_bins_to_plot);
r_cong_0_hist_rate = r_cong_0_hist/bin_width;
f_cong_0_hist = combined_forward_congruent_replay_hist(tbl_copy.laser_state==0,1:num_bins_to_plot);
f_cong_0_hist_rate = f_cong_0_hist/bin_width;
r_cong_1_hist = combined_reverse_congruent_replay_hist(tbl_copy.laser_state==1,1:num_bins_to_plot);
r_cong_1_hist_rate = r_cong_1_hist/bin_width;
f_cong_1_hist = combined_forward_congruent_replay_hist(tbl_copy.laser_state==1,1:num_bins_to_plot);
f_cong_1_hist_rate = f_cong_1_hist/bin_width;
r_cong_0_hist_smooth = nan(size(r_cong_0_hist));
r_cong_1_hist_smooth = nan(size(r_cong_1_hist));
f_cong_0_hist_smooth = nan(size(f_cong_0_hist));
f_cong_1_hist_smooth = nan(size(f_cong_1_hist));

r_incong_0_hist = combined_reverse_incongruent_replay_hist(tbl_copy.laser_state==0,1:num_bins_to_plot);
r_incong_0_hist_rate = r_incong_0_hist/bin_width;
f_incong_0_hist = combined_forward_incongruent_replay_hist(tbl_copy.laser_state==0,1:num_bins_to_plot);
f_incong_0_hist_rate = f_incong_0_hist/bin_width;
r_incong_1_hist = combined_reverse_incongruent_replay_hist(tbl_copy.laser_state==1,1:num_bins_to_plot);
r_incong_1_hist_rate = r_incong_1_hist/bin_width;
f_incong_1_hist = combined_forward_incongruent_replay_hist(tbl_copy.laser_state==1,1:num_bins_to_plot);
f_incong_1_hist_rate = f_incong_1_hist/bin_width;
r_incong_0_hist_smooth = nan(size(r_incong_0_hist));
r_incong_1_hist_smooth = nan(size(r_incong_1_hist));
f_incong_0_hist_smooth = nan(size(f_incong_0_hist));
f_incong_1_hist_smooth = nan(size(f_incong_1_hist));


cong_0_hist = combined_congruent_replay_hist(tbl_copy.laser_state==0,1:num_bins_to_plot);
cong_0_hist_rate = cong_0_hist/bin_width;
cong_1_hist = combined_congruent_replay_hist(tbl_copy.laser_state==1,1:num_bins_to_plot);
cong_1_hist_rate = cong_1_hist/bin_width;
incong_0_hist = combined_incongruent_replay_hist(tbl_copy.laser_state==0,1:num_bins_to_plot);
incong_0_hist_rate = incong_0_hist/bin_width;
incong_1_hist = combined_incongruent_replay_hist(tbl_copy.laser_state==1,1:num_bins_to_plot);
incong_1_hist_rate = incong_1_hist/bin_width;
cong_0_hist_smooth = nan(size(cong_0_hist));
cong_1_hist_smooth = nan(size(cong_1_hist));
incong_0_hist_smooth = nan(size(incong_0_hist));
incong_1_hist_smooth = nan(size(incong_1_hist));

sde_hist = combined_sde_hist(:,1:num_bins_to_plot);
sde_hist_rate = sde_hist/bin_width;
sde_hist_smooth = nan(size(sde_hist));

sde_0_hist = combined_sde_hist(tbl_copy.laser_state==0,1:num_bins_to_plot);
sde_0_hist_rate = sde_0_hist/bin_width;
sde_1_hist = combined_sde_hist(tbl_copy.laser_state==1,1:num_bins_to_plot);
sde_1_hist_rate = sde_1_hist/bin_width;
sde_0_hist_smooth = nan(size(sde_0_hist));
sde_1_hist_smooth = nan(size(sde_1_hist));

replay_0_hist = combined_replay_hist(tbl_copy.laser_state==0,1:num_bins_to_plot);
replay_0_hist_rate = replay_0_hist/bin_width;
replay_1_hist = combined_replay_hist(tbl_copy.laser_state==1,1:num_bins_to_plot);
replay_1_hist_rate = replay_1_hist/bin_width;
replay_0_hist_smooth = nan(size(replay_0_hist));
replay_1_hist_smooth = nan(size(replay_1_hist));

for i = 1:size(r_hist_rate,1)
    rate_sub = r_hist_rate(i,:); 
    no_nan_inds = find(~isnan(rate_sub));    
    rate_sub(isnan(rate_sub)) = [];    
    r_hist_smooth(i,no_nan_inds) = conv(rate_sub,setUp_gaussFilt([1 filter_length ],smoothing_sigma),'same');
end

for i = 1:size(f_hist_rate,1)
    rate_sub = f_hist_rate(i,:); 
        no_nan_inds = find(~isnan(rate_sub));  
        rate_sub(isnan(rate_sub)) = [];
    f_hist_smooth(i,no_nan_inds) = conv(rate_sub,setUp_gaussFilt([1 filter_length ],smoothing_sigma),'same');
end

for i = 1:size(r_0_hist_rate,1)
    rate_sub = r_0_hist_rate(i,:); 
    no_nan_inds = find(~isnan(rate_sub));    
    rate_sub(isnan(rate_sub)) = [];    
    r_0_hist_smooth(i,no_nan_inds) = conv(rate_sub,setUp_gaussFilt([1 filter_length ],smoothing_sigma),'same');
end

for i = 1:size(f_0_hist_rate,1)
    rate_sub = f_0_hist_rate(i,:); 
        no_nan_inds = find(~isnan(rate_sub));  
        rate_sub(isnan(rate_sub)) = [];
    f_0_hist_smooth(i,no_nan_inds) = conv(rate_sub,setUp_gaussFilt([1 filter_length ],smoothing_sigma),'same');
end

for i = 1:size(r_1_hist_rate,1)
    rate_sub = r_1_hist_rate(i,:); 
 no_nan_inds = find(~isnan(rate_sub));  
        rate_sub(isnan(rate_sub)) = [];
    r_1_hist_smooth(i,no_nan_inds) = conv(rate_sub,setUp_gaussFilt([1 filter_length ],smoothing_sigma),'same');
end

for i = 1:size(f_1_hist_rate,1)
    rate_sub = f_1_hist_rate(i,:); 
 no_nan_inds = find(~isnan(rate_sub));  
        rate_sub(isnan(rate_sub)) = [];
    f_1_hist_smooth(i,no_nan_inds) = conv(rate_sub,setUp_gaussFilt([1 filter_length ],smoothing_sigma),'same');
end

for i = 1:size(r_cong_0_hist_rate,1)
    rate_sub = r_cong_0_hist_rate(i,:); 
     no_nan_inds = find(~isnan(rate_sub));  
        rate_sub(isnan(rate_sub)) = [];
    r_cong_0_hist_smooth(i,no_nan_inds) = conv(rate_sub,setUp_gaussFilt([1 filter_length ],smoothing_sigma),'same');
end

for i = 1:size(f_cong_0_hist_rate,1)
    rate_sub = f_cong_0_hist_rate(i,:);
      no_nan_inds = find(~isnan(rate_sub));  
        rate_sub(isnan(rate_sub)) = [];
    f_cong_0_hist_smooth(i,no_nan_inds) = conv(rate_sub,setUp_gaussFilt([1 filter_length ],smoothing_sigma),'same');
end

for i = 1:size(r_cong_1_hist_rate,1)
    rate_sub = r_cong_1_hist_rate(i,:); 
     no_nan_inds = find(~isnan(rate_sub));  
        rate_sub(isnan(rate_sub)) = [];
    r_cong_1_hist_smooth(i,no_nan_inds) = conv(rate_sub,setUp_gaussFilt([1 filter_length ],smoothing_sigma),'same');
end

for i = 1:size(f_cong_1_hist_rate,1)
    rate_sub = f_cong_1_hist_rate(i,:); 
     no_nan_inds = find(~isnan(rate_sub));  
        rate_sub(isnan(rate_sub)) = [];
    f_cong_1_hist_smooth(i,no_nan_inds) = conv(rate_sub,setUp_gaussFilt([1 filter_length ],smoothing_sigma),'same');
end

for i = 1:size(r_incong_0_hist_rate,1)
    rate_sub = r_incong_0_hist_rate(i,:); 
         no_nan_inds = find(~isnan(rate_sub));  
        rate_sub(isnan(rate_sub)) = [];
    r_incong_0_hist_smooth(i,no_nan_inds) = conv(rate_sub,setUp_gaussFilt([1 filter_length ],smoothing_sigma),'same');
end

for i = 1:size(f_incong_0_hist_rate,1)
    rate_sub = f_incong_0_hist_rate(i,:); 
     no_nan_inds = find(~isnan(rate_sub));  
        rate_sub(isnan(rate_sub)) = [];
    f_incong_0_hist_smooth(i,no_nan_inds) = conv(rate_sub,setUp_gaussFilt([1 filter_length ],smoothing_sigma),'same');
end

for i = 1:size(r_incong_1_hist_rate,1)
    rate_sub = r_incong_1_hist_rate(i,:); 
     no_nan_inds = find(~isnan(rate_sub));  
        rate_sub(isnan(rate_sub)) = [];
    r_incong_1_hist_smooth(i,no_nan_inds) = conv(rate_sub,setUp_gaussFilt([1 filter_length ],smoothing_sigma),'same');
end

for i = 1:size(f_incong_1_hist_rate,1)
    rate_sub = f_incong_1_hist_rate(i,:); 
     no_nan_inds = find(~isnan(rate_sub));  
        rate_sub(isnan(rate_sub)) = [];
    f_incong_1_hist_smooth(i,no_nan_inds) = conv(rate_sub,setUp_gaussFilt([1 filter_length ],smoothing_sigma),'same');
end

for i = 1:size(cong_0_hist_rate,1)
    rate_sub = cong_0_hist_rate(i,:); 
     no_nan_inds = find(~isnan(rate_sub));  
        rate_sub(isnan(rate_sub)) = [];
    cong_0_hist_smooth(i,no_nan_inds) = conv(rate_sub,setUp_gaussFilt([1 filter_length ],smoothing_sigma),'same');
end

for i = 1:size(cong_1_hist_rate,1)
    rate_sub = cong_1_hist_rate(i,:); 
     no_nan_inds = find(~isnan(rate_sub));  
        rate_sub(isnan(rate_sub)) = [];
    cong_1_hist_smooth(i,no_nan_inds) = conv(rate_sub,setUp_gaussFilt([1 filter_length ],smoothing_sigma),'same');
end

for i = 1:size(incong_0_hist_rate,1)
    rate_sub = incong_0_hist_rate(i,:); 
     no_nan_inds = find(~isnan(rate_sub));  
        rate_sub(isnan(rate_sub)) = [];
    incong_0_hist_smooth(i,no_nan_inds) = conv(rate_sub,setUp_gaussFilt([1 filter_length ],smoothing_sigma),'same');
end

for i = 1:size(incong_1_hist_rate,1)
    rate_sub = incong_1_hist_rate(i,:); 
     no_nan_inds = find(~isnan(rate_sub));  
        rate_sub(isnan(rate_sub)) = [];
    incong_1_hist_smooth(i,no_nan_inds) = conv(rate_sub,setUp_gaussFilt([1 filter_length ],smoothing_sigma),'same');
end

for i = 1:size(sde_hist_rate,1)
    rate_sub = sde_hist_rate(i,:); 
     no_nan_inds = find(~isnan(rate_sub));  
        rate_sub(isnan(rate_sub)) = [];
    sde_hist_smooth(i,no_nan_inds) = conv(rate_sub,setUp_gaussFilt([1 filter_length ],smoothing_sigma),'same');
end

for i = 1:size(sde_0_hist_rate,1)
    rate_sub = sde_0_hist_rate(i,:); 
     no_nan_inds = find(~isnan(rate_sub));  
        rate_sub(isnan(rate_sub)) = [];
    sde_0_hist_smooth(i,no_nan_inds) = conv(rate_sub,setUp_gaussFilt([1 filter_length ],smoothing_sigma),'same');
end

for i = 1:size(sde_1_hist_rate,1)
    rate_sub = sde_1_hist_rate(i,:); 
     no_nan_inds = find(~isnan(rate_sub));  
        rate_sub(isnan(rate_sub)) = [];
    sde_1_hist_smooth(i,no_nan_inds) = conv(rate_sub,setUp_gaussFilt([1 filter_length ],smoothing_sigma),'same');
end


for i = 1:size(replay_0_hist_rate,1)
    rate_sub = replay_0_hist_rate(i,:);
     no_nan_inds = find(~isnan(rate_sub));  
        rate_sub(isnan(rate_sub)) = [];
    replay_0_hist_smooth(i,no_nan_inds) = conv(rate_sub,setUp_gaussFilt([1 filter_length ],smoothing_sigma),'same');
end

for i = 1:size(replay_1_hist_rate,1)
    rate_sub = replay_1_hist_rate(i,:); 
     no_nan_inds = find(~isnan(rate_sub));  
        rate_sub(isnan(rate_sub)) = [];
    replay_1_hist_smooth(i,no_nan_inds) = conv(rate_sub,setUp_gaussFilt([1 filter_length ],smoothing_sigma),'same');
end

if smooth_rate_plots == 1

    sde_hist_mean = nanmean(sde_hist_smooth);
    sde_hist_n = sum(~isnan(sde_hist_smooth),1);
    sde_hist_sem = nanstd(sde_hist_smooth)./sqrt(sde_hist_n);
    
    sde_0_hist_mean = nanmean(sde_0_hist_smooth);
    sde_0_hist_n = sum(~isnan(sde_0_hist_smooth),1);
    sde_0_hist_sem = nanstd(sde_0_hist_smooth)./sqrt(sde_0_hist_n);

    sde_1_hist_mean = nanmean(sde_1_hist_smooth);
    sde_1_hist_n = sum(~isnan(sde_1_hist_smooth),1);
    sde_1_hist_sem = nanstd(sde_1_hist_smooth)./sqrt(sde_1_hist_n);
   
    replay_0_hist_mean = nanmean(replay_0_hist_smooth);
    replay_0_hist_n = sum(~isnan(replay_0_hist_smooth),1);
    replay_0_hist_sem = nanstd(replay_0_hist_smooth)./sqrt(replay_0_hist_n);

    replay_1_hist_mean = nanmean(replay_1_hist_smooth);
    replay_1_hist_n = sum(~isnan(replay_1_hist_smooth),1);
    replay_1_hist_sem = nanstd(replay_1_hist_smooth)./sqrt(replay_1_hist_n);

    r_hist_mean = nanmean(r_hist_smooth);
    r_hist_n = sum(~isnan(r_hist_smooth),1);
    r_hist_sem = nanstd(r_hist_smooth)./sqrt(r_hist_n);
    f_hist_mean = nanmean(f_hist_smooth);
    f_hist_n = sum(~isnan(f_hist_smooth),1);
    f_hist_sem = nanstd(f_hist_smooth)./sqrt(f_hist_n);

    r_0_hist_mean = nanmean(r_0_hist_smooth);
    r_0_hist_n = sum(~isnan(r_0_hist_smooth),1);
    r_0_hist_sem = nanstd(r_0_hist_smooth)./sqrt(r_0_hist_n);
    f_0_hist_mean = nanmean(f_0_hist_smooth);
    f_0_hist_n = sum(~isnan(f_0_hist_smooth),1);
    f_0_hist_sem = nanstd(f_0_hist_smooth)./sqrt(f_0_hist_n);
    r_1_hist_mean = nanmean(r_1_hist_smooth);
    r_1_hist_n = sum(~isnan(r_1_hist_smooth),1);
    r_1_hist_sem = nanstd(r_1_hist_smooth)./sqrt(r_1_hist_n);
    f_1_hist_mean = nanmean(f_1_hist_smooth);
    f_1_hist_n = sum(~isnan(f_1_hist_smooth),1);
    f_1_hist_sem = nanstd(f_1_hist_smooth)./sqrt(f_1_hist_n);

    difference_in_replay_rate = r_hist_smooth - f_hist_smooth;    
    difference_in_replay_rate_laser_on = r_1_hist_smooth - f_1_hist_smooth;
    difference_in_replay_rate_laser_off = r_0_hist_smooth - f_0_hist_smooth;

    r_cong_0_hist_mean = nanmean(r_cong_0_hist_smooth);
    r_cong_0_hist_n = sum(~isnan(r_cong_0_hist_smooth),1);
    r_cong_0_hist_sem = nanstd(r_cong_0_hist_smooth)./sqrt(r_cong_0_hist_n);
    f_cong_0_hist_mean = nanmean(f_cong_0_hist_smooth);
    f_cong_0_hist_n = sum(~isnan(f_cong_0_hist_smooth),1);
    f_cong_0_hist_sem = nanstd(f_cong_0_hist_smooth)./sqrt(f_cong_0_hist_n);
    r_cong_1_hist_mean = nanmean(r_cong_1_hist_smooth);
    r_cong_1_hist_n = sum(~isnan(r_cong_1_hist_smooth),1);
    r_cong_1_hist_sem = nanstd(r_cong_1_hist_smooth)./sqrt(r_cong_1_hist_n);
    f_cong_1_hist_mean = nanmean(f_cong_1_hist_smooth);
    f_cong_1_hist_n = sum(~isnan(f_cong_1_hist_smooth),1);
    f_cong_1_hist_sem = nanstd(f_cong_1_hist_smooth)./sqrt(f_cong_1_hist_n);
    difference_in_cong_replay_rate_laser_on = r_cong_1_hist_smooth - f_cong_1_hist_smooth;
    difference_in_cong_replay_rate_laser_off = r_cong_0_hist_smooth - f_cong_0_hist_smooth;

    r_incong_0_hist_mean = nanmean(r_incong_0_hist_smooth);
    r_incong_0_hist_n = sum(~isnan(r_incong_0_hist_smooth),1);
    r_incong_0_hist_sem = nanstd(r_incong_0_hist_smooth)./sqrt(r_incong_0_hist_n);
    f_incong_0_hist_mean = nanmean(f_incong_0_hist_smooth);
    f_incong_0_hist_n = sum(~isnan(f_incong_0_hist_smooth),1);
    f_incong_0_hist_sem = nanstd(f_incong_0_hist_smooth)./sqrt(f_incong_0_hist_n);
    r_incong_1_hist_mean = nanmean(r_incong_1_hist_smooth);
    r_incong_1_hist_n = sum(~isnan(r_incong_1_hist_smooth),1);
    r_incong_1_hist_sem = nanstd(r_incong_1_hist_smooth)./sqrt(r_incong_1_hist_n);
    f_incong_1_hist_mean = nanmean(f_incong_1_hist_smooth);
    f_incong_1_hist_n = sum(~isnan(f_incong_1_hist_smooth),1);
    f_incong_1_hist_sem = nanstd(f_incong_1_hist_smooth)./sqrt(f_incong_1_hist_n);
    difference_in_incong_replay_rate_laser_on = r_incong_1_hist_smooth - f_incong_1_hist_smooth;
    difference_in_incong_replay_rate_laser_off = r_incong_0_hist_smooth - f_incong_0_hist_smooth;

    cong_0_hist_mean = nanmean(cong_0_hist_smooth);
    cong_0_hist_n = sum(~isnan(cong_0_hist_smooth),1);
    cong_0_hist_sem = nanstd(cong_0_hist_smooth)./sqrt(cong_0_hist_n);
    cong_1_hist_mean = nanmean(cong_1_hist_smooth);
    cong_1_hist_n = sum(~isnan(cong_1_hist_smooth),1);
    cong_1_hist_sem = nanstd(cong_1_hist_smooth)./sqrt(cong_1_hist_n); 
    incong_0_hist_mean = nanmean(incong_0_hist_smooth);
    incong_0_hist_n = sum(~isnan(incong_0_hist_smooth),1);
    incong_0_hist_sem = nanstd(incong_0_hist_smooth)./sqrt(incong_0_hist_n);
    incong_1_hist_mean = nanmean(incong_1_hist_smooth);
    incong_1_hist_n = sum(~isnan(incong_1_hist_smooth),1);
    incong_1_hist_sem = nanstd(incong_1_hist_smooth)./sqrt(incong_1_hist_n); 

else
   
    sde_hist_mean = nanmean(sde_hist);
    sde_hist_n = sum(~isnan(sde_hist),1);
    sde_hist_sem = nanstd(sde_hist)./sqrt(sde_hist_n);

    sde_0_hist_mean = nanmean(sde_0_hist);
    sde_0_hist_n = sum(~isnan(sde_0_hist),1);
    sde_0_hist_sem = nanstd(sde_0_hist)./sqrt(sde_0_hist_n);
    sde_1_hist_mean = nanmean(sde_1_hist);
    sde_1_hist_n = sum(~isnan(sde_1_hist),1);
    sde_1_hist_sem = nanstd(sde_1_hist)./sqrt(sde_1_hist_n);

    replay_0_hist_mean = nanmean(replay_0_hist);
    replay_0_hist_n = sum(~isnan(replay_0_hist),1);
    replay_0_hist_sem = nanstd(replay_0_hist)./sqrt(replay_0_hist_n);
    replay_1_hist_mean = nanmean(replay_1_hist);
    replay_1_hist_n = sum(~isnan(replay_1_hist),1);
    replay_1_hist_sem = nanstd(replay_1_hist)./sqrt(replay_1_hist_n);


 r_hist_mean = nanmean(r_hist);
    r_hist_n = sum(~isnan(r_hist),1);
    r_hist_sem = nanstd(r_hist)./sqrt(r_hist_n);
    f_hist_mean = nanmean(f_hist);
    f_hist_n = sum(~isnan(f_hist),1);
    f_hist_sem = nanstd(f_hist)./sqrt(f_hist_n);

    r_0_hist_mean = nanmean(r_0_hist);
    r_0_hist_n = sum(~isnan(r_0_hist),1);
    r_0_hist_sem = nanstd(r_0_hist)./sqrt(r_0_hist_n);
    f_0_hist_mean = nanmean(f_0_hist);
    f_0_hist_n = sum(~isnan(f_0_hist),1);
    f_0_hist_sem = nanstd(f_0_hist)./sqrt(f_0_hist_n);
    r_1_hist_mean = nanmean(r_1_hist);
    r_1_hist_n = sum(~isnan(r_1_hist),1);
    r_1_hist_sem = nanstd(r_1_hist)./sqrt(r_1_hist_n);
    f_1_hist_mean = nanmean(f_1_hist);
    f_1_hist_n = sum(~isnan(f_1_hist),1);
    f_1_hist_sem = nanstd(f_1_hist)./sqrt(f_1_hist_n);
    difference_in_replay_rate = r_hist - f_hist;
    difference_in_replay_rate_laser_on = r_1_hist - f_1_hist;
    difference_in_replay_rate_laser_off = r_0_hist - f_0_hist;

    r_cong_0_hist_mean = nanmean(r_cong_0_hist);
    r_cong_0_hist_n = sum(~isnan(r_cong_0_hist),1);
    r_cong_0_hist_sem = nanstd(r_cong_0_hist)./sqrt(r_cong_0_hist_n);
    f_cong_0_hist_mean = nanmean(f_cong_0_hist);
    f_cong_0_hist_n = sum(~isnan(f_cong_0_hist),1);
    f_cong_0_hist_sem = nanstd(f_cong_0_hist)./sqrt(f_cong_0_hist_n);
    r_cong_1_hist_mean = nanmean(r_cong_1_hist);
    r_cong_1_hist_n = sum(~isnan(r_cong_1_hist),1);
    r_cong_1_hist_sem = nanstd(r_cong_1_hist)./sqrt(r_cong_1_hist_n);
    f_cong_1_hist_mean = nanmean(f_cong_1_hist);
    f_cong_1_hist_n = sum(~isnan(f_cong_1_hist),1);
    f_cong_1_hist_sem = nanstd(f_cong_1_hist)./sqrt(f_cong_1_hist_n);
    difference_in_cong_replay_rate_laser_on = r_cong_1_hist - f_cong_1_hist;
    difference_in_cong_replay_rate_laser_off = r_cong_0_hist - f_cong_0_hist;

r_incong_0_hist_mean = nanmean(r_incong_0_hist);
r_incong_0_hist_n = sum(~isnan(r_incong_0_hist),1);
r_incong_0_hist_sem = nanstd(r_incong_0_hist)./sqrt(r_incong_0_hist_n);
f_incong_0_hist_mean = nanmean(f_incong_0_hist);
f_incong_0_hist_n = sum(~isnan(f_incong_0_hist),1);
f_incong_0_hist_sem = nanstd(f_incong_0_hist)./sqrt(f_incong_0_hist_n);
r_incong_1_hist_mean = nanmean(r_incong_1_hist);
r_incong_1_hist_n = sum(~isnan(r_incong_1_hist),1);
r_incong_1_hist_sem = nanstd(r_incong_1_hist)./sqrt(r_incong_1_hist_n);
f_incong_1_hist_mean = nanmean(f_incong_1_hist);
f_incong_1_hist_n = sum(~isnan(f_incong_1_hist),1);
f_incong_1_hist_sem = nanstd(f_incong_1_hist)./sqrt(f_incong_1_hist_n);
difference_in_incong_replay_rate_laser_on = r_incong_1_hist - f_incong_1_hist;
difference_in_incong_replay_rate_laser_off = r_incong_0_hist - f_incong_0_hist;

    cong_0_hist_mean = nanmean(cong_0_hist);
    cong_0_hist_n = sum(~isnan(cong_0_hist),1);
    cong_0_hist_sem = nanstd(cong_0_hist)./sqrt(cong_0_hist_n);
    cong_1_hist_mean = nanmean(cong_1_hist);
    cong_1_hist_n = sum(~isnan(cong_1_hist),1);
    cong_1_hist_sem = nanstd(cong_1_hist)./sqrt(cong_1_hist_n);
    incong_0_hist_mean = nanmean(incong_0_hist);
    incong_0_hist_n = sum(~isnan(incong_0_hist),1);
    incong_0_hist_sem = nanstd(incong_0_hist)./sqrt(incong_0_hist_n);
    incong_1_hist_mean = nanmean(incong_1_hist);
    incong_1_hist_n = sum(~isnan(incong_1_hist),1);
    incong_1_hist_sem = nanstd(incong_1_hist)./sqrt(incong_1_hist_n);

end

% plot smoothed difference reverse-forward replay rates
mean_diff_in_replay_rate = nanmean(difference_in_replay_rate);
sem_diff_in_replay_rate = nanstd(difference_in_replay_rate)./sqrt(r_hist_n);
pval_diff_in_replay_rate = nan(size(difference_in_replay_rate,2),1);
for i = 1:size(difference_in_replay_rate,2)
    [h, pval_diff_in_replay_rate(i)] = ttest(difference_in_replay_rate(:,i));
end


% plot smoothed difference reverse-forward replay rates
mean_diff_in_replay_rate_laser_on = nanmean(difference_in_replay_rate_laser_on);
sem_diff_in_replay_rate_laser_on = nanstd(difference_in_replay_rate_laser_on)./sqrt(r_1_hist_n);
mean_diff_in_replay_rate_laser_off = nanmean(difference_in_replay_rate_laser_off);
sem_diff_in_replay_rate_laser_off = nanstd(difference_in_replay_rate_laser_off)./sqrt(r_0_hist_n);
pval_diff_in_replay_rate_laser_off = nan(size(difference_in_replay_rate_laser_off,2),1);
for i = 1:size(difference_in_replay_rate_laser_off,2)
    [~,pval_diff_in_replay_rate_laser_off(i)] = ttest(difference_in_replay_rate_laser_off(:,i));
end

% plot smoothed difference reverse-forward congruent replay rates
mean_diff_in_cong_replay_rate_laser_on = nanmean(difference_in_cong_replay_rate_laser_on);
sem_diff_in_cong_replay_rate_laser_on = nanstd(difference_in_cong_replay_rate_laser_on)./sqrt(r_cong_1_hist_n);
mean_diff_in_cong_replay_rate_laser_off = nanmean(difference_in_cong_replay_rate_laser_off);
sem_diff_in_cong_replay_rate_laser_off = nanstd(difference_in_cong_replay_rate_laser_off)./sqrt(r_cong_0_hist_n);

% plot smoothed difference reverse-forward incongruent replay rates
mean_diff_in_incong_replay_rate_laser_on = nanmean(difference_in_incong_replay_rate_laser_on);
sem_diff_in_incong_replay_rate_laser_on = nanstd(difference_in_incong_replay_rate_laser_on)./sqrt(r_incong_1_hist_n);
mean_diff_in_incong_replay_rate_laser_off = nanmean(difference_in_incong_replay_rate_laser_off);
sem_diff_in_incong_replay_rate_laser_off = nanstd(difference_in_incong_replay_rate_laser_off)./sqrt(r_incong_0_hist_n);

%%


bin_centers = (hist_edges(1:end-1)+hist_edges(2:end))/2;

figure('Renderer', 'painters', 'Position', [20 20 600 900])
% plot smoothed congruent replay rate
subplot(2,1,1)
shadedErrorBar(bin_centers,cong_0_hist_mean,cong_0_hist_sem,'lineprops','k')
hold on
shadedErrorBar(bin_centers,cong_1_hist_mean,cong_1_hist_sem,'lineprops','r')
xlim([hist_edges(1) hist_edges(end)])

set(gca,'Xticklabel',[])
set(gca,'FontSize', 16)
ylabel('Congruent replays/s','FontSize', 16)
legend({'','laser off','','laser on'})
title('Congruent replays')
ax = gca; 
ax.FontSize = 16; 

% plot smoothed incongruent replay rate
subplot(2,1,2)
shadedErrorBar(bin_centers,incong_0_hist_mean,incong_0_hist_sem,'lineprops','k')
hold on
shadedErrorBar(bin_centers,incong_1_hist_mean,incong_1_hist_sem,'lineprops','r')
set(gca,'FontSize', 16)
ylabel('Inongruent replays/s','FontSize', 16)
legend({'','laser off','','laser on'})
title('Incongruent replays')
ylabel('Time till run onset (s)')
ax = gca; 
ax.FontSize = 16; 
%%


figure('Renderer', 'painters', 'Position', [20 20 600 900])

subplot(4,1,1)
% plot smoothed sde rate
shadedErrorBar(bin_centers,sde_hist_mean,sde_hist_sem,'lineprops','k')
set(gca,'Xticklabel',[])
set(gca,'FontSize', 16)
ylim([0 0.4])
ylabel('Events/s','FontSize', 16)
title('Spike density events')
ax = gca; 
ax.FontSize = 16; 


subplot(4,1,2)
% plot smoothed reverse replay rate
shadedErrorBar(bin_centers,r_hist_mean,r_hist_sem,'lineprops','k')
ylim([0 0.05])
set(gca,'Xticklabel',[])
set(gca,'FontSize', 16)
ylabel('Events/s','FontSize', 16)
title('Reverse replays')
ax = gca; 
ax.FontSize = 16; 


subplot(4,1,3)
% plot smoothed forward replay rate
shadedErrorBar(bin_centers,f_hist_mean,f_hist_sem,'lineprops','k')
ylim([0 0.05])
set(gca,'Xticklabel',[])
set(gca,'FontSize', 16)
ylabel('Events/s','FontSize', 16)
title('Forward replays')
ax = gca; 
ax.FontSize = 16; 


subplot(4,1,4)

% plot smoothed difference between reverse and forward replay rates

shadedErrorBar(bin_centers,mean_diff_in_replay_rate,sem_diff_in_replay_rate,'lineprops','k')
set(gca,'FontSize', 16)
yline(0)
xlabel('Time till run (s)','FontSize', 16)
ylabel('Events/s')
title('Reverse-forward replays','FontSize', 16)

saveas(gcf,fullfile(fig_path,['Replay rate aligned to run time laser state combined rats' num2str(rats)]),'jpg')

figure('Renderer', 'painters', 'Position', [20 20 600 900])

%%
subplot(4,1,1)
% plot smoothed sde rate
shadedErrorBar(bin_centers,sde_0_hist_mean,sde_0_hist_sem,'lineprops','k')
% hold on
% shadedErrorBar(bin_centers,sde_1_hist_mean,sde_1_hist_sem,'lineprops','r')
set(gca,'Xticklabel',[])
set(gca,'FontSize', 16)
ylim([0 0.4])
ylabel('Events/s','FontSize', 16)
% legend({'','laser off','','laser on'})
title('Spike density events')


subplot(4,1,2)
% plot smoothed reverse replay rate
shadedErrorBar(bin_centers,r_0_hist_mean,r_0_hist_sem,'lineprops','k')
% hold on
% shadedErrorBar(bin_centers,r_1_hist_mean,r_1_hist_sem,'lineprops','r')
ylim([0 0.05])
set(gca,'Xticklabel',[])
set(gca,'FontSize', 16)
ylabel('Events/s','FontSize', 16)
title('Reverse replays')


subplot(4,1,3)
% plot smoothed forward replay rate
shadedErrorBar(bin_centers,f_0_hist_mean,f_0_hist_sem,'lineprops','k')
% hold on
% shadedErrorBar(bin_centers,f_1_hist_mean,f_1_hist_sem,'lineprops','r')
ylim([0 0.05])
set(gca,'Xticklabel',[])
set(gca,'FontSize', 16)
ylabel('Events/s','FontSize', 16)
title('Forward replays')


subplot(4,1,4)

% plot smoothed difference between reverse and forward replay rates
shadedErrorBar(bin_centers,mean_diff_in_replay_rate_laser_off,sem_diff_in_replay_rate_laser_off,'lineprops','k')
% hold on
% shadedErrorBar(bin_centers,mean_diff_in_replay_rate_laser_on,sem_diff_in_replay_rate_laser_on,'lineprops','r')
set(gca,'FontSize', 16)
yline(0)
xlabel('Time till run (s)','FontSize', 16)
ylabel('Events/s')
title('Reverse-forward replays','FontSize', 16)

ylimit = gca().YLim; ylimit = ylimit(2);
astrisk_location = ylimit - 0.1*ylimit;
plot(bin_centers(pval_diff_in_replay_rate_laser_off<0.05),repmat(astrisk_location,size(find(pval_diff_in_replay_rate_laser_off<0.05))),'*k')


saveas(gcf,fullfile(fig_path,['Replay rate aligned to run time rats' num2str(rats)]),'jpg')
%%
figure('Renderer', 'painters', 'Position', [20 20 600 900])

subplot(4,1,1)
% plot smoothed sde rate
shadedErrorBar(bin_centers,sde_0_hist_mean,sde_0_hist_sem,'lineprops','k')
hold on
shadedErrorBar(bin_centers,sde_1_hist_mean,sde_1_hist_sem,'lineprops','r')
set(gca,'Xticklabel',[])
set(gca,'FontSize', 16)
ylabel('Events/s','FontSize', 16)
legend({'','laser off','','laser on'})
title('Spike density events')

subplot(4,1,2)
% plot smoothed reverse congruent replay rate
shadedErrorBar(bin_centers,r_cong_0_hist_mean,r_cong_0_hist_sem,'lineprops','k')
hold on
shadedErrorBar(bin_centers,r_cong_1_hist_mean,r_cong_1_hist_sem,'lineprops','r')
ylim([0 0.05])
set(gca,'Xticklabel',[])
set(gca,'FontSize', 16)
ylabel('Events/s','FontSize', 16)
legend({'','laser off','','laser on'})
title('Reverse congruent replays')

subplot(4,1,3)
% plot smoothed forward congruent replay rate
shadedErrorBar(bin_centers,f_cong_0_hist_mean,f_cong_0_hist_sem,'lineprops','k')
hold on
shadedErrorBar(bin_centers,f_cong_1_hist_mean,f_cong_1_hist_sem,'lineprops','r')
ylim([0 0.05])
set(gca,'Xticklabel',[])
set(gca,'FontSize', 16)
ylabel('Events/s','FontSize', 16)
title('Forward congruent replays')

subplot(4,1,4)
% plot smoothed difference between reverse and forward congruent replay rates
shadedErrorBar(bin_centers,mean_diff_in_cong_replay_rate_laser_off,sem_diff_in_cong_replay_rate_laser_off,'lineprops','k')
hold on
shadedErrorBar(bin_centers,mean_diff_in_cong_replay_rate_laser_on,sem_diff_in_cong_replay_rate_laser_on,'lineprops','r')
set(gca,'FontSize', 16)
title('Reverse-forward congruent replays')
xlabel('Time till run (s)','FontSize', 16)
ylabel('Events/s','FontSize', 16)
yline(0);

saveas(gcf,fullfile(fig_path,['Congruent replay rate aligned to run rats' num2str(rats)]),'jpg')

%%
figure('Renderer', 'painters', 'Position', [20 20 600 900])

subplot(4,1,1)
% plot smoothed sde rate
shadedErrorBar(bin_centers,sde_0_hist_mean,sde_0_hist_sem,'lineprops','k')
hold on
shadedErrorBar(bin_centers,sde_1_hist_mean,sde_1_hist_sem,'lineprops','r')
set(gca,'Xticklabel',[])
set(gca,'FontSize', 16)
ylabel('Events/s','FontSize', 16)
legend({'','laser off','','laser on'})
title('Spike density events')

subplot(4,1,2)
% plot smoothed reverse incongruent replay rate
shadedErrorBar(bin_centers,r_incong_0_hist_mean,r_incong_0_hist_sem,'lineprops','k')
hold on
shadedErrorBar(bin_centers,r_incong_1_hist_mean,r_incong_1_hist_sem,'lineprops','r')
ylim([0 0.01])
set(gca,'Xticklabel',[])
set(gca,'fontsize',16)
ylabel('Events/s','FontSize', 16)
legend({'','laser off','','laser on'})
title('Reverse incongruent replays')

subplot(4,1,3)
% plot smoothed forward incongruent replay rate
shadedErrorBar(bin_centers,f_incong_0_hist_mean,f_incong_0_hist_sem,'lineprops','k')
hold on
shadedErrorBar(bin_centers,f_incong_1_hist_mean,f_incong_1_hist_sem,'lineprops','r')
ylim([0 0.01])
set(gca,'Xticklabel',[])
set(gca,'fontsize',16)
ylabel('Events/s','FontSize',16)
title('Forward incongruent replays')

subplot(4,1,4)
% plot smoothed difference between reverse and forward incongruent replay rates
shadedErrorBar(bin_centers,mean_diff_in_incong_replay_rate_laser_off,sem_diff_in_incong_replay_rate_laser_off,'lineprops','k')
hold on
shadedErrorBar(bin_centers,mean_diff_in_incong_replay_rate_laser_on,sem_diff_in_incong_replay_rate_laser_on,'lineprops','r')
set(gca,'FontSize',16)
title('Reverse-forward incongruent replays')
xlabel('Time till run (s)','FontSize',16)
ylabel('Events/s','FontSize',16)
yline(0);

saveas(gcf,fullfile(fig_path,['incongruent replay rate aligned to run rats' num2str(rats)]),'jpg')

%%

num_laps_to_plot = 30;

rl_0_time_in_stopping_period_by_lap = nan(num_laps_to_plot,num_bins_to_plot);
rl_1_time_in_stopping_period_by_lap = nan(num_laps_to_plot,num_bins_to_plot);
r_0_time_in_stopping_period_by_lap = nan(num_laps_to_plot,num_bins_to_plot);
r_1_time_in_stopping_period_by_lap = nan(num_laps_to_plot,num_bins_to_plot);
fl_0_time_in_stopping_period_by_lap = nan(num_laps_to_plot,num_bins_to_plot);
fl_1_time_in_stopping_period_by_lap = nan(num_laps_to_plot,num_bins_to_plot);
f_0_time_in_stopping_period_by_lap = nan(num_laps_to_plot,num_bins_to_plot);
f_1_time_in_stopping_period_by_lap = nan(num_laps_to_plot,num_bins_to_plot);
for lap = 1:num_laps_to_plot
    inds_0 = find(tbl_copy.laser_state == 0 & tbl_copy.laser_state_stopping_period == lap);
    inds_1 = find(tbl_copy.laser_state == 1 & tbl_copy.laser_state_stopping_period == lap);  
    r_0_time_in_stopping_period_by_lap(lap,:) = nansum(combined_reverse_replay_hist(inds_0,1:num_bins_to_plot),1);  
    r_1_time_in_stopping_period_by_lap(lap,:) = nansum(combined_reverse_replay_hist(inds_1,1:num_bins_to_plot),1);
    f_0_time_in_stopping_period_by_lap(lap,:) = nansum(combined_forward_replay_hist(inds_0,1:num_bins_to_plot),1);
    f_1_time_in_stopping_period_by_lap(lap,:) = nansum(combined_forward_replay_hist(inds_1,1:num_bins_to_plot),1);
end

green_colormap = customcolormap(linspace(0,1,5), {'#014419','#1c7735','#5aae60','#a6db9d','#d7f1d6'});
purple_colormap = customcolormap(linspace(0,1,5), {'#410149','#762a84','#9b6fac','#c1a5cd','#e7d4e8'});



figure()
ax(1) = subplot(2,2,1)
imagesc(r_0_time_in_stopping_period_by_lap);
colormap(ax(1),purple_colormap)
caxis([0 2])
% xticks()
% xticks([1 plot_length/2 plot_length])
% xticklabels([0 time_to_plot/2 time_to_plot])
set(gca,'Xticklabel',[])
yticks([0:10:num_laps_to_plot])
ylabel('Lap number')
title('Laser OFF')
ax = gca; 
ax.FontSize = 10;

ax(2) = subplot(2,2,2)
imagesc(r_1_time_in_stopping_period_by_lap);
colormap(ax(2),purple_colormap)
caxis([0 2])
% xticks([1 plot_length/2 plot_length])
% xticklabels([0 time_to_plot/2 time_to_plot])
set(gca,'Xticklabel',[])
yticks([0:10:num_laps_to_plot])
set(gca,'YTickLabel',[])
title('Laser ON')
ax = gca; 
ax.FontSize = 10;


ax(3) = subplot(2,2,3)
imagesc(f_0_time_in_stopping_period_by_lap);
colormap(ax(3),green_colormap)
caxis([0 2])
% xticks([1 plot_length/2 plot_length])
% xticklabels([0 time_to_plot/2 time_to_plot])
xticks (linspace(1,length(bin_centers),6))
xticklabels({'-10'; '-8'; '-6'; '-4'; '-2'; '0'})
xlabel('Time till run (s)')
yticks([0:10:num_laps_to_plot])
ylabel('Lap number')
ax = gca; 
ax.FontSize = 10;




ax(4) = subplot(2,2,4)
imagesc(f_1_time_in_stopping_period_by_lap);
colormap(ax(4),green_colormap)
caxis([0 2])
xticks (linspace(1,length(bin_centers),6))
xticklabels({'-10'; '-8'; '-6'; '-4'; '-2'; '0'})
set(gca,'YTickLabel',[])
xlabel('Time till run (s)')
ax = gca; 
ax.FontSize = 10;


saveas(gcf,fullfile(fig_path,['Replays over stopping period rats' num2str(rats)]),'jpg')



%%

[rho_1,pval_1] = corr(tbl.lap(tbl.laser_state==0),reverse_replay_rate(tbl.laser_state==0));
[rho_2,pval_2] = corr(tbl.lap(tbl.laser_state==1),reverse_replay_rate(tbl.laser_state==1));
coeffs_1 = polyfit(tbl.lap(tbl.laser_state==0),reverse_replay_rate(tbl.laser_state==0),1);
coeffs_2 = polyfit(tbl.lap(tbl.laser_state==1),reverse_replay_rate(tbl.laser_state==1),1);

figure()
plot(tbl.lap(tbl.laser_state==0),reverse_replay_rate(tbl.laser_state==0),'ob')
hold on
plot(tbl.lap(tbl.laser_state==1),reverse_replay_rate(tbl.laser_state==1),'or')
x1 = 1:100;
y1 = polyval(coeffs_1, x1);
x2 = 1:100;
y2 = polyval(coeffs_2, x2);
hold on
plot(x1,y1,'b','LineWidth',2)
hold on
plot(x2,y2,'r','LineWidth',2)

xlabel('Pass number','Interpreter','none')
ylabel('reverse replay rate','Interpreter','none')

dim = [.2 .5 .3 .3];
str = {['Laser off p='  num2str(pval_1,2)]; ['Laser on p=' num2str(pval_2,2)]};
annotation('textbox',dim,'String',str,'FitBoxToText','on');



%%
[rho_1,pval_1] = corr(tbl.lap(tbl.laser_state==0),forward_replay_rate(tbl.laser_state==0));
[rho_2,pval_2] = corr(tbl.lap(tbl.laser_state==1),forward_replay_rate(tbl.laser_state==1));
coeffs_1 = polyfit(tbl.lap(tbl.laser_state==0),forward_replay_rate(tbl.laser_state==0),1);
coeffs_2 = polyfit(tbl.lap(tbl.laser_state==1),forward_replay_rate(tbl.laser_state==1),1);

figure()
plot(tbl.lap(tbl.laser_state==0),forward_replay_rate(tbl.laser_state==0),'ob')
hold on
plot(tbl.lap(tbl.laser_state==1),forward_replay_rate(tbl.laser_state==1),'or')
x1 = 1:100;
y1 = polyval(coeffs_1, x1);
x2 = 1:100;
y2 = polyval(coeffs_2, x2);
hold on
plot(x1,y1,'b','LineWidth',2)
hold on
plot(x2,y2,'r','LineWidth',2)

xlabel('Pass number','Interpreter','none')
ylabel('forward replay rate','Interpreter','none')

dim = [.2 .5 .3 .3];
str = {['Laser off p='  num2str(pval_1,2)]; ['Laser on p=' num2str(pval_2,2)]};
annotation('textbox',dim,'String',str,'FitBoxToText','on');



