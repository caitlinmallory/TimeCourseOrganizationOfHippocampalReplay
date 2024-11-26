rats = [1:11];

novel_range = [10.1 inf]; %inclusive range
flags_to_include = [11];
flags_to_exclude = [0 1002 1003];
flags.flags_to_include = flags_to_include;
flags.flags_to_exclude = flags_to_exclude;
must_include_all_flags = 1; % 1: requires that all flags in Flags_to_include are present. 0: requires that at least one flag in Flags_to_include are present.
data_tbl = table();

Rat_Names = {'Clover','Bo','MEC1','Bolt','Dash','CM1','Janni','Harpy','Imp','W18','W19'};

for rat_num = 1:length(rats)
    rat = rats(rat_num);
    load_linear_track_session_list

    candidate_event_count = 1;

    for day = 1:length(dayFiles)
        display([Rat_Names{rat},' Day ', dayFiles{day}])
        cd(fullfile(directory,dayFiles{day}))

        load Experiment_Information
        load Behavior_Data

        sessions_that_meet_criterion_day = check_for_relevant_sessions(flags,1,Experiment_Information,must_include_all_flags,1,novel_range);

        if ~isempty(sessions_that_meet_criterion_day)

            for session_count = 1:length(sessions_that_meet_criterion_day)
             
                sub_session_num = sessions_that_meet_criterion_day(session_count);
                load cluster_participation_in_replays.mat
                data_tbl = [data_tbl; cluster_participation_table];

            end
        end
    end
end

%%
bin_width = 0.5;
time_vec = bin_width/2:bin_width:10-bin_width/2;
define_by_place_fields=1;
define_by_place_fields_and_modulation=0;
define_by_modulation=0;
a = data_tbl.clusters_left_field_only & data_tbl.Well_Isolated==1;
b = data_tbl.clusters_right_field_only & data_tbl.Well_Isolated==1;
c = data_tbl.clusters_right_and_left_field & data_tbl.Well_Isolated==1;

% define_by_place_fields=0;
% define_by_place_fields_and_modulation=1;
% define_by_modulation=0;
% a = data_tbl.clusters_left_field_only & data_tbl.lr_firing_rate_modulation>0.5 & data_tbl.Well_Isolated==1;
% b = data_tbl.clusters_right_field_only & data_tbl.lr_firing_rate_modulation>0.5 & data_tbl.Well_Isolated==1;
% c = data_tbl.clusters_right_and_left_field & data_tbl.Well_Isolated==1;

% define_by_place_fields=0;
% define_by_place_fields_and_modulation=0;
% define_by_modulation=1;
% a = data_tbl.lr_firing_rate_modulation>0.5 & data_tbl.clusters_left_max_rate_higher & data_tbl.Well_Isolated==1;
% b = data_tbl.lr_firing_rate_modulation>0.5 & data_tbl.clusters_right_max_rate_higher & data_tbl.Well_Isolated==1;
% c = data_tbl.lr_firing_rate_modulation<=0.5 & data_tbl.Well_Isolated==1;

figure()
cells_active_past = [data_tbl.percent_participation_left(a==1,:); data_tbl.percent_participation_right(b==1,:)];
cells_active_future = [data_tbl.percent_participation_left(b==1,:); data_tbl.percent_participation_right(a==1,:)];
cells_active_both = [data_tbl.percent_participation_left(c==1,:); data_tbl.percent_participation_right(c==1,:)];


past_cells_future_window = nanmean(cells_active_past(:,time_vec<=3),2);
past_cells_past_window = nanmean(cells_active_past(:,time_vec>3),2);
[p,h,z]= signrank(past_cells_future_window,past_cells_past_window)

future_cells_future_window = nanmean(cells_active_future(:,time_vec<=3),2);
future_cells_past_window = nanmean(cells_active_future(:,time_vec>3),2);
[p,h,z]= signrank(future_cells_future_window,future_cells_past_window)

bi_cells_future_window = nanmean(cells_active_both(:,time_vec<=3),2);
bi_cells_past_window = nanmean(cells_active_both(:,time_vec>3),2);
[p,h,z]= signrank(bi_cells_future_window,bi_cells_past_window)


past_mean = nanmean(cells_active_past);
past_sem = nanstd(cells_active_past)./sqrt(length(cells_active_past));

future_mean = nanmean(cells_active_future);
future_sem = nanstd(cells_active_future)./sqrt(length(cells_active_future));

both_mean = nanmean(cells_active_both);
both_sem = nanstd(cells_active_both)./sqrt(length(cells_active_both));

sum_mean = nanmean(cells_active_past) + nanmean(cells_active_future);

figure('Position',[150 150 500 500])
time_vec = bin_width/2:bin_width:10-bin_width/2;
shadedErrorBar(time_vec,past_mean,past_sem,'lineProps','-b');
hold on
shadedErrorBar(time_vec,future_mean,future_sem,'lineProps','-r');
hold on
shadedErrorBar(time_vec,both_mean,both_sem,'lineProps','-m');
hold on
plot(time_vec,sum_mean,'color',[0 0 0])
xlabel('Time since stopping')
ylabel('Percentage of replays participating in');
legend({'active past lap','active future lap','active both'})
ylim([-0 0.6])

set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
if define_by_place_fields==1
    saveas(gcf,'/home/caitlin/Data/Processed_Data/Manuscript_Figures/CellParticipationInReplays_UsingFields','pdf');
elseif define_by_modulation==1
saveas(gcf,'/home/caitlin/Data/Processed_Data/Manuscript_Figures/CellParticipationInReplays_UsingFiringRateModulation','pdf');
elseif define_by_place_fields_and_modulation==1
saveas(gcf,'/home/caitlin/Data/Processed_Data/Manuscript_Figures/CellParticipationInReplays_UsingFieldsAndFiringRateModulation','pdf');
end    