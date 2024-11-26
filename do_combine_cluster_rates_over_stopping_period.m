
rats = [1:11];
novel_range = [10.1 inf]; %inclusive range
flags_to_include = [11];
flags_to_exclude = [0 1002 1003];
flags.flags_to_include = flags_to_include;
flags.flags_to_exclude = flags_to_exclude;
must_include_all_flags = 1; % 1: requires that all flags in Flags_to_include are present. 0: requires that at least one flag in Flags_to_include are present.
long_replays_only = 1;
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
                load cluster_rates_during_stops_v1.mat

                data_tbl = [data_tbl; cluster_rate_table];
            end
        end
    end
end
%%
plot_cell_firing_rates_over_time_since_stopping