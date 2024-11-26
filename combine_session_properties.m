
rats = [1:11];
rats = [1 4 6 12 13 14]
% rats = 1:14

include_linear = 0;
include_open_field = 1;
combine_control_epoch_data = 0;
must_include_all_flags = 1; % 1: requires that all flags in Flags_to_include are present. 0: requires that at least one flag in Flags_to_include are present.
novel_range = [10.2 10.9];
flags.flags_to_include = [11];
% 18 = dreadds; 19 = saline; 13 = laser off; 14 = laser on
flags.flags_to_exclude = [1003 1002 0  101 14 18];
remove_epochs_from_unchanged_reward_end = 0;
remove_epochs_from_changed_reward_end = 0;
hist_bin_width = 0.5; %seconds
reward_size = [0 inf];
% Session quality criterion:
min_num_replays_per_session = 0; % this will remove sessions with fewer than min_num_replays_per_session in each epoch (Laser on and Laser off) from analysis
pnct_directional_decoding_correct_thr = 0.70;
median_decoding_error_thr = 5;

smooth_rate_plots = 1;
smoothing_sigma = 1; % bins

Rat_Names = {'Clover','Bo','MEC1','Bolt','Dash','CM1','Janni','Harpy','Imp','W18','W19','Billy3','Curly2','Goethe2'};
experimental_rats = [1 2 3 6]; % opto; with Jaws
control_rats1 = [4 5 ]; % opto, gfp only
control_rats2 = [7 8 9 10 11]; % no laser, no injections
session_id = 0;

session_num_trials = [];
error_table = table();
session_rat = {};
session_rat_num = [];
session_unique_session_id = [];
session_label = {};
session_flag_list = [];
session_decoding_error = [];
session_pcnt_correct_directional_decoding = [];
session_num_place_cells = [];

for rat = rats
    clear dayFiles1
    clear dayFiles2
    clear dayFiles
    if include_linear==1
        load_linear_track_session_list
        try
            dayFiles1 = fullfile(directory,dayFiles);
        catch
            continue
        end
    else
        dayFiles1 = {};
    end
    if include_open_field==1
        load_open_field_session_list
        try
            dayFiles2 = fullfile(directory,dayFiles);
        catch
            continue
        end
    else
        dayFiles2 = {};
    end
    dayFiles = [dayFiles1 dayFiles2]


    for day = 1:length(dayFiles)
        display([Rat_Names{rat},' Day ', dayFiles{day}])
        cd(dayFiles{day})

        load Experiment_Information.mat
        sessions_that_meet_criterion_day = check_for_relevant_sessions(flags,1,Experiment_Information,must_include_all_flags,1,novel_range);
     
        % Pull out all unique run sessions in this day:
        run_segments = [];
        for i = 1:length(Experiment_Information.Segments)
            if ismember(11,Experiment_Information.Segments(i).Flags)
                run_segments = [run_segments; i];
            end
        end

        if ~isempty(sessions_that_meet_criterion_day)

            for session_count = 1:length(sessions_that_meet_criterion_day)
                session_id = session_id + 1;
                sub_session_num = sessions_that_meet_criterion_day(session_count);
                clear median_decoding_error
                clear percent_correct_directional_assigment
                sub_session_ind = find(run_segments == sub_session_num)

                if Experiment_Information.spatialDim==1
                    load session_wide_properties.mat
                    session_decoding_error = [session_decoding_error; median_decoding_error];
                    session_pcnt_correct_directional_decoding = [session_pcnt_correct_directional_decoding; percent_correct_directional_assigment];
                    session_num_place_cells = [session_num_place_cells; number_place_cells];

                    load Behavior_Data
                    session_num_trials = [session_num_trials; length(Reward_Epoch_Time_Boundaries)];
                else
                    session_decoding_error = [session_decoding_error; Experiment_Information.Segments(sub_session_num).decodingError];
                    load clusters
                    session_num_place_cells = [session_num_place_cells; length(clusters)];
                    session_pcnt_correct_directional_decoding = [session_pcnt_correct_directional_decoding; nan];

                    load true_drink_periods
                    session_num_trials = [session_num_trials; length(true_drink_periods_summary(sub_session_ind).true_drink_periods)]
                end


                sessionNum_decoder = Experiment_Information.Segments(sub_session_num).Decoder;
                session_flags = Experiment_Information.Segments(sub_session_num).Flags;
                % check the range of session numbers that are to be considered in terms
                % of novelty
                [~,index_of_novelty_flag] = find(floor(session_flags) == 10);
                novel_session_number = session_flags(index_of_novelty_flag) - 10;

                % Save a number of session-wide properties
                session_rat = [session_rat; [Rat_Names{rat}]];
                session_rat_num = [session_rat_num; rat];
                session_unique_session_id = [session_unique_session_id; session_id];
                session_label = [session_label; [Rat_Names{rat},' Day ', dayFiles{day}]];
                session_flag_list = [session_flag_list; {session_flags}];


            end
        end
    end

end

data = session_num_trials
good_sessions = ones(length(session_decoding_error),1);
good_sessions(session_decoding_error > 5) = 0;
good_sessions(session_pcnt_correct_directional_decoding < 0.7) = 0;
good_sessions=logical(good_sessions)

data = session_decoding_error(good_sessions);
% data = session_pcnt_correct_directional_decoding;
data = session_num_place_cells;
%data = session_num_trials;

min(data)
max(data)

nanmean(data)
nanstd(data)/(sqrt(sum(~isnan(data))))