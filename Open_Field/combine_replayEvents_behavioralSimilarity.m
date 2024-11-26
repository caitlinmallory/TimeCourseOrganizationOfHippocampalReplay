Rat_Names = {'Clover','Bo','MEC1','Bolt','Dash','CM1','Janni','Harpy','Imp','W18','W19','Billy3','Curly2','Goethe2'};
experimental_rats = [1 2 3 6]; % opto; with Jaws
control_rats1 = [4 5 ]; % opto, gfp only
control_rats2 = [7 8 9 10 11 12 13 14]; % no laser, no injections
binSize = 2;
data_tbl = table();

rats=[12 13 14];
load_shuffles=0;

session_id = 0;
for rat_num = 1:length(rats)
    rat = rats(rat_num);
    load_open_field_session_list

    for day = 1:length(dayFiles)
        display([Rat_Names{rat},' Day ', dayFiles{day}])
        cd(fullfile(directory,dayFiles{day}))

        load Experiment_Information

        sessions_that_meet_criterion_day = check_for_relevant_sessions(flags,1,Experiment_Information,must_include_all_flags,1,novel_range);

        if ~isempty(sessions_that_meet_criterion_day)

            for session_count = 1:length(sessions_that_meet_criterion_day)
                %for session_count = 1
                session_id = session_id + 1;
                sub_session_num = sessions_that_meet_criterion_day(session_count);

                sessionNum = sub_session_num;
                sessionNum_decoder = Experiment_Information.Segments(sessionNum).Decoder;
                sessionDecoding_error = Experiment_Information.Segments(sessionNum).decodingError;

                if sessionDecoding_error < session_decoding_accuracy_thr

                    if look_at_speed_thresholded_events == 1
                        load replayEvents
                        find_exact_stopping_period = 1;
                    else
                        load replayEvents_no_speed_thr.mat
                        find_exact_stopping_period = 0;
                    end

                    % Sometimes values weren't assigned to entries of
                    % replayEvents- we need to replace these with nans.
                    properties_in_decoder_replay = fieldnames(decoder_replay(sessionNum_decoder).replayEvents);
                    for p = 1:length(properties_in_decoder_replay)
                        maskEmptyId = arrayfun(@(a) isempty(a.(properties_in_decoder_replay{p})), decoder_replay(sessionNum_decoder).replayEvents);
                        [decoder_replay(sessionNum_decoder).replayEvents(maskEmptyId).(properties_in_decoder_replay{p})] = deal(nan);
                    end

                    t_sub = struct2table(decoder_replay(sessionNum_decoder).replayEvents);

                    varNames_data_tbl = data_tbl.Properties.VariableNames;
                    varNames_t_sub = t_sub.Properties.VariableNames;

                    % In general we will only include variables that have
                    % been generated for all sessions. The exceptions for
                    % now are time_since_real_drink_onset and
                    % time_between_anticipatory_locking_and_drink_start
                    % (the latter is only relevant to John's data).
                    if ~ismember('time_since_real_drink_onset',varNames_t_sub)
                        t_sub.time_since_real_drink_onset = nan([height(t_sub),1]);
                        varNames_t_sub = [varNames_t_sub {'time_since_real_drink_onset'}];
                    end


                    if ~ismember('meanAngDisplacement_HD',varNames_t_sub)
                        t_sub.meanAngDisplacement_HD = nan([height(t_sub),1]);
                        varNames_t_sub = [varNames_t_sub {'meanAngDisplacement_HD'}];
                    end

                    if ~ismember('meanAngDisplacement_HD_and_past_path',varNames_t_sub)
                        t_sub.meanAngDisplacement_HD_and_past_path = nan([height(t_sub),1]);
                        varNames_t_sub = [varNames_t_sub {'meanAngDisplacement_HD_past_path'}];
                    end

                    if ~ismember('meanAngDisplacement_HD_and_future_path',varNames_t_sub)
                        t_sub.meanAngDisplacement_HD_and_future_path = nan([height(t_sub),1]);
                        varNames_t_sub = [varNames_t_sub {'meanAngDisplacement_HD_future_path'}];
                   end

                    if ~ismember('time_between_anticipatory_licking_and_drink_start',varNames_t_sub)
                        t_sub.time_between_anticipatory_licking_and_drink_start = zeros([height(t_sub),1]);
                        varNames_t_sub = [varNames_t_sub {'time_between_anticipatory_licking_and_drink_start'}];
                    end

                    if ~ismember('percent_participation_excitatory',varNames_t_sub)
                        t_sub.percent_participation_excitatory = zeros([height(t_sub),1]);
                        varNames_t_sub = [varNames_t_sub {'percent_participation_excitatory'}];
                    end
                    if ~ismember('percent_participation_inhibitory',varNames_t_sub)
                        t_sub.percent_participation_inhibitory = zeros([height(t_sub),1]);
                        varNames_t_sub = [varNames_t_sub {'percent_participation_inhibitory'}];
                    end


                    if isempty(data_tbl)
                        vars_to_keep = varNames_t_sub;
                    else
                        vars_to_keep = intersect(varNames_t_sub,varNames_data_tbl);
                    end

                    vars_to_keep(strcmp(vars_to_keep,'time_since_stopping_period_onset')==1) = [];

                    t_sub = t_sub(:,vars_to_keep);


                    if load_shuffles == 0
                        % Remove shuffles from table- because these haven't
                        % been calculated in the same way for all sessions.
                        if sum(strcmp(t_sub.Properties.VariableNames,'all_angularDisplacement_past_shuffled'))==1
                            t_sub.all_angularDisplacement_past_shuffled = [];
                        end
                        if sum(strcmp(t_sub.Properties.VariableNames,'all_angularDisplacement_future_shuffled'))==1
                            t_sub.all_angularDisplacement_future_shuffled = [];
                        end
                        if sum(strcmp(t_sub.Properties.VariableNames,'replay_shuffle_past'))==1
                            t_sub.replay_shuffle_past = [];
                        end
                        if sum(strcmp(t_sub.Properties.VariableNames,'replay_shuffle_future'))==1
                            t_sub.replay_shuffle_future = [];
                        end
                    end

                    timePoints = mean(t_sub.timePoints,2);
                    ind_session = ismember(timePoints,compute_dataTemporalConcatenation(timePoints,Experiment_Information.Segments(sessionNum).Times));
                    t_sub = t_sub(ind_session,:);

                    % data cleanup: weird cases where there was no future
                    % or past path due to the replay being right at the
                    % start or end of a session.
                    bad_inds = [];
                    for i = 1:height(t_sub)
                        if prod(size(t_sub.angDisplacement_futPath{i}) == size(t_sub.angDisplacement_pastPath{i})) ~= 1
                            bad_inds = [bad_inds; i];

                        end
                    end
                    t_sub(bad_inds,:) = [];

                    session_flags = Experiment_Information.Segments(sub_session_num).Flags;
                    % check the range of session numbers that are to be considered in terms
                    % of novelty
                    [~,index_of_novelty_flag] = find(floor(session_flags) == 10);
                    novel_session_number = session_flags(index_of_novelty_flag) - 10;

                    t_sub.rat_label(:) = repmat(rat,height(t_sub),1);
                    t_sub.day_label(:) = repmat(day,height(t_sub),1);
                    t_sub.segment_label(:) = repmat(sub_session_num,height(t_sub),1);
                    t_sub.unique_session_id(:) = repmat(session_id, height(t_sub),1);
                    t_sub.session_str(:) = repmat({fullfile(directory,dayFiles{day})},height(t_sub),1);
                    t_sub.session_flags(:) = repmat({session_flags},height(t_sub),1);

                    if iscell(t_sub.all_angles_between_past_future_trajectory)
                        for i = 1:height(t_sub)
                            if isnan(t_sub.all_angles_between_past_future_trajectory{i})
                                t_sub.all_angles_between_past_future_trajectory{i} = nan(1,90) ;
                            end
                        end
                        t_sub.all_angles_between_past_future_trajectory = cell2mat(t_sub.all_angles_between_past_future_trajectory);
                    end
                    if iscell(t_sub.all_angles_between_past_future_trajectory_inner_90)
                        for i = 1:height(t_sub)
                            if isnan(t_sub.all_angles_between_past_future_trajectory_inner_90{i})
                                t_sub.all_angles_between_past_future_trajectory_inner_90{i} = nan(1,2);
                            end
                        end
                        t_sub.all_angles_between_past_future_trajectory_inner_90 = cell2mat(t_sub.all_angles_between_past_future_trajectory_inner_90);
                    end
                    if iscell(t_sub.all_angles_between_past_future_trajectory_inner_80)
                        for i = 1:height(t_sub)
                            if isnan(t_sub.all_angles_between_past_future_trajectory_inner_80{i})
                                t_sub.all_angles_between_past_future_trajectory_inner_80{i} = nan(1,2);
                            end
                        end
                        t_sub.all_angles_between_past_future_trajectory_inner_80 = cell2mat(t_sub.all_angles_between_past_future_trajectory_inner_80);
                    end
                    data_tbl = [data_tbl; t_sub];

                    % Additionally, pull out all crossings

                end
            end
        end
    end
end


%%

for i = 1:height(data_tbl)
    data_tbl.laser_state_binary(i) = mode(data_tbl.laser_state{i});
end
if ismember(14,flags_to_include)
    data_tbl(data_tbl.laser_state_binary==0,:) = [];
end
