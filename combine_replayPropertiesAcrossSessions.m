





%%
% Save a number of session-wide properties
session_rat = {};
session_rat_num = [];
session_unique_session_id = [];
session_label = {};
session_flag_list = {};
session_directional_field_corr = [];
session_firing_rate_modulation = [];
session_decoding_error = [];
session_pcnt_correct_directional_decoding = [];
session_num_place_cells  = [];
session_ripple_rate  = [];
session_ripple_rate_0  = [];
session_ripple_rate_1 = [];
session_sde_rate = [];
session_sde_rate_0 = [];
session_sde_rate_1 = [];
session_novel_num  = [];
session_num_passes = [];
session_num_passes_2 = [];
session_time_stationary_at_endzone = [];
session_time_stationary_at_endzone_0 = [];
session_time_stationary_at_endzone_1 = [];
session_mean_stopping_period_duration = [];
session_median_stopping_period_duration = [];
session_mean_drink_duration = [];
session_median_drink_duration = [];
session_mean_stopping_period_hesitation = [];
session_median_stopping_period_hesitation = [];
session_reward_size = [];


data_tbl = table();

session_id = 0;
if reload_data == 1
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
                session_id = session_id + 1;
                sub_session_num = sessions_that_meet_criterion_day(session_count);
                load session_wide_properties.mat

                % check that the session has enough laps
                num_passes_in_session = compute_dataTemporalConcatenation(Pass_Info,Experiment_Information.Segments(sub_session_num).Times);
                num_passes_in_session = length(unique(num_passes_in_session(:,2)));
                disp(['num passes =' num2str(num_passes_in_session)])

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
                session_directional_field_corr = [session_directional_field_corr; nanmedian(directional_correlation)];
                session_firing_rate_modulation = [session_firing_rate_modulation; nanmedian(firing_rate_modulation)];
                session_decoding_error = [session_decoding_error; median_decoding_error];
                session_pcnt_correct_directional_decoding = [session_pcnt_correct_directional_decoding; percent_correct_directional_assigment];
                session_num_place_cells = [session_num_place_cells; number_place_cells];


                if distance_from_end_tight_threshold == 5
                    column = 1;
                elseif distance_from_end_tight_threshold == 10
                    column = 2;
                elseif distance_from_end_tight_threshold == 30
                    column = 3;
                else
                    keyboard
                end

                session_ripple_rate = [session_ripple_rate; ripple_rate(column)];
                session_ripple_rate_0 = [session_ripple_rate_0; ripple_rate_0(column)];
                session_ripple_rate_1 = [session_ripple_rate_1; ripple_rate_1(column)];
                session_sde_rate = [session_sde_rate; sde_rate(column)];
                session_sde_rate_0 = [session_sde_rate_0; sde_rate_0(column)];
                session_sde_rate_1 = [session_sde_rate_1; sde_rate_1(column)];
                session_time_stationary_at_endzone = [session_time_stationary_at_endzone; amount_of_time_stationary_at_endzone(column)];
                session_time_stationary_at_endzone_0 = [session_time_stationary_at_endzone_0; amount_of_time_stationary_at_endzone_0(column)];
                session_time_stationary_at_endzone_1 = [session_time_stationary_at_endzone_1; amount_of_time_stationary_at_endzone_1(column)];


                session_novel_num = [session_novel_num; novel_session_number];
                session_num_passes = [session_num_passes; num_passes_in_session];
                session_num_passes_2 = [session_num_passes_2; num_passes];
                session_mean_stopping_period_duration = [session_mean_stopping_period_duration; mean_stopping_period_duration];
                session_median_stopping_period_duration = [session_median_stopping_period_duration; median_stopping_period_duration];
                session_mean_drink_duration = [session_mean_drink_duration; mean_drink_duration];
                session_median_drink_duration = [session_median_drink_duration; median_drink_duration];
                session_mean_stopping_period_hesitation = [session_mean_stopping_period_hesitation; mean_stopping_period_hesitation];
                session_median_stopping_period_hesitation = [session_median_stopping_period_hesitation; median_stopping_period_hesitation];
                session_reward_size = [session_reward_size; Experiment_Information.reward_size];

                load('decoder_candidateEvents.mat')
                decoder_events = decoder_events(event_choice);



                if event_choice == 4
                    % this removes events that were detected by the filtering
                    % method in one map, when the majority of the posterior was
                    % actually in the toher map
                    inds_to_remove = [];
                    for ind = 1:length(decoder_events.replayEvents)
                        if isempty(decoder_events.replayEvents(ind).duration_og)
                            inds_to_remove = [inds_to_remove; ind];
                        end
                    end
                    decoder_events.replayEvents(inds_to_remove) = [];


                    inds_to_remove = find(([decoder_events(sessionNum_decoder).replayEvents.total_posterior_left] > [decoder_events(sessionNum_decoder).replayEvents.total_posterior_right] & [decoder_events(sessionNum_decoder).replayEvents.best_map] == 2)...
                        | ([decoder_events(sessionNum_decoder).replayEvents.total_posterior_left] < [decoder_events(sessionNum_decoder).replayEvents.total_posterior_right] & [decoder_events(sessionNum_decoder).replayEvents.best_map] == 1));

                    decoder_events(sessionNum_decoder).replayEvents(inds_to_remove) = [];
                end



                % get the indices of the candidate events for the session of interest
                t_sub = struct2table(decoder_events(sessionNum_decoder).replayEvents);


                % Flag 19 means this was a saline session. change all
                % 'laser states' to 0. Flag 18 means this was a CNO session.
                % change all 'laser states' to 1.
                if sum(ismember(session_flags,19)>0)
                    t_sub.laser_state(:) = zeros(height(t_sub),1);
                end
                if sum(ismember(session_flags,18)>0)
                    t_sub.laser_state(:) = ones(height(t_sub),1);
                end

                ind_session = ismember(t_sub.timePoints_og,compute_dataTemporalConcatenation(t_sub.timePoints_og,Experiment_Information.Segments(sub_session_num).Times),'rows');

                % t_sub.weighted_r(:) = abs(t_sub.weighted_r);

                % Limit table to events that occured during the session
                % of interest
%                 if iscell(t_sub.dispersion)
%                 t_sub.dispersion(cellfun(@isempty,t_sub.dispersion)) = {nan};
%                 t_sub.dispersion = cell2mat(t_sub.dispersion);
%                 end


                t_sub = t_sub(ind_session==1,:);
                t_sub.rat_label(:) = repmat(rat,height(t_sub),1);
                t_sub.day_label(:) = repmat(day,height(t_sub),1);
                t_sub.segment_label(:) = repmat(sub_session_num,height(t_sub),1);
                t_sub.unique_session_id(:) = repmat(session_id, height(t_sub),1);
                t_sub.session_str(:) = repmat({fullfile(directory,dayFiles{day})},height(t_sub),1);
                t_sub.number_passes_in_session(:) = repmat(num_passes_in_session,height(t_sub),1);
                t_sub.session_flags(:) = repmat({session_flags},height(t_sub),1);
                t_sub.session_percent_correct_directional_assignment = repmat(percent_correct_directional_assigment,height(t_sub),1);
                t_sub.session_median_decoding_error = repmat(median_decoding_error,height(t_sub),1);


                % Make sure t_sub has the same number of fields as
                % data_tbl;
                a = t_sub.Properties.VariableNames;
                b = data_tbl.Properties.VariableNames;

                if length(a)~=length(b) &&  length(b)>0
                    keyboard
                end

                setdiff(b,a)
                % Add this session to the growing table
                data_tbl = [data_tbl; t_sub];
            end

        end
    end
    cd ..
    end
end
%%
% Make a session_table;
% Save a number of session-wide properties

t_session = table;
t_session.rat = session_rat;
t_session.rat_num = session_rat_num;
t_session.session_flags = session_flag_list;
t_session.unique_id = session_unique_session_id;
t_session.label = session_label;
t_session.directional_field_corr = session_directional_field_corr;
t_session.firing_rate_modulation = session_firing_rate_modulation;
t_session.decoding_error = session_decoding_error;
t_session.pcnt_correct_directional_decoding = session_pcnt_correct_directional_decoding;
t_session.num_place_cells = session_num_place_cells;
t_session.ripple_rate = session_ripple_rate;
t_session.ripple_rate_0 = session_ripple_rate_0;
t_session.ripple_rate_1 = session_ripple_rate_1;
t_session.sde_rate = session_sde_rate;
t_session.sde_rate_0 = session_sde_rate_0;
t_session.sde_rate_1 = session_sde_rate_1;
t_session.novel_num  = session_novel_num;
t_session.num_passes = session_num_passes;
t_session.num_passes_2 = session_num_passes_2;
t_session.time_stationary_at_endzone = session_time_stationary_at_endzone;
t_session.time_stationary_at_endzone_0 = session_time_stationary_at_endzone_0;
t_session.time_stationary_at_endzone_1 = session_time_stationary_at_endzone_1;
t_session.mean_stopping_period_duration = session_mean_stopping_period_duration;
t_session.median_stopping_period_duration = session_median_stopping_period_duration;
t_session.mean_drink_duration = session_mean_drink_duration;
t_session.median_drink_duration = session_median_drink_duration;
t_session.mean_stopping_period_hesitation = session_mean_stopping_period_hesitation;
t_session.median_stopping_period_hesitation = session_median_stopping_period_hesitation;
t_session.reward_size = session_reward_size;



