function sessions_that_meet_criterion_day = check_for_relevant_sessions(groups_to_compare,num_groups_to_compare,Experiment_Information,must_include_all_flags,segment_flag,novel_range)
%Check if any of the sessions in this day's folder meet the criterion to
%examine:

%If segment_flag is 1, look through the Segments structures to find
%relevant pieces of the session/day. If segment flag is 0, look through the
% Session_Runs structures to find relevant sessions from the day.



% sessions_that_meet_criterion_day = cell(num_groups_to_compare,1);

% for group_to_compare = 1:num_groups_to_compare

% determine which segments from this day contain the appropriate
% flags for each session type to be compared:
%     flags_to_include = groups_to_compare(group_to_compare).flags_to_include;
%     flags_to_exclude = groups_to_compare(group_to_compare).flags_to_exclude;
flags_to_include = groups_to_compare.flags_to_include;
flags_to_exclude = groups_to_compare.flags_to_exclude;

sessions_that_meet_criterion_sub = [];

if segment_flag == 1
    for segment = 1:length(Experiment_Information.Segments)
        meets_criterion = check_session_flags(Experiment_Information.Segments(segment).Flags,flags_to_include,must_include_all_flags,flags_to_exclude,novel_range);
        if meets_criterion == 1
            sessions_that_meet_criterion_sub = [sessions_that_meet_criterion_sub; segment];
            disp(Experiment_Information.Segments(segment).Flags);
        end
    end
elseif segment_flag == 0
    for segment = 1:length(Experiment_Information.Session_Runs)
        meets_criterion = check_session_flags(Experiment_Information.Session_Runs(segment).Flags,flags_to_include,must_include_all_flags,flags_to_exclude,novel_range);
        if meets_criterion == 1
            sessions_that_meet_criterion_sub = [sessions_that_meet_criterion_sub; segment];
            disp(Experiment_Information.Session_Runs(segment).Flags);
        end
    end
end

%    sessions_that_meet_criterion_day{group_to_compare} = sessions_that_meet_criterion_sub;
sessions_that_meet_criterion_day = sessions_that_meet_criterion_sub;

