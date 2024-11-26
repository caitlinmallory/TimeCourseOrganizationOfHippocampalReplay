function meets_criterion = check_session_flags(session_flags,flags_to_include,must_include_all_flags,flags_to_exclude,novel_range)

    meets_criterion = 0;

    desired_flags_met = zeros(length(flags_to_include),1);
    for i = 1:length(desired_flags_met)
        if ismember(flags_to_include(i),session_flags)
            desired_flags_met(i) = 1;
        end
    end

    if ismember(10,flags_to_include)
        % check the range of session numbers that are to be considered in terms
        % of novelty
        [~,index_of_novelty_flag] = find(floor(session_flags) == 10);
        if session_flags(index_of_novelty_flag) >= novel_range(1) && session_flags(index_of_novelty_flag) <= novel_range(2)
            desired_flags_met(flags_to_include==10) = 1;
        end
    end


    if must_include_all_flags == 1 % in this case, the session must contain all the flags in the set flags_to_include
        if sum(desired_flags_met) == length(flags_to_include)
            meets_criterion = 1;
        end
    elseif must_include_all_flags == 0 % % in this case, the session must contain any of the flags in the set flags_to_include
        if sum(desired_flags_met) > 0
                meets_criterion = 1;
        end
    end
    
    % If the session contains any of the flags in the set flags_to_exclude,
    % it does not meet criterion:
    
    if sum(ismember(flags_to_exclude,session_flags)) > 0
        meets_criterion = 0;
    end