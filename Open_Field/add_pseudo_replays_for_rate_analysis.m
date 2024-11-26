replay_with_pseudo_replays_for_trial = replay;
[C,ia, ic] = unique(replay_with_pseudo_replays_for_trial(:,{'unique_session_id','drink_period_number','home_event','rat_label','drink_period_time','laser_state_binary'}));
unique_rats_sessions = unique(replay_with_pseudo_replays_for_trial(:,{'rat_label','day_label','segment_label','session_str'}))
sessions = unique_rats_sessions.session_str;
replay_with_pseudo_replays_for_trial.fake_replay_included_just_for_trial_entry(:) = 0;


%%
for n=1:length(sessions)
    segment_label = unique_rats_sessions.segment_label(n);
    session_path=sessions{n};
    session_path = strrep(session_path,'Insync/caitlinmallory@berkeley.edu/Google Drive','Data');
    cd(session_path)
    load Experiment_Information
    run_segments= [];
    for seg=1:length(Experiment_Information.Segments)
        if ismember(11,Experiment_Information.Segments(seg).Flags)
            run_segments = [run_segments; seg];
        end
    end

    run_segment_ind = find(run_segments==segment_label);
    if ismember(unique_rats_sessions.rat_label(n),[1,4,6])
        load Behavior_Analysis.mat
        trials = 1:length(Behavior_Analysis(run_segment_ind ).goal_well_list);
        trial_times = (Behavior_Analysis(run_segment_ind ).entry_exit_times(:,2)-Behavior_Analysis(run_segment_ind ).entry_exit_times(:,1))./30000;
        replay_sub = replay_with_pseudo_replays_for_trial(strcmp(replay_with_pseudo_replays_for_trial.session_str,sessions{n})&replay_with_pseudo_replays_for_trial.segment_label==segment_label,:);
        trials_with_replays = unique(replay_sub.drink_period_number);
        trials_without_replays = setdiff(trials,trials_with_replays);
        new_entries = repmat(replay_sub(end,:),[length(trials_without_replays),1]);
        new_entries.drink_period_number(:) = trials_without_replays;
        new_entries.drink_period_time(:) = trial_times(trials_without_replays);
        home = mode(Behavior_Analysis(run_segment_ind).goal_well_list(:,1));
        for new_entry = 1:length(trials_without_replays)
            if Behavior_Analysis(run_segment_ind ).goal_well_list(trials_without_replays(new_entry),1)==home
                new_entries.home_event(new_entry) = 1;
                new_entries.away_event(new_entry) = 0;
            else
                new_entries.home_event(new_entry) = 0;
                new_entries.away_event(new_entry) = 1;
            end
            new_entries.future(new_entry)=0;
            new_entries.past(new_entry)=0;
            new_entries.fake_replay_included_just_for_trial_entry(new_entry)=1;
            
        end
    elseif ismember(unique_rats_sessions.rat_label(n),[12, 13, 14])
        load mazeTimes
        trials = 1:length([session(run_segment_ind).trials.rewardLoc_ID]);
        trial_times = [session(run_segment_ind).trials.trialDuration];
        home_event = [session(run_segment_ind).trials.home];
        replay_sub = replay_with_pseudo_replays_for_trial(strcmp(replay_with_pseudo_replays_for_trial.session_str,sessions{n})&replay_with_pseudo_replays_for_trial.segment_label==segment_label,:);
        trials_with_replays = unique(replay_sub.drink_period_number);
        trials_without_replays = setdiff(trials,trials_with_replays);
        new_entries = repmat(replay_sub(end,:),[length(trials_without_replays),1]);
        new_entries.drink_period_number(:) = trials_without_replays;
        new_entries.drink_period_time(:) = trial_times(trials_without_replays);
        new_entries.home_event(:) = home_event(trials_without_replays);
        for new_entry = 1:length(trials_without_replays)
            if new_entries.home_event(new_entry)==1
                new_entries.away_event(new_entry) = 0;
            else
                new_entries.away_event(new_entry) = 1;
            end
            new_entries.future(new_entry)=0;
            new_entries.past(new_entry)=0;
            new_entries.fake_replay_included_just_for_trial_entry(new_entry)=1;
        end
    end
    replay_with_pseudo_replays_for_trial = [replay_with_pseudo_replays_for_trial; new_entries];
end
