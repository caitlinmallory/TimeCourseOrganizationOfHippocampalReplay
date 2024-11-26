load_allFieldVecs_replayEvents

ind_dispersion = (dispersion>(replay_dispersionThr/binSize));
ind_duration = (duration>replay_durationThr);
ind_session = ismember(timePoints,compute_dataTemporalConcatenation(timePoints,load_timeBounds(Experiment_Information.Run_Times{sessionNum})),'rows');

ind_replay = find(ind_duration==1 & ind_dispersion==1 & ind_session==1);