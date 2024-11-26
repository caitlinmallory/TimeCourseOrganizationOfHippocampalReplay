function [ind_replay_laser_off, ind_replay_laser_on, ind_replay] = load_ind_replayEvents_cm(replay_dispersionThr,replay_durationThr,binSize,decoder_replay,session_times,sessionNum_decoder)

[ratLoc, ratSpeed, ratHD, timePoints, duration, dispersion, replayLaserState, replayLaserStateBinary] = load_allFieldVecs_replayEvents_cm(decoder_replay(sessionNum_decoder).replayEvents);

ind_dispersion = (dispersion>(replay_dispersionThr/binSize));
ind_duration = (duration>replay_durationThr);
ind_session = ismember(timePoints,compute_dataTemporalConcatenation(timePoints,load_timeBounds(session_times)),'rows');
ind_laser_on = replayLaserStateBinary == 1;
ind_laser_off = replayLaserStateBinary == 0;


ind_replay = find(ind_duration & ind_dispersion & ind_session);
ind_replay_laser_off = find(ind_duration & ind_dispersion & ind_session & ind_laser_off);
ind_replay_laser_on = find(ind_duration & ind_dispersion & ind_session & ind_laser_on);

