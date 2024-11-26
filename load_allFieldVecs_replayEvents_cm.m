function [ratPos, ratSpeed, ratHD, timePoints, duration, dispersion, replayLaserState, replayLaserStateBinary] = load_allFieldVecs_replayEvents_cm(replayEvents)
% 
% % replayEvents = [decoder_events(sessionNum_decoder).replayEvents];
try
ratPos = load_fieldVec(replayEvents,'ratPos',2);
catch
    ratPos = load_fieldVec(replayEvents,'ratLoc',2);
end
ratSpeed = load_fieldVec(replayEvents,'ratSpeed',1);
ratHD = load_fieldVec(replayEvents,'ratHD',1);

timePoints = load_fieldVec(replayEvents,'timePoints',2);
duration = load_fieldVec(replayEvents,'duration',1); 
dispersion = load_fieldVec(replayEvents,'dispersion',1); 

replayLaserState = arrayfun(@(x) x.laser_state, replayEvents, 'UniformOutput', false); 
replayLaserStateBinary = (cellfun(@(x) any(x == 1), replayLaserState))';

