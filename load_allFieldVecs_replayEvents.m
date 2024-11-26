replayEvents = [decoder_replay(sessionNum_decoder).replayEvents];

ratPos = load_fieldVec(replayEvents,'ratPos',2);
ratSpeed = load_fieldVec(replayEvents,'ratSpeed',1);
ratHD = load_fieldVec(replayEvents,'ratHD',1);

timePoints = load_fieldVec(replayEvents,'timePoints',2);
duration = load_fieldVec(replayEvents,'duration',1); 
dispersion = load_fieldVec(replayEvents,'dispersion',1); 

replayLaserState = arrayfun(@(x) x.laser_state, replayEvents, 'UniformOutput', false); 
replayLaserStateBinary = (cellfun(@(x) any(x == 1), replayLaserState))';


% distance = load_fieldVec(replayEvents,'distance',1); 
% linCorr = load_fieldVec(replayEvents,'linCorr',1); 
% linearity = load_fieldVec(replayEvents,'linearity',1); 
% curvature = load_fieldVec(replayEvents,'curvature',1); 

% mvl = load_fieldVec(replayEvents,'mvl',1); 
% alpha = load_fieldVec(replayEvents,'alpha',1); 
% diffusionCoef = load_fieldVec(replayEvents,'diffusionCoef',1); 

% maxJump_NaNremoved = load_fieldVec(replayEvents,'maxJump_NaNremoved',1); 

%meanAngDisplacement_futPath = load_fieldVec(replayEvents,'meanAngDisplacement_futPath',1); 
%meanAngDisplacement_pastPath = load_fieldVec(replayEvents,'meanAngDisplacement_pastPath',1); 

%spikeDensityPeak = load_fieldVec(replayEvents,'spikeDensityPeak',1); 
% spikeDensityMean = load_fieldVec(replayEvents,'spikeDensityMean',1); 
% ripplePowerPeak = load_fieldVec(replayEvents,'ripplePowerPeak',1); 
% ripplePowerMean = load_fieldVec(replayEvents,'ripplePowerMean',1); 
% thetaPowerPeak = load_fieldVec(replayEvents,'thetaPowerPeak',1); 
% thetaPowerMean = load_fieldVec(replayEvents,'thetaPowerMean',1); 
% posteriorSpreadMean = load_fieldVec(replayEvents,'posteriorSpreadMean',1);  
% posteriorSpreadMax = load_fieldVec(replayEvents,'posteriorSpreadMax',1);  

%barrier avoidance score
% for i = 1:length(replayEvents)
%     barrierAvoidance_allBarriers = replayEvents(i).barrierAvoidance_allBarriers;
%     if size(barrierAvoidance_allBarriers,1)>size(barrierAvoidance_allBarriers,2)
%     replayEvents(i).barrierAvoidance_allBarriers = replayEvents(i).barrierAvoidance_allBarriers';
%     end
% end
% barrierAvoidance_allBarriers = load_fieldVec(replayEvents,'barrierAvoidance_allBarriers',numConfigs);


% home = load_fieldVec(replayEvents,'home',1);
% trialDuration = load_fieldVec(replayEvents,'trialDuration',1);

% replay_startSector = NaN(length(replayEvents),1);
% replay_endSector = NaN(length(replayEvents),1);
% ratLoc_sector = NaN(length(replayEvents),1);
% numBarrierCrossings = NaN(length(replayEvents),1);
% numBarrierCrossings_randBarriers = NaN(length(replayEvents),1);
% barrierCrossings = NaN(length(replayEvents),1);
% barrierCrossings_full = NaN(length(replayEvents),1);
% barrierCrossings_prev = NaN(length(replayEvents),1);
% barrierCrossings_randBarriers = NaN(length(replayEvents),1);
% barrierCrossings_pValue = NaN(length(replayEvents),1);
% barrierOverlap = NaN(length(replayEvents),1);
% barrierOverlap_full = NaN(length(replayEvents),1);
% barrierOverlap_prev = NaN(length(replayEvents),1);
% barrierOverlap_randBarriers = NaN(length(replayEvents),1);
% barrierOverlap_pValue = NaN(length(replayEvents),1);
% home = NaN(length(replayEvents),1);
% trialDuration = NaN(length(replayEvents),1);
% for i = 1:length(replayEvents)
%     if ~isempty(replayEvents(i).replay_startSector)
%         home(i) = replayEvents(i).home;
%         trialDuration(i) = replayEvents(i).trialDuration;
%         replay_startSector(i) = replayEvents(i).replay_startSector; 
%         replay_endSector(i) = replayEvents(i).replay_endSector;  
%         ratLoc_sector(i) = replayEvents(i).ratLoc_sector;  
% %         barrierCrossings(i) = replayEvents(i).barrierCrossings;
% %         barrierCrossings_full(i) = replayEvents(i).barrierCrossings_full;
% %         barrierCrossings_prev(i) = replayEvents(i).barrierCrossings_prev;
% %         barrierCrossings_randBarriers(i) = replayEvents(i).barrierCrossings_randBarriers;
% %         barrierCrossings_pValue(i) = replayEvents(i).barrierCrossings_pValue;
% %         barrierOverlap(i) = replayEvents(i).barrierOverlap;
% %         barrierOverlap_full(i) = replayEvents(i).barrierOverlap_full;
% %         barrierOverlap_prev(i) = replayEvents(i).barrierOverlap_prev;
% %         barrierOverlap_randBarriers(i) = replayEvents(i).barrierOverlap_randBarriers;
% %         barrierOverlap_pValue(i) = replayEvents(i).barrierOverlap_pValue;
% 
%     end
% end

%extract data from struct 
%     timePoints = load_fieldVec(replayEvents,'timePoints',3);
%         
%     spikeDensityPeak = load_fieldVec(replayEvents,'spikeDensityPeak',1); 
%     spikeDensityMean = load_fieldVec(replayEvents,'spikeDensityMean',1); 
%     ripplePowerPeak = load_fieldVec(replayEvents,'ripplePowerPeak',1); 
%     ripplePowerMean = load_fieldVec(replayEvents,'ripplePowerMean',1); 
%     thetaPowerPeak = load_fieldVec(replayEvents,'thetaPowerPeak',1); 
%     thetaPowerMean = load_fieldVec(replayEvents,'thetaPowerMean',1); 
%     
%     if strcmp(sessionString(1),'Run') == 1
%         ratLoc = load_fieldVec(replayEvents,'ratLoc',2);
%         ratSpeed = load_fieldVec(replayEvents,'ratSpeed',1);
%         ratHD = load_fieldVec(replayEvents,'ratHD',1);
%         ratAccMag = load_fieldVec(replayEvents,'ratAccMag',1);
%         ratAccHD = load_fieldVec(replayEvents,'ratAccHD',1);
% 
% %         if strcmp(Track_Type,'binary tree maze')~=1 && strcmp(Track_Type,'triangle maze')~=1 
% % %             trial = load_fieldVec(replayEvents,'trial',1);
% %             futureWell = load_fieldVec(replayEvents,'futureWell',1);
% %             futureWellTimeDiff = load_fieldVec(replayEvents,'futureWellTimeDiff',1);
% %             futureWellDistance = load_fieldVec(replayEvents,'futureWellDistance',1);
% %             futureWellAngle = load_fieldVec(replayEvents,'futureWellAngle',1);
% %             pastWell = load_fieldVec(replayEvents,'pastWell',1);
% %             pastWellTimeDiff = load_fieldVec(replayEvents,'pastWellTimeDiff',1);
% %             pastWellDistance = load_fieldVec(replayEvents,'pastWellDistance',1);
% %             pastWellAngle = load_fieldVec(replayEvents,'pastWellAngle',1);
% %             nearestWell = load_fieldVec(replayEvents,'nearestWell',1);
% %             nearestWellDistance = load_fieldVec(replayEvents,'nearestWellDistance',1); 
% %             nearestWellAngle = load_fieldVec(replayEvents,'nearestWellAngle',1); 
% %         end
%     end
%     
%     distance = load_fieldVec(replayEvents,'distance',1); 
%     duration = load_fieldVec(replayEvents,'duration',1); 
%     linCorr = load_fieldVec(replayEvents,'linCorr',1); 
%     linearity = load_fieldVec(replayEvents,'linearity',1); 
%     spread = load_fieldVec(replayEvents,'spread',1); 
%     meanAcceleration = load_fieldVec(replayEvents,'meanAcceleration',1); 
%     numDataBins = load_fieldVec(replayEvents,'numDataBins',1); 
%     maxJumpSize = load_fieldVec(replayEvents,'maxJumpSize',1); 
% 
%     %largest jump distance across frames of replay Event:
%     replayJumps = zeros(length(replayEvents),1);
%     for i = 1:length(replayEvents)
%         replayJumps(i) = max(replayEvents(i).replayJumps);
%     end
%     
%     %all jumps, concatenated
%     replayJumpsConcat = [];
%     for i = 1:length(replayEvents)
%         replayJumpsConcat = [replayJumpsConcat;replayEvents(i).replayJumps];
%     end
%     
%     %sectorDecoding
%     if strcmp(Track_Type,'triangle maze')==1
%         replay_sector_decodingCoverage = load_fieldVec(replayEvents,'replay_sector_decodingCoverage',6);
%         replayPosterior_directionalBias = load_fieldVec(replayEvents,'replayPosterior_directionalBias',2);
%         replayPosterior_trackDecodingCoverage = load_fieldVec(replayEvents,'replayPosterior_trackDecodingCoverage',3);
%         replayPosterior_trackDecodingCoverage_dir = load_fieldVec(replayEvents,'replayPosterior_trackDecodingCoverage_dir',6); 
%     end
%     
%     %numClusters
%         numClusters = load_fieldVec(replayEvents,'numClusters',1);
%     
%     %all cluster IDs, concatenated
%     clusterIDsAll = [];    %mean across frames of replay Event:
%     numSpks = zeros(length(replayEvents),1);
%     numFrames = zeros(length(replayEvents),1);
%     for i = 1:length(replayEvents)
%         numSpks(i) = mean(replayEvents(i).numSpks);
%         numFrames(i) = length(replayEvents(i).numSpks);
%     end
%     for i = 1:length(replayEvents)
%         for j = 1:length(replayEvents(i).clusterIDs)
%             clusterIDsAll = [clusterIDsAll;replayEvents(i).clusterIDs{j}];
%         end
%     end
   
%     %matrix of all extended events
%     spikeDensity_extendedEvent_mat = zeros(length(replayEvents),length(replayEvents(10).spikeDensity_extendedEvent));
%     ripplePower_extendedEvent_mat = zeros(length(replayEvents),length(replayEvents(10).ripplePower_extendedEvent));
%     thetaPower_extendedEvent_mat = zeros(length(replayEvents),length(replayEvents(10).thetaPower_extendedEvent));
%     for i = 1:length(replayEvents)
%         if length(replayEvents(i).spikeDensity_extendedEvent) == length(replayEvents(10).spikeDensity_extendedEvent)
%             spikeDensity_extendedEvent_mat(i,:) = replayEvents(i).spikeDensity_extendedEvent;
%             ripplePower_extendedEvent_mat(i,:) = replayEvents(i).ripplePower_extendedEvent;
%             thetaPower_extendedEvent_mat(i,:) = replayEvents(i).thetaPower_extendedEvent;
%         end
%     end


%     boutEndPoint = load_fieldVec(replayEvents,'boutEndPoint',2);
%     boutDuration = load_fieldVec(replayEvents,'boutDuration',1);   
%     boutElapsedFrac = load_fieldVec(replayEvents,'boutElapsedFrac',1);
%     boutRewardWellEndPoint = load_fieldVec(replayEvents,'boutRewardWellEndPoint',2);
%     boutRewardWellDuration = load_fieldVec(replayEvents,'boutRewardWellDuration',1);   
%     boutRewardWellElapsedFrac = load_fieldVec(replayEvents,'boutRewardWellElapsedFrac',1);
%     boutHomeWellEndPoint = load_fieldVec(replayEvents,'boutHomeWellEndPoint',2);
%     boutHomeWellDuration = load_fieldVec(replayEvents,'boutHomeWellDuration',1);   
%     boutHomeWellElapsedFrac = load_fieldVec(replayEvents,'boutHomeWellElapsedFrac',1);
    