load Experiment_Information
load Analysis_Information
load clusters
load Position_Data
load spikeDensity
load laser_state

if exist('candidateEvents.mat')==2
    load('candidateEvents.mat');
end

Run_Times = Experiment_Information.Run_Times;
Sleep_Times = Experiment_Information.Sleep_Times;
spikeSampRate = Experiment_Information.spikeSampRate;
spatialDim = Experiment_Information.spatialDim;

if spatialDim == 1
    directionalDecoding = 1;
    numSpatialBins = [1 numSpatialBins(2)];
    load binDecoding_02
else
    load binDecoding_08
end

times_list = cat(2, (cat(2, Run_Times{:})), Sleep_Times{:});
Times_day = [min(times_list) max(times_list)];


%candidate event criteria
if spatialDim==1
sequence_spikeDensityThr = -inf;
sequence_posteriorSpreadThr = 20;
% this essentially sets the max jump distance within a replay event;
sequence_delxThr = 0.4*(Experiment_Information.maze_size); % Jump threshold in cm 
sequence_jumpThr = 0.4*(Experiment_Information.maze_size); % Jump threshold in cm 
sequence_deltThr = 0.05;
sequence_durationThr = 0.05; %sec
replay_durationThr = 0.05;
else
    % John's criterion:
% sequence_spikeDensityThr = -inf;
% sequence_posteriorSpreadThr = 10;
% sequence_delxThr = 20;
% sequence_deltThr = .05;
% sequence_jumpThr = 20;
% sequence_durationThr = 0.05; %sec
% replay_durationThr = 0.05;

% Caitlin's old criterion:
%sequence_posteriorSpreadThr = 20;
%sequence_delxThr = 20;
%sequence_jumpThr = 20;
%sequence_posteriorSpreadThr = (Experiment_Information.maze_size^2)*0.0025; % equivalent to John's 20 cm.

sequence_spikeDensityThr = -inf;
sequence_posteriorSpreadThr = (Experiment_Information.maze_size^2)*0.0012; % equivalent to John's 10 cm.
sequence_delxThr = 0.2*(Experiment_Information.maze_size); % Jump threshold in cm 
sequence_jumpThr = 0.2*(Experiment_Information.maze_size); % Jump threshold in cm 
sequence_deltThr = .05;
sequence_durationThr = 0.05; %sec
replay_durationThr = 0.05;
end


speedThr = inf;
speedRangeThr = [0 speedThr];
% speedRangeThr = [20 inf];

%positions
times = load_timeBins_cm(Times_day,decoder_binDecoding(1).shiftSizeDecoding*spikeSampRate,decoder_binDecoding(1).windowSizeDecoding*spikeSampRate);

if size(Position_Data,2)==7 % Caitlin's position data had 7 columns, number 4 being head direction
    circ_columns = [4];
else
    circ_columns = [4 7 9]; % John's position data had 9 columns, with 4 7 and 9 related to head direction
end

Position_Data_full = load_positions_full(Run_Times,Sleep_Times,times,Position_Data,circ_columns);
Position_Data_full(isnan(Position_Data_full(:,5)),5) = 0;

% pull out the start times for each segment of the day- this will be used
% to compute the time into the session that each replay occured.
segment_start_times = [];
for i = 1:length(Experiment_Information.Segments)
    segment_start_times = [segment_start_times; Experiment_Information.Segments(i).Times(1)];
end

decoder_replay = struct;

for sessionNum_decoder = 1:size(Run_Times,1)
    
    %filter bin decoding
    if spatialDim == 2
        [x_NaN,x,timeBins,posteriorSpread,spikeDensity_sub,Position_Data_sub,laser_state_sub,~,posteriorPeak] = compute_filtering_binDecoding_cm(decoder_binDecoding(sessionNum_decoder),...
            sequence_posteriorSpreadThr,speedRangeThr,sequence_spikeDensityThr,Position_Data_full,spikeDensity,spikeDensityStepSize,sequence_jumpThr,Times_day,binSize,spikeSampRate, laser_state);
        
        %extract sequences
        %find sequences bounded by NaNs
        [boundaries,lengths] = compute_allSequences_NaNseparated(x_NaN(:,2));
        
        %merge sequences if gap is delx<delxThr and delt<deltThr
        [boundaries,lengths] = compute_allSequences_NaNseparated_merge(x_NaN,boundaries,sequence_delxThr,sequence_deltThr*spikeSampRate,binSize);
        
        %remove short sequences (< sequence_durationThr)
        ind_remove = find(lengths<sequence_durationThr/decoder_binDecoding(1).shiftSizeDecoding);
        boundaries(ind_remove,:) = [];
        lengths(ind_remove,:) = [];
        
    elseif spatialDim == 1
       
        % for now, creating two structures (one for each map/runing direction) that hold the decoded position
        % for all timepoints using the decoder of interest. Because this is
        % 1D, the second column of the posterior COM is set to all ones.
        decoder_binDecoding_map1 = decoder_binDecoding(sessionNum_decoder);
        decoder_binDecoding_map1.posteriorCOM(:,1) = decoder_binDecoding(sessionNum_decoder).posteriorCOM(:,1);
        decoder_binDecoding_map1.posteriorCOM(:,2) = ones(size(decoder_binDecoding(sessionNum_decoder).posteriorCOM(:,1)));
        decoder_binDecoding_map1.posteriorSpread = decoder_binDecoding(sessionNum_decoder).posteriorSpread(:,1);
        decoder_binDecoding_map1.posteriorPeak = decoder_binDecoding(sessionNum_decoder).posteriorPeak(:,1);
        
        decoder_binDecoding_map2 = decoder_binDecoding(sessionNum_decoder);
        decoder_binDecoding_map2.posteriorCOM(:,1) = decoder_binDecoding(sessionNum_decoder).posteriorCOM(:,2);
        decoder_binDecoding_map2.posteriorCOM(:,2) = ones(size(decoder_binDecoding(sessionNum_decoder).posteriorCOM(:,2)));
        decoder_binDecoding_map2.posteriorSpread = decoder_binDecoding(sessionNum_decoder).posteriorSpread(:,2);
        decoder_binDecoding_map2.posteriorPeak = decoder_binDecoding(sessionNum_decoder).posteriorPeak(:,2);
        
        % this filters the decoded position using John's filtering method
        % (which Caitlin adapted for 1D/directional decoding).
        % For now, each map is filtered separately. The second column of x_NaN_direction_1 is
        % the COM at each decoding timepoint, but with times 
        % when the rat was either moving too fast, or the posteriorSpread
        % too large, set to NaN.

        [x_NaN_direction_1,x_direction_1,~,posteriorSpread_direction_1,~,~,~,~,posteriorPeak_direction_1,jumps_direction_1]...
            = compute_filtering_binDecoding_cm(decoder_binDecoding_map1,...
            sequence_posteriorSpreadThr,speedRangeThr,sequence_spikeDensityThr,Position_Data_full,spikeDensity,spikeDensityStepSize,sequence_jumpThr,Times_day,binSize,spikeSampRate,laser_state);
        [x_NaN_direction_2,x_direction_2,timeBins,posteriorSpread_direction_2,spikeDensity_sub,Position_Data_sub,laser_state_sub,~,posteriorPeak_direction_2,jumps_direction_2]...
            = compute_filtering_binDecoding_cm(decoder_binDecoding_map2,...
            sequence_posteriorSpreadThr,speedRangeThr,sequence_spikeDensityThr,Position_Data_full,spikeDensity,spikeDensityStepSize,sequence_jumpThr,Times_day,binSize,spikeSampRate,laser_state);
        decoding_time_bin_centers = x_direction_1(:,1);

        %find sequences bounded by NaNs:
        % boundaries gives the indices of each segment pulled out in
        % compute_filtering_binDecoding_cm
        [boundaries_direction_1,lengths_direction_1] = compute_allSequences_NaNseparated(x_NaN_direction_1(:,2));
        [boundaries_direction_2,lengths_direction_2] = compute_allSequences_NaNseparated(x_NaN_direction_2(:,2));
        
        %merge sequences if the time gap between them (delt) is less than deltThr, and if the physical jump between them (delx) is below the delxThr 
        % this effectively merges sequences that were interupted (briefly)
        % by poor decoding quality 
        [boundaries_direction_1,lengths_direction_1] = compute_allSequences_NaNseparated_merge(x_NaN_direction_1,boundaries_direction_1,sequence_delxThr,sequence_deltThr*spikeSampRate,binSize);
        [boundaries_direction_2,lengths_direction_2] = compute_allSequences_NaNseparated_merge(x_NaN_direction_2,boundaries_direction_2,sequence_delxThr,sequence_deltThr*spikeSampRate,binSize);
        
        %remove short sequences (< sequence_durationThr)
        ind_remove_direction_1 = find(lengths_direction_1<sequence_durationThr/decoder_binDecoding(1).shiftSizeDecoding);
        boundaries_direction_1(ind_remove_direction_1,:) = [];
        lengths_direction_1(ind_remove_direction_1,:) = [];
        
        ind_remove_direction_2 = find(lengths_direction_2<sequence_durationThr/decoder_binDecoding(1).shiftSizeDecoding);
        boundaries_direction_2(ind_remove_direction_2,:) = [];
        lengths_direction_2(ind_remove_direction_2,:) = [];
     
        % previous version:
%          boundaries_direction_1 = timeBins(boundaries_direction_1);
%          boundaries_direction_2 = timeBins(boundaries_direction_2);
        
%         boundaries_direction_1 = [decoding_time_bin_centers(boundaries_direction_1(:,1)) decoding_time_bin_centers(boundaries_direction_1(:,2))];
%         boundaries_direction_2 = [decoding_time_bin_centers(boundaries_direction_2(:,1)) decoding_time_bin_centers(boundaries_direction_2(:,2))];


%        long_segs = find(lengths_direction_2 > 20);
%         for n = 101:200
%             figure()
%             plot(x_direction_2(boundaries_direction_2(long_segs(n),1):boundaries_direction_2(long_segs(n),2),2))
%             hold on
%             plot(x_NaN_direction_2(boundaries_direction_2(long_segs(n),1):boundaries_direction_2(long_segs(n),2),2))
%         end
        
        %At this point, filtering has been performed separately on the two
        %maps. This means you can sometimes wind up with a segment of
        %smooth decoding in the "wrong" map (i.e., the one with less
        %overall posterior). I take care of this later when pulling out
        %good replays.

    end
    
    Position_Data_sub_scaled = compute_locsToBins(Position_Data_sub,numSpatialBins,x_edges,y_edges);

    if spatialDim == 2
        %compute properties of candidate events and load struct
        for i = 1:size(boundaries,1)
            indSub = boundaries(i,1):boundaries(i,2);
  
            event_timeBins = timeBins(indSub,:);
            event_timepoints = mean(event_timeBins,2);
            
            %compute the time of the replay since the start of the
            %session
            [time_in_session, session_ind] = min(event_timepoints(1) - segment_start_times(event_timepoints(1) - segment_start_times > 0));
            time_in_session = time_in_session/spikeSampRate; % sec
            decoder_replay(sessionNum_decoder).replayEvents(i).time_in_session = time_in_session;
            
            ratPos = [Position_Data_sub_scaled(indSub(1),2),Position_Data_sub_scaled(indSub(1),3)];
            decoder_replay(sessionNum_decoder).replayEvents(i).ratPos = ratPos;
            
            indNaN = find(isnan(x_NaN(indSub,2)));
            replay = x(indSub,2:3);
            replay_NaN = x_NaN(indSub,2:3);
            replay_NaNremoved = replay_NaN; replay_NaNremoved(indNaN,:) = [];
            times_NaNremoved = x(indSub,1); times_NaNremoved(indNaN) = [];
            decoder_replay(sessionNum_decoder).replayEvents(i).indNaN = indNaN;
            decoder_replay(sessionNum_decoder).replayEvents(i).indData = indSub(setdiff(1:length(indSub),indNaN))';
            decoder_replay(sessionNum_decoder).replayEvents(i).timeBins = timeBins(indSub,:);
            decoder_replay(sessionNum_decoder).replayEvents(i).timePoints = [x(indSub(1),1) x(indSub(end),1)];
            decoder_replay(sessionNum_decoder).replayEvents(i).duration = (x(indSub(end),1)-x(indSub(1),1))/spikeSampRate;
            decoder_replay(sessionNum_decoder).replayEvents(i).maxJump_NaN = max(compute_sequenceJumps(replay_NaN));
            decoder_replay(sessionNum_decoder).replayEvents(i).maxJump_NaNremoved = max(compute_sequenceJumps(replay_NaNremoved));
            decoder_replay(sessionNum_decoder).replayEvents(i).maxJump_NaNremoved_time = max(diff(times_NaNremoved));
            decoder_replay(sessionNum_decoder).replayEvents(i).replay = replay;
            decoder_replay(sessionNum_decoder).replayEvents(i).startDistFromRat = sqrt((ratPos(1)-replay(1,1))^2 + (ratPos(2)-replay(1,2))^2);
            decoder_replay(sessionNum_decoder).replayEvents(i).distance = compute_sequenceDistance(replay_NaN);
            decoder_replay(sessionNum_decoder).replayEvents(i).dispersion = compute_sequenceDispersion(replay_NaN);
            decoder_replay(sessionNum_decoder).replayEvents(i).ratSpeed = Position_Data_sub(indSub(1),5);
            decoder_replay(sessionNum_decoder).replayEvents(i).ratHD = Position_Data_sub(indSub(1),4);
            decoder_replay(sessionNum_decoder).replayEvents(i).laser_state = laser_state_sub(indSub,2);
            
        end
        
        save('replayEvents_no_speed_thr.mat','decoder_replay')
        
    elseif spatialDim == 1
        % for linear track sessions, save the boundaries times. Will
        % compute all metrics later in 'load_candidate_events_cm'
        candidateEvents.filtering(sessionNum_decoder).timeBins = timeBins;        
        candidateEvents.filtering(sessionNum_decoder).boundaries_direction_1 = boundaries_direction_1;
        candidateEvents.filtering(sessionNum_decoder).boundaries_direction_2 = boundaries_direction_2;
        candidateEvents.filtering(sessionNum_decoder).x_NaN_direction_1 = x_NaN_direction_1;
        candidateEvents.filtering(sessionNum_decoder).x_NaN_direction_2 = x_NaN_direction_2;   
        candidateEvents.filtering(sessionNum_decoder).x_direction_1 = x_direction_1;
        candidateEvents.filtering(sessionNum_decoder).x_direction_2 = x_direction_2;           
        candidateEvents.filtering(sessionNum_decoder).posteriorSpread_direction_1 = posteriorSpread_direction_1;
        candidateEvents.filtering(sessionNum_decoder).posteriorSpread_direction_2 = posteriorSpread_direction_2;
        candidateEvents.filtering(sessionNum_decoder).posteriorPeak_direction_1 = posteriorPeak_direction_1;
        candidateEvents.filtering(sessionNum_decoder).posteriorPeak_direction_2 = posteriorPeak_direction_2;
               
        
        if exist('candidateEvents_no_speed_thr.mat') == 2
            save('candidateEvents_no_speed_thr.mat','candidateEvents','-append')
        else
            save('candidateEvents_no_speed_thr.mat','candidateEvents')
        end
    end
    
end 
    
    
