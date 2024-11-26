load Experiment_Information
load Analysis_Information
load clusters
load Position_Data
load replayEvents
load binDecoding_08
load laser_state

color_future = [ 0.3922 0.8314 0.0745]; color_past=   [0.6392 0.0118 0.6392];

%load positions
Position_Data_scaled = compute_locsToBins(Position_Data,numSpatialBins,x_edges,y_edges);
%load integrated path length
pathLength = compute_sequenceDistance_cumsum(Position_Data_scaled(:,2:3)); pathLength = [0; pathLength];

%plot spikeDensity/SWR amp
plot_spikeDensityPeak = 0;
replay_durationThr = 0.00;

replay_dispersionThr = 12;
replay_dispersionThr = 5;
        
Run_Times = Experiment_Information.Run_Times;
Sleep_Times = Experiment_Information.Sleep_Times;


spikeSampRate = Experiment_Information.spikeSampRate;
spatialDim = Experiment_Information.spatialDim;

unique_Session_Times = vertcat(Experiment_Information.Segments.Times);

Times_day = [min(unique_Session_Times(:)) max(unique_Session_Times(:))];

times = load_timeBins_cm(Times_day,decoder_binDecoding(1).shiftSizeDecoding*spikeSampRate,decoder_binDecoding(1).windowSizeDecoding*spikeSampRate);
if ~isempty(Sleep_Times)
Position_Data = load_positions_full(Run_Times,Sleep_Times,times,Position_Data,4);
end
Position_Data(isnan(Position_Data(:,5)),5) = 0;

if plot_spikeDensityPeak==1
    %load spike density
    load spikeDensity
    spikeDensity_smoothed = conv(spikeDensity(:,2),setUp_gaussFilt([1 1000],windowSizeDecoding_replay/spikeDensityStepSize),'same');
    spikeDensity = [spikeDensity(:,1),spikeDensity_smoothed];
    
    %load LFPs
    
    ripple_tetrode = Experiment_Information.ripple_refs;
    lfp_file_denoised = ['LFP_Data' num2str(ripple_tetrode) '_denoised.mat'];
    lfp_file = ['LFP_Data' num2str(ripple_tetrode) '.mat'];
    
    load(lfp_file_denoised);
    load(lfp_file);
    
    LFP = LFP_Data;
    clear LFP_data
    LFPSampRate = 1500;
    
    [~,LFP_filt_SWRAmp,~] = compute_filteredLFP(SWRFreqRange,LFP(:,2),LFPSampRate);
    LFP = [LFP,LFP_filt_SWRAmp,LFP_filt_SWRAmp];
    
    LFP_smoothed = conv(LFP(:,4),setUp_gaussFilt([1 1000],decoder_binDecoding(1).windowSizeDecoding/(1/LFPSampRate)),'same');
    LFP(:,4) = LFP_smoothed;
end
%%
for sessionNum = 1
% for sessionNum = 1:length(unique_Session_Times)
   
    sessionNum_decoder = Experiment_Information.Segments(sessionNum).Decoder;

    
    figure_page_count = 0;
    
    if plot_spikeDensityPeak==1
        % position:
        Position_Data_sub = compute_dataTemporalConcatenation(Position_Data,load_timeBounds(unique_Session_Times(sessionNum,:)));
        
        %load spike density within session and zscore
        spikeDensity_sub = compute_dataTemporalConcatenation(spikeDensity,load_timeBounds(unique_Session_Times(sessionNum,:)));
        Position_Data_sub = compute_dataInterpolation(Position_Data_sub,spikeDensity_sub(:,1),[]);
        moving_ind = find(Position_Data_sub(:,5) > speedThr);
        spikeDensity_sub(moving_ind,2) = nan;
        spikeDensity_sub(:,2) = compute_zscore(spikeDensity_sub(:,2));
        
        %load LFP ripple amp within session and zscore
        LFP_sub = compute_dataTemporalConcatenation(LFP,load_timeBounds(unique_Session_Times(sessionNum,:)));
        LFP_Data_denoised_sub = compute_dataTemporalConcatenation(LFP_Data_denoised,load_timeBounds(unique_Session_Times(sessionNum,:)));
        
        % position:
        Position_Data_sub = compute_dataInterpolation(Position_Data_sub, LFP_sub(:,1), []);
        
        %limit to stationary periods and mask out artifacts prior to z-scoring:
        artifact_ind = find(isnan(LFP_Data_denoised_sub(:,2)));
        moving_ind = find(Position_Data_sub(:,5) > speedThr);
        bad_inds = unique([artifact_ind; moving_ind]);
        
        LFP_sub(bad_inds,4) = nan;
        LFP_sub(:,4) = compute_zscore(LFP_sub(:,4));
        
    end
    
    %load Position_Data 
        [~,Position_Data_sub] = compute_dataTemporalConcatenation(Position_Data,unique_Session_Times(sessionNum,:));
        Position_Data_sub = compute_locsToBins(Position_Data_sub,numSpatialBins,x_edges,y_edges);
        
    %load decoder
        [ind_cluster,rateMap_smoothed,~,rateMap_smoothed_NaN] = load_rateMapsForDecoding_2D_cm(clusters,sessionNum_decoder,Run_Times,numSpatialBins,meanFiringRateThr,spatialInfoThr,noiseOverlapThr,isolationThr,peakSnrThr,numSpatialBins_coarse);
        rateMap_sub = rateMap_smoothed;
        
    %filter sequences for replays

        load_ind_replayEvents_all_session_types
        ind_sorted = sortrows([dispersion(ind_replay) ind_replay],'descend');
           
    %plot replay mosiaic
        numPanels = 50;
        subplot_xDim = 10;
        subplot_yDim = (ceil((numPanels)/subplot_xDim));
        for i = 1:length(ind_replay)

            %for plotting moaic
            if mod(i,numPanels)==1
                ha = setUp_tightSubplot(subplot_yDim,subplot_xDim);
                %set(gcf,'renderer','Painters')
            end
            
            %load replay properties
                indNaN = decoder_replay(sessionNum_decoder).replayEvents(ind_replay(i)).indNaN;
                replay = decoder_replay(sessionNum_decoder).replayEvents(ind_replay(i)).replay;
                replay_NaN = replay; replay_NaN(indNaN,:) = NaN;
                replay_NaNremoved = replay; replay_NaNremoved(indNaN,:) = [];
                timeBins = decoder_replay(sessionNum_decoder).replayEvents(ind_replay(i)).timeBins;
                timePoints = decoder_replay(sessionNum_decoder).replayEvents(ind_replay(i)).timePoints;
                                
                maxJump_NaN = decoder_replay(sessionNum_decoder).replayEvents(ind_replay(i)).maxJump_NaN;
                maxJump_NaNremoved = decoder_replay(sessionNum_decoder).replayEvents(ind_replay(i)).maxJump_NaNremoved;
                maxJump_NaNremoved_time = decoder_replay(sessionNum_decoder).replayEvents(ind_replay(i)).maxJump_NaNremoved_time;

                ratPos = decoder_replay(sessionNum_decoder).replayEvents(ind_replay(i)).ratPos;
                ratSpeed = decoder_replay(sessionNum_decoder).replayEvents(ind_replay(i)).ratSpeed;
                ratHD = decoder_replay(sessionNum_decoder).replayEvents(ind_replay(i)).ratHD;
                replay_laser_state = decoder_replay(sessionNum_decoder).replayEvents(ind_replay(i)).laser_state;
                dispersion = decoder_replay(sessionNum_decoder).replayEvents(ind_replay(i)).dispersion;

             %   time_since_stopping = decoder_replay(sessionNum_decoder).replayEvents(ind_replay(i)).time_since_stopping_period_onset;
              %  angular_distance_future = decoder_replay(sessionNum_decoder).replayEvents(ind_replay(i)).meanAngDisplacement_futPath;
               % angular_distance_past = decoder_replay(sessionNum_decoder).replayEvents(ind_replay(i)).meanAngDisplacement_pastPath;
%                 angular_distance_future(isnan(angular_distance_future)) = [];
%                 angular_distance_past(isnan(angular_distance_past)) = [];
angular_distance_future = nan;
angular_distance_past = nan;
                 
                laser_on_ind = find([0; diff(replay_laser_state)] == 1);
                if replay_laser_state(1) == 1
                    laser_on_ind = [1; laser_on_ind];
                end
                laser_off_ind = find([0; diff(replay_laser_state)] == -1);
                
                %alpha = decoder_replay(sessionNum_decoder).replayEvents(ind_replay(i)).alpha;
                %speed = decoder_replay(sessionNum_decoder).replayEvents(ind_replay(i)).diffusionCoef;
                
                %linearity = decoder_replay(sessionNum_decoder).replayEvents(ind_replay(i)).linearity;
                %linCorr = decoder_replay(sessionNum_decoder).replayEvents(ind_replay(i)).linCorr;
                             
                if plot_spikeDensityPeak==1
                    spikeDensity_replay = compute_dataTemporalConcatenation(spikeDensity_sub,timePoints);
                    ripplePower_replay = compute_dataTemporalConcatenation(LFP_sub,timePoints);

                    spikeDensityPeak = max(spikeDensity_replay(:,end));
                    ripplePowerPeak = max(ripplePower_replay(:,end));
                end
                                
%                 spikeDensityPeak = decoder_replay(sessionNum_decoder).replayEvents(ind_replay(i)).spikeDensityPeak;
%                 spikeDensityMean = decoder_replay(sessionNum_decoder).replayEvents(ind_replay(i)).spikeDensityMean;
%                 ripplePowerPeak = decoder_replay(sessionNum_decoder).replayEvents(ind_replay(i)).ripplePowerPeak;
%                 ripplePowerMean = decoder_replay(sessionNum_decoder).replayEvents(ind_replay(i)).ripplePowerMean;
%                 posteriorSpreadMean = decoder_replay(sessionNum_decoder).replayEvents(ind_replay(i)).posteriorSpreadMean;
%                 posteriorSpreadMax = decoder_replay(sessionNum_decoder).replayEvents(ind_replay(i)).posteriorSpreadMax;

%                 plot(spikeDensity_replay(:,1),spikeDensity_replay(:,2))
%                 hold on, plot(ripplePower_replay(:,1),ripplePower_replay(:,4)), hold off

% past/future path
% Brad methods: look at position data in the past 10 seconds, or 50 cm,
% whichever is greater.
  % Find the past path that will be considered:
        % pastPath will be the path segment beginning 10 seconds prior, or 50 cm prior- whichever leads
        % to the larger path.
        [~,currentInd] = min(abs(Position_Data_scaled(:,1)-mean(timePoints)));

        pastInd_dist = find(pathLength(currentInd)-pathLength >= 50,1,'last');
        [val,pastInd_time] = min(abs(Position_Data_scaled(:,1) - (mean(timePoints) - 10*30000)));
        pastInd_time_distance_traveled = pathLength(currentInd)-pathLength(pastInd_time);
        if isempty(pastInd_dist)
            pastInd = pastInd_time;
        elseif pastInd_time_distance_traveled >= 50
            pastInd = pastInd_time;
        else
            pastInd = pastInd_dist;
        end
        pastPath = Position_Data_scaled(pastInd:currentInd,:);
        pastPath = pastPath(end:-1:1,:); % reverse the direction so that the first index of past path is the most recent timepoint.
        pastPathLength = pathLength(currentInd)-pathLength(pastInd);
%         if pastPathLength < 50
%             keyboard
%         end

        % Find the furture path that will be considered:
        % futurePath will be the path segment beginning 10 seconds after the current timepoint, or 50 cm ahead- whichever leads
        % to the larger path.
        futureInd_dist = find(pathLength - pathLength(currentInd) >= 50,1,'first');
        [val,futureInd_time] = min(abs(Position_Data_scaled(:,1) - (mean(timePoints) + 10*30000)));
        futureInd_time_distance_traveled = pathLength(futureInd_time) - pathLength(currentInd);
        if isempty(futureInd_dist)
            futureInd = futureInd_time;
        elseif futureInd_time_distance_traveled >= 50
            futureInd = futureInd_time;
        else
            futureInd = futureInd_dist;
        end
        futurePath = Position_Data_scaled(currentInd:futureInd,:);
        futurePathLength = pathLength(futureInd)-pathLength(currentInd);
%         if futurePathLength < 50
%             keyboard
%         end


            %replay decoding
                numSpks = load_numSpks_timeBins(timeBins,clusters,ind_cluster,decoder_binDecoding(sessionNum_decoder).shiftSizeDecoding*spikeSampRate,decoder_binDecoding(sessionNum_decoder).windowSizeDecoding*spikeSampRate);
                [posterior,posteriorCOM,posteriorSpread,posteriorPeak] =  compute_BayesianDecoding_fast(numSpks,rateMap_sub,numSpatialBins,spatialDim,decoder_binDecoding(sessionNum_decoder).windowSizeDecoding,1);
%                     jumps = compute_sequenceJumps(posteriorCOM); jumps = [jumps; jumps(end)];
%                     indNaN = find(posteriorSpread>sequence_posteriorSpreadThr);% & jumps>sequence_jumpThr);
                    posterior(indNaN,:) = NaN;
                    posteriorCOM(indNaN,:) = NaN;
%                     posteriorCOM_NaNremoved = posteriorCOM; posteriorCOM_NaNremoved(indNaN,:) = [];
                [~,posterior_color,posterior_sum,posterior_sum_binarized] = compute_flattenedDecoding(posterior,spatialDim); 
                posterior_color = reshape(posterior_color',numSpatialBins(1),numSpatialBins(2),3);
                
                               
            %plot mosaic
                set(gcf,'color','w')
                if mod(i,numPanels)~=0
                    axes(ha(mod(i,numPanels)));
                else
                    axes(ha(numPanels));
                end
                
                image(posterior_color,'alphaData',sum(posterior_color,3)~=3), set(gca,'ydir','normal')
                hold on, plot(replay_NaNremoved(:,1),replay_NaNremoved(:,2),'color',0.5*ones(1,3),'linewidth',1), hold off
                hold on, plot(replay_NaN(:,1),replay_NaN(:,2),'k','linewidth',1), hold off
                
                hold on
                if ~isempty(laser_on_ind)
                    hold on, plot(replay(laser_on_ind,1),replay(laser_on_ind,2),'or','markerFaceColor','r','markersize', 6,'linewidth',1), hold off
                end
                
                if ~isempty(laser_off_ind)
                    hold on, plot(replay(laser_off_ind,1),replay(laser_off_ind,2),'ok','markerFaceColor','k','markersize',6,'linewidth',1), hold off
                end

                plot_ratLoc
                plot_mazeProperties_cm
                if ~isempty(futurePath), hold on, p1 = plot(futurePath(:,2),futurePath(:,3),'color',color_future,'linewidth',3); p1.Color(4) = 0.8; hold off, end
                if ~isempty(pastPath), hold on, p2 = plot(pastPath(:,2),pastPath(:,3),'color',color_past,'linewidth',3); p2.Color(4) = 0.8; hold off, end
                
                time = {floor(((timePoints(1)/spikeSampRate)-(unique_Session_Times(sessionNum,1)/spikeSampRate))/60),round(rem(((timePoints(1)/spikeSampRate)-(unique_Session_Times(sessionNum,1)/spikeSampRate))/60,1)*60)};

                if time{2}>=10    
                    time_string = strcat(num2str(time{1}),':',num2str(time{2}));
                else
                    time_string = strcat(num2str(time{1}),':0',num2str(time{2}));
                end
                text(1,numSpatialBins(1),time_string,'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Color','K','fontsize',12)
                text(1,numSpatialBins(1)-4,num2str(dispersion,2),'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Color','K','fontsize',12)
%                 text(1,numSpatialBins(1)-4,num2str(angular_distance_future,2),'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Color','K','fontsize',12)
%                 text(1,numSpatialBins(1)-8,num2str(angular_distance_past,2),'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Color','K','fontsize',12)

                
                %                 text(numSpatialBins(2),numSpatialBins(1),num2str(compute_round(duration(ind_replay(i)),100)),'HorizontalAlignment', 'right', 'VerticalAlignment', 'top','Color','K','fontsize',6)
%                 text(numSpatialBins(2),numSpatialBins(1)-4,num2str(compute_round(binSize*dispersion(ind_replay(i)),10)),'HorizontalAlignment', 'right', 'VerticalAlignment', 'top','Color','K','fontsize',6)
%                 text(numSpatialBins(2),numSpatialBins(1)-8,num2str(compute_round(ratSpeed,10)),'HorizontalAlignment', 'right', 'VerticalAlignment', 'top','Color','K','fontsize',1)
%                 
                if plot_spikeDensityPeak==1
                    text(numSpatialBins(2),5,num2str(compute_round(spikeDensityPeak,10)),'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','Color','K','fontsize',6)
                    text(numSpatialBins(2),1,num2str(compute_round(ripplePowerPeak,10)),'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','Color','K','fontsize',6)
                end
                
            %drawnow

             if mod(i,50) == 0 || i == length(ind_replay)
          
                 figure_page_count = figure_page_count + 1;
                 figure_name = ['Decoder' num2str(sessionNum_decoder) '_Session' num2str(sessionNum) '_page' num2str(figure_page_count)];
                 %saveas(gcf,figure_name,'tiff')
    
                 export_fig(figure_name,'-jpeg')

                 saveas(gcf,figure_name)
                 %print(gcf,figure_name,'-tiff','-bestfit');
             end
                
        end
end




