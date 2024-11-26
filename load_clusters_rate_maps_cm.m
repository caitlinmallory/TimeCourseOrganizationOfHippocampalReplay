clearvars -except dayFiles day directory rat windows hand_clustered_only

load Experiment_Information
load Analysis_Information
if exist('clusters.mat')==2
load clusters
else
    return
end
load Position_Data

Minimum_Place_Field_Peak_Firing_Rate = 1;
Minimum_Place_Field_Firing_Rate_Fraction=0.2;  % Brad set this to 0.2.
Minimum_Contiguous_Place_Field_Bins = 15; % Brad had this set to 20.
% search_dist = 2.5;
% search_minPts = 50;
% num_reSamples = 2500;

if  Experiment_Information.spatialDim == 1
    edges_spatialBins = {x_edges};
    numSpatialBins = length(x_edges)-1;
    numSpatialBins_coarse = length(x_edges_coarse)-1;
end

spatialDim = Experiment_Information.spatialDim;
posSampRate = Experiment_Information.posSampRate;
spikeSampRate = Experiment_Information.spikeSampRate;
Run_Times = Experiment_Information.Run_Times;
trialsPerRun = Experiment_Information.trialsPerRun;

for i = 1:length(clusters)
    [i,length(clusters)]
    
    for j = 1:length(Experiment_Information.Run_Times)
        
        if isempty(clusters(i).spkTime)
            clusters(i).runs(j).rateMap = nan(numSpatialBins);
            clusters(i).runs(j).rateMap_vectorized = reshape(nan(numSpatialBins),1,prod(numSpatialBins));
            clusters(i).runs(j).rateMap_smoothed = nan(numSpatialBins);
            clusters(i).runs(j).rateMap_smoothed_vectorized = reshape(nan(numSpatialBins),1,prod(numSpatialBins));
            clusters(i).runs(j).rateMap_smoothed_NaN = nan(numSpatialBins);
            clusters(i).runs(j).rateMap_smoothed_NaN_vectorized = reshape(nan(numSpatialBins),1,prod(numSpatialBins));
            
            clusters(i).runs(j).rateMap_coarse = nan(numSpatialBins_coarse);
            clusters(i).runs(j).rateMap_coarse_vectorized = reshape(nan(numSpatialBins_coarse),1,prod(numSpatialBins_coarse));
            clusters(i).runs(j).rateMap_coarse_smoothed = nan(numSpatialBins_coarse);
            clusters(i).runs(j).rateMap_coarse_smoothed_vectorized = reshape(nan(numSpatialBins_coarse),1,prod(numSpatialBins_coarse));
            clusters(i).runs(j).rateMap_coarse_smoothed_NaN = nan(numSpatialBins_coarse);
            clusters(i).runs(j).rateMap_coarse_smoothed_NaN_vectorized = reshape(nan(numSpatialBins_coarse),1,prod(numSpatialBins_coarse));
            
            clusters(i).runs(j).meanFiringRate = nan;
            clusters(i).runs(j).spatialInfo = nan;
            clusters(i).runs(j).spatialCoherence = nan;
            clusters(i).runs(j).spatialSpread = nan;
            clusters(i).runs(j).maxSpatialFiringRate = nan;
            clusters(i).runs(j).meanSpatialFiringRate = nan;
            clusters(i).runs(j).maxSpatialFiringBinLoc = [nan nan];
            clusters(i).runs(j).COMSpatialFiringBinLoc = [nan nan];
            clusters(i).runs(j).speedTuningCurve = nan;
            
            clusters(i).runs(j).num_subFields = nan;
            
            continue
        end
        
        times = Run_Times{j};
        Position_Data_sub = compute_dataTemporalConcatenation(Position_Data,times);
        
        spkData = compute_dataTemporalConcatenation([clusters(i).spkTime,clusters(i).spkPos,clusters(i).spkHD],times);
        spkSpeed = compute_dataTemporalConcatenation([clusters(i).spkTime,clusters(i).spkSpeed],times); spkSpeed = spkSpeed(:,2);
        %spkData = spkData(spkSpeed>speedThr,:);
        spkTime_filt = spkData(spkSpeed>speedThr,1);
        spkPos_filt = spkData(spkSpeed>speedThr,2:2+spatialDim-1);
        spkHD_filt = spkData(spkSpeed>speedThr,3);
        
        %compute spatial tuning curve
        [rateMap, rateMap_smoothed, rateMap_smoothed_NaN, posDensity, ~] = compute_rateMap(spkPos_filt,Position_Data_sub,numSpatialBins,numSpatialBins_smoothing,...
            posSampRate,spatialDim,{x_edges,y_edges},speedThr,binSize);
        [rateMap_coarse, rateMap_coarse_smoothed, rateMap_coarse_smoothed_NaN, ~, ~] = compute_rateMap(spkPos_filt,Position_Data_sub,numSpatialBins_coarse,...
            numSpatialBins_smoothing,posSampRate,spatialDim,{x_edges_coarse,y_edges_coarse},speedThr,binSize_coarse);
        
        %compute spatial spread
        spatialSpread = compute_imageSpread(rateMap_smoothed,2);
        
        %compute spatial info
        spatialInfo = compute_spatialInfo(rateMap_smoothed,posDensity);
        
        %compute spatial coherence
        % SLOW: commenting out for now.
        %             spatialCoherenceNumNearestBins = 10;
        %             spatialCoherence = compute_spatialCoherence(rateMap_smoothed,spatialCoherenceNumNearestBins);
        
        %compute mean firing rate
        meanFiringRate = size(spkPos_filt,1)/sum(diff(times')/spikeSampRate);
        
        %compute mean/max spatial firing rate
        maxSpatialFiringRate = max(rateMap_smoothed(:));
        meanSpatialFiringRate = nanmean(rateMap_smoothed(:));
        
        %compute max firing location
        [maxSpatialFiringBinLoc,COMSpatialFiringBinLoc] = compute_maxSpatialFiringLoc(rateMap_smoothed,spatialDim);
        
        %compute HD tuning
        [rateMap_HD, rateMap_HD_smoothed, rateMap_HD_smoothed_NaN, HDTuningPref, HDTuningSpec, angle_centers] = compute_HDTuning(spkHD_filt,Position_Data_sub,speedThr,posSampRate,spatialDim);
        
        %compute speed tuning
        [speedTuningCurve, speedBins] = compute_speedTuning(times,clusters(i).spkTime,clusters(i).spkSpeed,Position_Data_sub,posSampRate,speedThr);
        
        %number of subfields-- only works for 2D- John's code. Using Brad's
        %method (below) instead for now, as it works for both 1D and 2D.
        %         if spatialDim == 2
        %             if ~isempty(find(~isnan(rateMap_smoothed)))
        %                 [S_fields,~,S_COMs] = compute_clusters_densityBasedClusters(rateMap_smoothed,num_reSamples,search_dist,search_minPts,numSpatialBins,spatialDim,numSpatialBins_smoothing,binSize);
        %                 num_subFields = size(S_fields,1);
        %             else
        %                 S_fields = NaN; S_COMs = NaN;
        %                 num_subFields = NaN;
        %             end
        %         end
        
        % Caitlin added the following, based off Brad Pfeifer's code
        % 'IFRS_CALCULATE_UNIMODAL_AND_BIMODAL_PLACE_FIELD_PROPERTIES
        
        if maxSpatialFiringRate >= Minimum_Place_Field_Peak_Firing_Rate
            Place_Field = double(rateMap_smoothed>=maxSpatialFiringRate*Minimum_Place_Field_Firing_Rate_Fraction);
            num_subFields=0;
            Contiguous_Place_Fields=zeros(size(Place_Field,1),size(Place_Field,2));
            Place_Field_Sizes=0;
            [Y_Index,X_Index]=find(Place_Field==1);
            for M=1:length(Y_Index)
                if Contiguous_Place_Fields(Y_Index(M),X_Index(M))==0 && Place_Field(Y_Index(M),X_Index(M))==1
                    This_Place_Field=grayconnected(Place_Field,Y_Index(M),X_Index(M),0);
                    if sum(sum(This_Place_Field))>=Minimum_Contiguous_Place_Field_Bins
                        num_subFields=num_subFields+1;
                        Contiguous_Place_Fields(This_Place_Field)=num_subFields;
                        Place_Field_Sizes=[Place_Field_Sizes;sum(sum(This_Place_Field))];
                    end
                end
                clear This_Place_Field;
            end
            [Y,X]=find(Contiguous_Place_Fields>0);
            Mean_InField_Firing_Rate=mean(mean(rateMap_smoothed(Y,X)));
            clear M;
            clear X;
            clear Y;
            clear X_Index;
            clear Y_Index;
            if length(Place_Field_Sizes)>1
                Place_Field_Sizes=Place_Field_Sizes(2:end);
            end
            
        else
            Place_Field_Sizes = nan;
            Mean_InField_Firing_Rate = nan;
            num_subFields = 0;
            Place_Field = nan;
            Contiguous_Place_Fields = nan;
        end
        
        
        
        %compute 'burstiness' measure: fraction of spikes whos nearest
        %spike on either side is less than 6 ms (see Wang and Pfeiffer 2021).
        clusters(i).runs(j).burstiness = compute_burstiness(spkData(:,1),spikeSampRate);
        
        
        %compute ISI distribution
        %             edges = linspace(0,0.05,500);
        %             h = histc(diff(spkTime),edges);
        
        %load into struct
        clusters(i).runs(j).rateMap = rateMap;
        clusters(i).runs(j).rateMap_vectorized = reshape(rateMap,1,prod(numSpatialBins));
        clusters(i).runs(j).rateMap_smoothed = rateMap_smoothed;
        clusters(i).runs(j).rateMap_smoothed_vectorized = reshape(rateMap_smoothed,1,prod(numSpatialBins));
        clusters(i).runs(j).rateMap_smoothed_NaN = rateMap_smoothed_NaN;
        clusters(i).runs(j).rateMap_smoothed_NaN_vectorized = reshape(rateMap_smoothed_NaN,1,prod(numSpatialBins));
        
        clusters(i).runs(j).rateMap_coarse = rateMap_coarse;
        clusters(i).runs(j).rateMap_coarse_vectorized = reshape(rateMap_coarse,1,prod(numSpatialBins_coarse));
        clusters(i).runs(j).rateMap_coarse_smoothed = rateMap_coarse_smoothed;
        clusters(i).runs(j).rateMap_coarse_smoothed_vectorized = reshape(rateMap_coarse_smoothed,1,prod(numSpatialBins_coarse));
        clusters(i).runs(j).rateMap_coarse_smoothed_NaN = rateMap_coarse_smoothed_NaN;
        clusters(i).runs(j).rateMap_coarse_smoothed_NaN_vectorized = reshape(rateMap_coarse_smoothed_NaN,1,prod(numSpatialBins_coarse));
        
        clusters(i).runs(j).meanFiringRate = meanFiringRate;
        clusters(i).runs(j).spatialInfo = spatialInfo;
        %clusters(i).runs(j).spatialCoherence = spatialCoherence;
        clusters(i).runs(j).spatialSpread = spatialSpread;
        clusters(i).runs(j).maxSpatialFiringRate = maxSpatialFiringRate;
        clusters(i).runs(j).meanSpatialFiringRate = meanSpatialFiringRate;
        clusters(i).runs(j).maxSpatialFiringBinLoc = maxSpatialFiringBinLoc;
        
        clusters(i).runs(j).COMSpatialFiringBinLoc = COMSpatialFiringBinLoc;
        clusters(i).runs(j).speedTuningCurve = speedTuningCurve;
        
        clusters(i).runs(j).rateMap_HD = rateMap_HD;
        clusters(i).runs(j).rateMap_HD_smoothed = rateMap_HD_smoothed;
        clusters(i).runs(j).rateMap_HD_smoothed_NaN = rateMap_HD_smoothed_NaN;
        
        clusters(i).runs(j).Place_Field_Sizes = Place_Field_Sizes;
        clusters(i).runs(j).Mean_InField_Firing_Rate = Mean_InField_Firing_Rate;
        clusters(i).runs(j).num_subFields = num_subFields;
        clusters(i).runs(j).Place_Field = Contiguous_Place_Fields;
        

        if spatialDim == 1
            %compute directional spatial tuning curves across trials
            
            directional_place_field_sizes = nan;
            directional_mean_in_field_firing_rate = [nan nan];
            directional_mean_firing_rate = [nan nan];
            directional_peak_firing_rate = [nan nan];
            
            for k = 1:2 %direction; 1=left, 2=right
                times = Run_Times{j};
                
                Position_Data_sub = compute_dataTemporalConcatenation(Position_Data,times);
                spkTime = compute_dataTemporalConcatenation(clusters(i).spkTime,times);
                
                spkPosition_Data = compute_dataInterpolation(Position_Data_sub,spkTime,[]);
                spkDir = spkPosition_Data(:,6);
                
                spkPos = compute_dataTemporalConcatenation([clusters(i).spkTime,clusters(i).spkPos],times); spkPos = spkPos(:,2:2+spatialDim-1);
                spkSpeed = compute_dataTemporalConcatenation([clusters(i).spkTime,clusters(i).spkSpeed],times); spkSpeed = spkSpeed(:,2);
                spkHD = compute_dataTemporalConcatenation([clusters(i).spkTime,clusters(i).spkHD],times); spkHD = spkHD(:,2:end);
                
                if k==1
                    spkPos_filt = spkPos(spkSpeed>speedThr & spkDir < 0,:);
                    Position_Data_sub = Position_Data_sub(Position_Data_sub(:,6) < 0,:);
                else
                    spkPos_filt = spkPos(spkSpeed>speedThr & spkDir > 0,:);
                    Position_Data_sub = Position_Data_sub(Position_Data_sub(:,6) > 0,:);
                end
                
                [rateMap, rateMap_smoothed, rateMap_smoothed_NaN, posDensity, ~] = compute_rateMap(spkPos_filt,Position_Data_sub,numSpatialBins,numSpatialBins_smoothing,posSampRate,spatialDim,{x_edges,y_edges},speedThr,binSize);
                
                
                maxSpatialFiringRate = max(rateMap_smoothed(:));
                meanSpatialFiringRate = mean(rateMap_smoothed(:));
                [maxSpatialFiringBinLoc,COMSpatialFiringBinLoc] = compute_maxSpatialFiringLoc(rateMap_smoothed,spatialDim);
                meanFiringRate = size(spkPos_filt,1)/sum(diff(times')/spikeSampRate);
                
                % Caitlin added the following, based off Brad Pfeifer's code
                % 'IFRS_CALCULATE_UNIMODAL_AND_BIMODAL_PLACE_FIELD_PROPERTIES
                
                if maxSpatialFiringRate >= Minimum_Place_Field_Peak_Firing_Rate
                    Place_Field = double(rateMap_smoothed>=maxSpatialFiringRate*Minimum_Place_Field_Firing_Rate_Fraction);
                    num_subFields=0;
                    Contiguous_Place_Fields=zeros(size(Place_Field,1),size(Place_Field,2));
                    Place_Field_Sizes=0;
                    [Y_Index,X_Index]=find(Place_Field==1);
                    for M=1:length(Y_Index)
                        if Contiguous_Place_Fields(Y_Index(M),X_Index(M))==0 && Place_Field(Y_Index(M),X_Index(M))==1
                            This_Place_Field=grayconnected(Place_Field,Y_Index(M),X_Index(M),0);
                            if sum(sum(This_Place_Field))>=Minimum_Contiguous_Place_Field_Bins
                                num_subFields=num_subFields+1;
                                Contiguous_Place_Fields(This_Place_Field)=num_subFields;
                                Place_Field_Sizes=[Place_Field_Sizes;sum(sum(This_Place_Field))];
                            end
                        end
                        clear This_Place_Field;
                    end
                    [Y,X]=find(Contiguous_Place_Fields>0);
                    Mean_InField_Firing_Rate=mean(mean(rateMap_smoothed(Y,X)));
                    clear M;
                    clear X;
                    clear Y;
                    clear X_Index;
                    clear Y_Index;
                    if length(Place_Field_Sizes)>1
                        Place_Field_Sizes=Place_Field_Sizes(2:end);
                    end
                    
                    Place_Field_Sizes(Place_Field_Sizes == 0) = [];
                else
                    Place_Field_Sizes = nan;
                    Mean_InField_Firing_Rate = nan;
                    num_subFields = 0;
                    Contiguous_Place_Fields = nan;
                end
                
                directional_place_field_sizes = [directional_place_field_sizes; Place_Field_Sizes];
                directional_mean_in_field_firing_rate(k) = Mean_InField_Firing_Rate;
                directional_peak_firing_rate(k) = maxSpatialFiringRate;
                directional_mean_firing_rate(k) = meanSpatialFiringRate;
                
                
                
                clusters(i).runs(j).directions(k).rateMap = rateMap;
                clusters(i).runs(j).directions(k).rateMap_vectorized = reshape(rateMap,1,prod(numSpatialBins));
                clusters(i).runs(j).directions(k).rateMap_smoothed = rateMap_smoothed;
                clusters(i).runs(j).directions(k).rateMap_smoothed_vectorized = reshape(rateMap_smoothed,1,prod(numSpatialBins));
                clusters(i).runs(j).directions(k).maxSpatialFiringRate = maxSpatialFiringRate;
                clusters(i).runs(j).directions(k).meanSpatialFiringRate = meanSpatialFiringRate;
                clusters(i).runs(j).directions(k).maxSpatialFiringBinLoc = maxSpatialFiringBinLoc;
                clusters(i).runs(j).directions(k).COMSpatialFiringBinLoc = COMSpatialFiringBinLoc;
                clusters(i).runs(j).directions(k).meanFiringRate = meanFiringRate;
                
                
                clusters(i).runs(j).directions(k).Place_Field_Sizes = Place_Field_Sizes;
                clusters(i).runs(j).directions(k).Mean_InField_Firing_Rate = Mean_InField_Firing_Rate;
                clusters(i).runs(j).directions(k).num_subFields = num_subFields;
                clusters(i).runs(j).directions(k).Place_Field = Contiguous_Place_Fields;
            end
            
            % Now that place fields have been found for each running
            % direction separately, you can recompute
            % Mean_InField_Firing_Rate, Place_Field_Sizes, and
            % num_subFields
            
            clusters(i).runs(j).directional_place_field_sizes = nanmean(directional_place_field_sizes);
            clusters(i).runs(j).directional_peak_firing_rate = nanmean(directional_peak_firing_rate);
            clusters(i).runs(j).directional_mean_firing_rate = nanmean(directional_mean_firing_rate);
            clusters(i).runs(j).directional_mean_in_field_firing_rate = nanmean(directional_mean_in_field_firing_rate);
            clusters(i).runs(j).num_subFields = nansum([clusters(i).runs(j).directions(1).num_subFields; clusters(i).runs(j).directions(2).num_subFields]);   
        end
    end
end


save('clusters.mat','clusters','-append')

