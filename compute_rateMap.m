function [rateMap, rateMap_smoothed, rateMap_smoothed_NaN, posDensity,spikeDensity] = ...
    compute_rateMap(spkPos_filt,Position_Data,numSpatialBins,numSpatialBins_smoothing,posSampRate,spatialDim,edges_spatialBins,speedThr,binSize)


min_place_field_activation = 1e-5;
min_time_spent_in_bin = 0.05;

dt = 1/posSampRate; % in seconds


if spatialDim == 2
    
    %binning for calculating spatial tuning curves
    x_edges = edges_spatialBins{1};
    y_edges = edges_spatialBins{2};
                                        
    %smoothed with gaussian filt
    gaussFilt = setUp_gaussFilt(numSpatialBins,numSpatialBins_smoothing/binSize);


    %histogram of animal position data
    if size(Position_Data,2)>3
        posDensity = hist3([Position_Data(Position_Data(:,5)>speedThr,3) Position_Data(Position_Data(:,5)>speedThr,2)],'Edges',{y_edges x_edges}); posDensity = posDensity(1:end-1,1:end-1);
    else
        posDensity = hist3([Position_Data(:,3) Position_Data(:,2)],'Edges',{y_edges x_edges}); posDensity = posDensity(1:end-1,1:end-1);
    end
    
    if isempty(spkPos_filt)==0
                
        %histogram of spike position data
        spikeDensity = hist3([spkPos_filt(:,2) spkPos_filt(:,1)],'Edges',{y_edges x_edges}); spikeDensity = spikeDensity(1:end-1,1:end-1);
        
        % Create an occupancy mask: positions for which the animal spent very
        % little time will essentially get blanked out of decoding
        occupancy_cutoff = median(posDensity(posDensity>0)) * min_time_spent_in_bin;
        occupancy_mask = posDensity < occupancy_cutoff; 
    
        %calculate spatial tuning curv
        rateMap = spikeDensity./(posDensity*dt);
        
        %if spike density = 0, set rateMap = 0
        rateMap(isnan(rateMap(:))) = 0;
        rateMap(isinf(rateMap(:))) = 0;
        
        %smoothed
        % Caitlin changed the smoothing on 1/18/22 to match Brad's
        rateMap_smoothed = conv2(rateMap,gaussFilt,'same');
        
        % smooth place fields with a Gaussian filter
        % rateMap_smoothed2 = imgaussfilt(rateMap,1);

        %replace unvisited bins with NaNs
        rateMap_smoothed_NaN = rateMap_smoothed;
        rateMap_smoothed_NaN(posDensity==0) = NaN;
        
        %Caitlin added the following:
        %rateMap(occupancy_mask) = min_place_field_activation;
        rateMap(rateMap == 0) = min_place_field_activation;
        %rateMap_smoothed(occupancy_mask) = min_place_field_activation;
        rateMap_smoothed(rateMap_smoothed == 0) = min_place_field_activation;
        
        %interpolate at missing values, convolve, then set back to nan
%             ind_NaN = find(isnan(rateMap_NaN));
%             ind_good = find(~isnan(rateMap_NaN));
%             [X,Y] = meshgrid(1:numSpatialBins(2),1:numSpatialBins(1)); 
%             F = scatteredInterpolant(X(ind_good),Y(ind_good),rateMap_NaN(ind_good),'nearest');
%             F_NaN = F(X(ind_NaN),Y(ind_NaN));
%             
%             rateMap_smoothed = rateMap_NaN;
%             rateMap_smoothed(ind_NaN) = F_NaN; 
%             
%             rateMap_smoothed = conv2(rateMap_smoothed,gaussFilt,'same');
%             
%             rateMap_NaN_smoothed = rateMap_smoothed;
%             rateMap_NaN_smoothed(ind_NaN) = nan;
                        
        
        %plot spikes on top of spatial tuning curve
%             subplot(221)
%                 imagesc(rateMap), set(gca,'ydir','normal')
%                 spkPos_filt_bins = compute_locsToBins(spkPos_filt,numSpatialBins,x_edges,y_edges);
%                 hold on, plot(spkPos_filt_bins(:,1),spkPos_filt_bins(:,2),'r.'), hold off, axis square
%                 
%             subplot(222)
%                 imagesc(rateMap_smoothed), set(gca,'ydir','normal')
%                 hold on, plot(spkPos_filt_bins(:,1),spkPos_filt_bins(:,2),'r.'), hold off, axis square
%                 
%             subplot(223)
%                 imagesc(log(posDensity)), set(gca,'ydir','normal'), axis square
%                 positions_bins = compute_locsToBins(positions,numSpatialBins,x_edges,y_edges);
%                 hold on, plot(positions_bins(:,2),positions_bins(:,3),'r.'), hold off, axis square
%                 
%             subplot(224)
%                 imagesc(spikeDensity), set(gca,'ydir','normal')
%                 hold on, plot(spkPos_filt_bins(:,1),spkPos_filt_bins(:,2),'r.'), hold off, axis square
% 
%             keyboard
            
    else
        rateMap = min_place_field_activation*ones(size(posDensity));
        rateMap_smoothed = min_place_field_activation*ones(size(posDensity));
        rateMap_smoothed_NaN = min_place_field_activation*ones(size(posDensity));
        spikeDensity = min_place_field_activation*ones(size(posDensity));
    end
    
elseif spatialDim == 1
    
    x_edges = edges_spatialBins{1};
%     [x_edges, ~] = setUp_spatialBinning(positions,numSpatialBins);
    
    %histogram of animal position data
    posDensity = histc(Position_Data(Position_Data(:,5)>speedThr,2),x_edges); posDensity = posDensity(1:end-1);
    
    % convert 
    occupancy_cutoff = median(posDensity(posDensity>0)) * min_time_spent_in_bin;
    occupancy_mask = posDensity < occupancy_cutoff;
 

    %gauss filt
%     gaussFilt = setUp_gaussFilt([1 length(x_edges)-1],numSpatialBins_smoothing);
    
    % Caitlin changed the smoothing to be Brad's version:
    % smooth place fields with a Gaussian filter
   
    N = 9; sigma = 1;
    alpha = (N - 1)/(2*sigma);
    w = gausswin(N,alpha)./sum(gausswin(N,alpha));


    %histogram of spike position data
    if isempty(spkPos_filt)==0
        
        spikeDensity = histc(spkPos_filt,x_edges); spikeDensity = spikeDensity(1:end-1);
        [a b] = size(spikeDensity);
        if a<b
            spikeDensity = spikeDensity';
        end

        %calculate spatial tuning curve
        rateMap = spikeDensity./posDensity/dt;
        rateMap(isnan(rateMap)) = 0;
        rateMap(isinf(rateMap)) = 0;

        %smoothed
        %rateMap_smoothed = conv(rateMap,gaussFilt,'same');
        % Caitlin changed to be Brad's version of smoothing:
     
        rateMap_smoothed = filtfilt(w,1,rateMap);
        rateMap_smoothed(rateMap_smoothed < 0) = 0;
        
        %replace unvisited bins with NaNs
        rateMap_smoothed_NaN = rateMap_smoothed;
        rateMap_smoothed_NaN(posDensity==0) = NaN;
        
        %Caitlin added the following:
        rateMap(occupancy_mask,:) = min_place_field_activation;
        rateMap_smoothed(occupancy_mask,:) = min_place_field_activation;
        rateMap(rateMap==0,:) = min_place_field_activation;
        rateMap_smoothed(rateMap_smoothed==0,:) = min_place_field_activation;
        
        %plot spikes on top of spatial tuning curve
%             subplot(311)
%                 plot(posDensity),
%             subplot(312)
%                 pos_bins = compute_locsToBins(spkPos_filt,numSpatialBins,x_edges);
%                 plot(pos_bins(:,1),0,'ko')
%                 hold on, plot(spikeDensity), hold off
%             subplot(313)
%                 plot(pos_bins(:,1),0,'ko')
%                 hold on, plot(rateMap_smoothed), hold off
%             keyboard

    else
        rateMap = min_place_field_activation*ones(size(posDensity));
        rateMap_smoothed = min_place_field_activation*ones(size(posDensity));
        rateMap_smoothed_NaN = min_place_field_activation*ones(size(posDensity));
        spikeDensity = min_place_field_activation*ones(size(posDensity));
    end

    
end