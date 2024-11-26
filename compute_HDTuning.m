function [rateMap, rateMap_smoothed, rateMap_smoothed_NaN, HDTuningPref, HDTuningSpec, angle_centers, posDensity_smoothed_NaN] = ...
    compute_HDTuning(spkHD,positions,speedThr,posSampRate,spatialDim)

dt = 1/posSampRate;
gaussFilt = setUp_gaussFilt([1 100],2);

numBins = 32;
if spatialDim==2
    angle_bins = linspace(-pi,pi,numBins+1)';
else
    angle_bins = [-1.5,0,1.5];
end
angle_centers = angle_bins(1:end-1)+mean(diff(angle_bins))/2;

%compute rate map
    spikeDensity = histcounts(spkHD,angle_bins); %spikeDensity = spikeDensity(1:end-1);
    posDensity = histcounts(positions(:,4),angle_bins); %posDensity = posDensity(1:end-1);
    
    rateMap = spikeDensity./(posDensity*dt);

    rateMap(isnan(rateMap(:))) = 0;
    rateMap(isinf(rateMap(:))) = 0;
    
%     if size(rateMap,2)~=size(angle_bins,2)-1
%         rateMap = rateMap';
%         posDensity = posDensity';
%     end
    rateMap = rateMap';
    posDensity = posDensity';
    
%smoothed
    rateMap_smoothed = conv([rateMap; rateMap; rateMap],gaussFilt,'same');
    rateMap_smoothed = rateMap_smoothed(end/3+1:2*end/3);

    posDensity_smoothed = conv([posDensity; posDensity; posDensity],gaussFilt,'same');
    posDensity_smoothed = posDensity_smoothed(end/3+1:2*end/3);

%replace unvisited bins with NaNs
    rateMap_smoothed_NaN = rateMap_smoothed;
    rateMap_smoothed_NaN(posDensity==0) = NaN;

    posDensity_smoothed_NaN = posDensity_smoothed;
    posDensity_smoothed_NaN(posDensity==0) = NaN;

if spatialDim==2
%     rateMap(end) = rateMap(1);
        
    [HDTuningPref,HDTuningSpec] = compute_meanVectorLength(rateMap,angle_centers);
    
elseif spatialDim==1

    HDTuningPref = NaN;
    HDTuningSpec = NaN;
    
end



        
