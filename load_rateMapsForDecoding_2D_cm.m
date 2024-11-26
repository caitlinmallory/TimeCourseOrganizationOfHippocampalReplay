function [ind_cluster,rateMap_smoothed,rateMap,rateMap_smoothed_NaN,rateMap_smoothed_all,rateMap_all,rateMap_smoothed_NaN_all,meanFiringRate,spatialInfo,noiseOverlap,isolation,peakSnr,rateMap_coarse_smoothed,rateMap_coarse,rateMap_coarse_smoothed_NaN,rateMap_coarse_smoothed_all,stability,num_subFields,meanSpatialFiringRate,rateOverlap,rateOverlap_meanSpatialFiringRate,cellBarrierSim] = ...
    load_rateMapsForDecoding_2D_cm(clusters,sessionNum_decoder,Run_Times,numSpatialBins,meanFiringRateThr,spatialInfoThr,noiseOverlapThr,isolationThr,peakSnrThr,numSpatialBins_coarse)


isolationThr = 0;
noiseOverlapThr = inf;
peakSnrThr = 0;

%select which clusters to use for decoding

clusters_runs = [clusters.runs]; clusters_runs = clusters_runs(sessionNum_decoder:length(Run_Times):end);
spatialInfo = load_fieldVec(clusters_runs,'spatialInfo',1);
meanFiringRate = load_fieldVec(clusters_runs,'meanFiringRate',1);
meanSpatialFiringRate = load_fieldVec(clusters_runs,'meanSpatialFiringRate',1);
num_subFields = load_fieldVec(clusters_runs,'num_subFields',1);

%     stability = load_fieldVec(clusters,'stability',size(compute_pairs(size(Run_Times,1)),1));
%     rateOverlap = load_fieldVec(clusters,'rateOverlap',size(compute_pairs(size(Run_Times,1)),1));
%     rateOverlap_meanSpatialFiringRate = load_fieldVec(clusters,'rateOverlap_meanSpatialRate',size(compute_pairs(size(Run_Times,1)),1));
%     cellBarrierSim = load_fieldVec(clusters,'cellBarrierSim',size(compute_pairs(size(Run_Times,1)),1));

if ~isfield(clusters,'tag') %MClust
    ind_cluster = find(meanFiringRate>meanFiringRateThr & spatialInfo>spatialInfoThr);
    noiseOverlap = NaN;
    isolation=NaN;
    peakSnr = NaN;

else %mountainsort
    tag = load_fieldVec(clusters,'tag_abbrev',1);
    noiseOverlap = load_fieldVec(clusters,'noise_overlap',1);
    isolation = load_fieldVec(clusters,'isolation',1);
    peakSnr = load_fieldVec(clusters,'peak_snr',1);

    if isfield(clusters,'tag')
        ind_cluster = find(noiseOverlap<noiseOverlapThr & isolation>isolationThr & peakSnr>peakSnrThr & meanFiringRate>meanFiringRateThr & spatialInfo>spatialInfoThr & (strcmp(string(tag),'A')|strcmp(string(tag),'B')));
    else
        ind_cluster = find(noiseOverlap<noiseOverlapThr & isolation>isolationThr & peakSnr>peakSnrThr & meanFiringRate>meanFiringRateThr & spatialInfo>spatialInfoThr);
    end
end

%load rate maps
rateMap_all = load_fieldVec(clusters_runs,'rateMap_vectorized',prod(numSpatialBins));
rateMap = rateMap_all(ind_cluster,:);

rateMap_smoothed_all = load_fieldVec(clusters_runs,'rateMap_smoothed_vectorized',prod(numSpatialBins));
rateMap_smoothed = rateMap_smoothed_all(ind_cluster,:);

rateMap_smoothed_NaN_all = load_fieldVec(clusters_runs,'rateMap_smoothed_NaN_vectorized',prod(numSpatialBins));
rateMap_smoothed_NaN = rateMap_smoothed_NaN_all(ind_cluster,:);

rateMap_coarse_all = load_fieldVec(clusters_runs,'rateMap_coarse_vectorized',prod(numSpatialBins_coarse));
rateMap_coarse = rateMap_coarse_all(ind_cluster,:);

rateMap_coarse_smoothed_all = load_fieldVec(clusters_runs,'rateMap_coarse_smoothed_vectorized',prod(numSpatialBins_coarse));
rateMap_coarse_smoothed = rateMap_coarse_smoothed_all(ind_cluster,:);

rateMap_coarse_smoothed_NaN_all = load_fieldVec(clusters_runs,'rateMap_coarse_smoothed_NaN_vectorized',prod(numSpatialBins_coarse));
rateMap_coarse_smoothed_NaN = rateMap_coarse_smoothed_NaN_all(ind_cluster,:);

