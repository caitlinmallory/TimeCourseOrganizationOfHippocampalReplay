function [maxSpatialFiringBinLoc,COMSpatialFiringBinLoc] = compute_maxSpatialFiringLoc(spatialTuningCurve,spatialDim)

if spatialDim == 2
    if nansum(spatialTuningCurve(:))~=0
        [J,I] = ind2sub(size(spatialTuningCurve),find(spatialTuningCurve(:)==max(spatialTuningCurve(:))));
        maxSpatialFiringBinLoc = [I(1) J(1)];
        COMSpatialFiringBinLoc = compute_centerOfMass(spatialTuningCurve);
    else
        maxSpatialFiringBinLoc = [NaN NaN];
        COMSpatialFiringBinLoc = [NaN NaN];
    end
elseif spatialDim == 1
    if nansum(spatialTuningCurve(:))~=0 && sum(spatialTuningCurve(:))~=inf
        maxSpatialFiringBinLoc = find(spatialTuningCurve==max(spatialTuningCurve)); maxSpatialFiringBinLoc = maxSpatialFiringBinLoc(1);
        COMSpatialFiringBinLoc = compute_centerOfMass(spatialTuningCurve); COMSpatialFiringBinLoc = COMSpatialFiringBinLoc(2);
    else
        maxSpatialFiringBinLoc = NaN;
        COMSpatialFiringBinLoc = NaN;
    end
end