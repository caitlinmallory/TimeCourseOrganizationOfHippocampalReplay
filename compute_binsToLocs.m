function [pos] = compute_binsToLocs(pos,numSpatialBins,x_edges,y_edges)

%rescales x,y, followed by flipping y axis
if size(pos,2)>2
    pos(:,2) = pos(:,2)*(max(x_edges)-min(x_edges))/numSpatialBins(2) + min(x_edges);
    pos(:,3) = pos(:,3)*(max(y_edges)-min(y_edges))/numSpatialBins(1) + min(y_edges);
else
    pos(:,1) = pos(:,1)*(max(x_edges)-min(x_edges))/numSpatialBins(2) + min(x_edges);
    pos(:,2) = pos(:,2)*(max(y_edges)-min(y_edges))/numSpatialBins(1) + min(y_edges);
end

