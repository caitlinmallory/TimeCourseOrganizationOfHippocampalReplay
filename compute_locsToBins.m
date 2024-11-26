function [pos] = compute_locsToBins(pos,numSpatialBins,x_edges,y_edges)

%rescales x,y
% if size(pos,2)>2
% 
%     pos(:,2) = (pos(:,2)-min(x_edges))/(max(x_edges)-min(x_edges))*numSpatialBins(2)+0.5;
%     pos(:,3) = (pos(:,3)-min(y_edges))/(max(y_edges)-min(y_edges))*numSpatialBins(1)+0.5;
% 
% elseif size(pos,2)==2
% 
%     pos(:,1) = (pos(:,1)-min(x_edges))/(max(x_edges)-min(x_edges))*numSpatialBins(2)+0.5;
%     pos(:,2) = (pos(:,2)-min(y_edges))/(max(y_edges)-min(y_edges))*numSpatialBins(1)+0.5;
% 
% else
% 
%     pos(:,1) = (pos(:,1)-min(x_edges))/(max(x_edges)-min(x_edges))*numSpatialBins(1)+0.5;
% 
% end

% Caitlin editted on 7/29/22 to remove the +0.5 Not sure what that was for
% but it was making compute_binsToLocs and compute_locsToBins asymmetrical.
if size(pos,2)>2

    pos(:,2) = (pos(:,2)-min(x_edges))/(max(x_edges)-min(x_edges))*numSpatialBins(2);
    pos(:,3) = (pos(:,3)-min(y_edges))/(max(y_edges)-min(y_edges))*numSpatialBins(1);

elseif size(pos,2)==2

    pos(:,1) = (pos(:,1)-min(x_edges))/(max(x_edges)-min(x_edges))*numSpatialBins(2);
    pos(:,2) = (pos(:,2)-min(y_edges))/(max(y_edges)-min(y_edges))*numSpatialBins(1);

else

    pos(:,1) = (pos(:,1)-min(x_edges))/(max(x_edges)-min(x_edges))*numSpatialBins(1);

end
