function [gaussFilt] = setUp_gaussFilt(filterDim,filterWidth)

if min(filterDim)==1 %1d filter
    [X,Y] = meshgrid(1:filterDim(2),1:filterDim(1));
    gaussFilt = mvnpdf([X(:) Y(:)],[filterDim(2)/2 filterDim(1)/2]+0.5,filterWidth*eye(2));
    gaussFilt = reshape(gaussFilt,filterDim(2),filterDim(1));

else    %2d filter
    filterDim_1 = max(filterDim);
    [X,Y] = meshgrid(1:filterDim_1,1:filterDim_1);
    gaussFilt = mvnpdf([X(:) Y(:)],[filterDim_1/2 filterDim_1/2]+0.5,filterWidth*eye(2));
    gaussFilt = reshape(gaussFilt,filterDim_1,filterDim_1);
end

% gaussFilt = gaussFilt/sum(gaussFilt(:));

gaussFilt = gaussFilt/trapz(gaussFilt(:));

