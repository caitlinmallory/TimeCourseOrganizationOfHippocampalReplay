function [spread,imageMoment] = compute_imageSpread(X,exponent)

[numBinsY,numBinsX] = size(X);
if numBinsY > 1
comX = dot(nansum(X,1),1:numBinsX)/nansum(nansum(X,1))+0.5/numBinsX;
comY = dot(nansum(X,2),1:numBinsY)/nansum(nansum(X,2))+0.5/numBinsY;

imageMoment = nan(numBinsX,numBinsY);
[XX,YY] = meshgrid(1:numBinsX,1:numBinsY);
for i = 1:numBinsY
    for j = 1:numBinsX
        imageMoment(i,j) = abs(XX(i,j)-comX)^exponent*abs(YY(i,j)-comY)^exponent*X(i,j);
    end
end
imageMoment = nansum(imageMoment(:))/nansum(X(:));
spread = sqrt(imageMoment);

else

% Caitlin added the following: compute spread for 1D case.
comX = compute_centerOfMass_1D(X);

imageMoment = nan(numBinsX,1);
[XX] = meshgrid(1:numBinsX,1);
for j = 1:numBinsX
        imageMoment(j) = abs(XX(j)-comX)^exponent*X(j);
end
imageMoment = nansum(imageMoment(:))/nansum(X(:));
spread = sqrt(imageMoment);


end

