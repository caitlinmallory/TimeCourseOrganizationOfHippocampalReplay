function [c] = compute_centerOfMass(X)

[numBinsY,numBinsX] = size(X);

comX = dot(nansum(X,1),1:numBinsX)/nansum(nansum(X,1))+0.5/numBinsX;
comY = dot(nansum(X,2),1:numBinsY)/nansum(nansum(X,2))+0.5/numBinsY;

c = [comX,comY];