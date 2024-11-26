function [comX] = compute_centerOfMass_1D(X)

[numBinsX] = size(X,2);

comX = nansum([1:numBinsX].*X)/nansum(X);



