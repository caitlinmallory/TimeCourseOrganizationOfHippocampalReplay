function [rho,pval] = nancorr(X,Y)

if size(X,1) > size(X,2)
bad_inds = unique([find(isnan(X));find(isnan(Y))]);
else
    bad_inds = unique([find(isnan(X)) find(isnan(Y))]);
end

X(bad_inds) = [];
Y(bad_inds) = [];

[rho,pval] = corr(X,Y);