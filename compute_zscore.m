function [Z] = compute_zscore(X)

Z = nan(size(X));
for i = 1:size(X,2)
    x = X(:,i);
    if isempty((isnan(x)))==1
        z = zscore(x);
    else
        z = (x-nanmean(x))/nanstd(x);
    end
    Z(:,i) = z;
end