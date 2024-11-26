function [d] = compute_sequenceDistance_cumsum(X)

d = sqrt(diff(X(:,1)).^2 + diff(X(:,2)).^2);
d = cumsum(d);

