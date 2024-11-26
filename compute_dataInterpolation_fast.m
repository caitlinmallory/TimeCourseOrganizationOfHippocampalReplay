function data_interp = compute_dataInterpolation_fast(data,times)

% data should be a matrix where the first column is sample times
% times is the a vector containing the times you are interested in
% (i.e., spike times).
% this will simply pull out the nearest point in data and return the
% row of data at that point.


nearest_ind = nearestpoint(times,data(:,1));
data_interp = data(nearest_ind,:);
