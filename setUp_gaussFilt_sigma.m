function w = setUp_gaussFilt_sigma(smoothing_sigma,bin_width)

% smoothing sigma is the desired standard deviation of the gaussian, (in
% seconds for time signals, could be in cm for spatial signals);
% bin_width: the width each bin in the signal (units much match smoothing
% sigma)

% N = filter_length_bins; 
% alpha = (filter_length_bins-1)/(2*sigma_bins);


smoothing_sigma_bins = smoothing_sigma/bin_width;

smoothing_window_length = smoothing_sigma_bins*5;


alpha = (smoothing_window_length-1)/(2*smoothing_sigma_bins);


w = gausswin(smoothing_window_length, alpha);
w = w./sum(w);


