function [timeBins] = load_timeBins_cm(Times,windowShift,windowSize)

decoding_start = min(Times(:));
decoding_end = max(Times(:));
n_decoding_bins = floor((decoding_end - decoding_start)/windowShift);

decoding_bin_start = (linspace(decoding_start, decoding_end, n_decoding_bins))';
decoding_bin_end = decoding_bin_start + windowSize;

timeBins = [decoding_bin_start decoding_bin_end];
