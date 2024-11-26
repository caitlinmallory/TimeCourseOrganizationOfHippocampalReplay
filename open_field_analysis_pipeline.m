load_spikes_behavior_cm
% will take spike times for each cell, and find the animnal's position at
% each time. Updates Clusters.mat with this information


load_clusters_rate_maps_cm
% Generate rate maps for each session and attach them to the clusters
% struct;

load_excitatory_designation
% adds on whether a cell is excitatory or inhibitory to Clusters.mat

load_spikeDensity_cm
% computes spike density for the whole session

load_binDecoding_cm(0.08,0.005);
% decodes the entire session in 80 ms bins shifting by 5 ms
% stores decoder_binDecoding, which has the timestamps for each decoded
% frame, and some information about the posterior (i.e., center of mass and
% spread). The output is called binDecoding_08.mat. If you don't see this
% in any of the folders, you can run the above again to generate it.

load_replayEvents_cm
% John's method of replay detection- finds segments of time with smoothly
% changing posterior. Saves these segments in a structure called
% 'decoder_replay'. The replay posterior isn't saved here because it would
% be too large, but using the times or inds of the replay, you can easily
% generate it again.

plot_replayEvents_mosaic_cm
% plots all replay events that criterion

plot_rate_maps
% plots all rate maps for a session

plot_nonlocalReplay
% plots all replay, taking out any segments where the posterior is close to
% the animal's current location

plot_nonlocalPosterior
% this is basically a less restrictive version of the above. Instead of
% looking for smooth replay sequences, it plots the posterior (decoded
% position of the rat) for all time frames that meet certain criterion.

