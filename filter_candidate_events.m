function candidateEvent = filter_candidate_events(Experiment_Information,decoder_binDecoding,sessionNum_decoder,indSub)

plot_figure = 0;
% 1) determine which map has more posterior- all calculations should be done on
% this map.
% 2) using this map- look for segments with smoothly changing posterior,
% while still allowing for small blips.
% 3) only consider the longest segment within the candidate event

binSize = 2;
spikeSampRate = 30000;

% jumps based on percentage of track:
jumpThr = 0.4;
jumpThr = jumpThr*Experiment_Information.maze_size/binSize;

% jumps based on absolute distance:
% jumpThr = 75; % cm
% jumpThr = jumpThr/binSize; % number of bins that can be jumped

sequence_durationThr = 0.03;
sequence_deltThr = 0.05;
num_mask_bins = 10;


candidate_event_posterior_full = decoder_binDecoding(sessionNum_decoder).posterior(indSub,:);
num_position_bins = size(candidate_event_posterior_full,2);
posterior_left = candidate_event_posterior_full(:,1:size(candidate_event_posterior_full,2)/2);
posterior_right = candidate_event_posterior_full(:,size(candidate_event_posterior_full,2)/2+1:size(candidate_event_posterior_full,2));
candidate_event_timeBins = decoder_binDecoding(sessionNum_decoder).timeBins(indSub,:);

good_posterior_thr = 5*(1/num_position_bins);
%good_posterior_thr = 0;

% now truncate posterior to the map of interest
% here I am masking the 10 cm on either end of the track when summing the
% total posterior in each map!
posterior_left_masked = posterior_left(:,num_mask_bins+1:end-num_mask_bins);
posterior_right_masked = posterior_right(:,num_mask_bins+1:end-num_mask_bins);
posterior_masked_sum = sum(sum(posterior_left_masked)) + sum(sum(posterior_right_masked));

if sum(sum(posterior_left_masked)) > sum(sum(posterior_right_masked))
    best_map = 1;
    candidate_event_posterior = posterior_left;
    
else
    best_map = 2;
    candidate_event_posterior = posterior_right;
end

% Here, posterior diff is being calculated before the good sequence is
% identified. not sure if this is appropriate.
posterior_diff = abs(sum(sum(posterior_left_masked))/posterior_masked_sum - sum(sum(posterior_right_masked))/posterior_masked_sum);


[candidate_event_peakPosterior,candidate_event_peakLoc] = max(candidate_event_posterior,[],2);

%decoded position
x =  candidate_event_peakLoc;
timeBins =  candidate_event_timeBins;
windowSizeDecoding = decoder_binDecoding.windowSizeDecoding;
shiftSizeDecoding = decoder_binDecoding.shiftSizeDecoding;
x = [mean(timeBins,2), x];

%decoded position jumps

jumps = abs(diff(x(:,2))); jumps = [0;jumps];

% filter bins based on jump size
%ind = find(jumps<(jumpThr*Experiment_Information.maze_size/binSize) & candidate_event_peakPosterior > good_posterior_thr);
ind = find(jumps<jumpThr & candidate_event_peakPosterior > good_posterior_thr);



if isempty(ind)

    candidateEvent.timeBins = timeBins;
    candidateEvent.boundaries = [];
    candidateEvent.x_NaN = [x(:,1) nan(size(x(:,1)))];
    candidateEvent.x = x;
    candidateEvent.best_map = nan;
    candidateEvent.boundaries = [];
    candidateEvent.posterior_diff = nan;
    return
else

    x_NaN = x;
    x_NaN(:,2) = NaN;
    x_NaN(ind,:) = x(ind,:);

    timeBins_NaNremoved = timeBins(ind,:);


    %merge overlapping time bins
    ind = [0; find(diff(timeBins_NaNremoved(:,1))>windowSizeDecoding*spikeSampRate);length(timeBins_NaNremoved(:,1))];
    timeBins_NaNremoved_merge = [timeBins_NaNremoved(ind(1:end-1)+1,1) timeBins_NaNremoved(ind(2:end),2)];
end

[boundaries_pre_merge,lengths] = compute_allSequences_NaNseparated_cm(x_NaN(:,2));

x_NaN = x; x_NaN(:,2) = NaN;
for segment = 1:size(boundaries_pre_merge)
    x_NaN(boundaries_pre_merge(segment,1):boundaries_pre_merge(segment,2),2) = x(boundaries_pre_merge(segment,1):boundaries_pre_merge(segment,2),2);
end

%merge sequences if gap is delx<delxThr and delt<deltThr
[boundaries,lengths] = compute_allSequences_NaNseparated_merge(x_NaN,boundaries_pre_merge,jumpThr,sequence_deltThr*spikeSampRate,binSize);


%remove short sequences (< sequence_durationThr)
ind_remove = find(lengths<sequence_durationThr/decoder_binDecoding(1).shiftSizeDecoding);
boundaries(ind_remove,:) = [];
lengths(ind_remove,:) = [];

% Only keep the longest sequence!

[~,ind_keep] = max(lengths);

if length(ind_keep)>0
    boundaries = boundaries(ind_keep,:);
    lengths = lengths(ind_keep);
    ind_remove = setdiff(1:length(indSub),boundaries(1):boundaries(2));
else
    ind_remove = 1:length(indSub);
end
% anything outside the selected segment should be nan'ed out
% filter bins based on jump size

x_NaN(ind_remove,2) = NaN;

% anything remaining with too low posterior should be nan'ed out
x_NaN(candidate_event_peakPosterior < good_posterior_thr,2) = nan;


boundaries = [find(~isnan(x_NaN(:,2)),1,'first') find(~isnan(x_NaN(:,2)),1,'last')];
if ~isempty(boundaries)
    boundaries = indSub(boundaries);
else
    boundaries = nan;
end

posterior_left_sub = posterior_left(~isnan(x_NaN(:,2)),:);
posterior_right_sub = posterior_right(~isnan(x_NaN(:,2)),:);

% Compute the overall difference in posterior of the sub-sequence
pcnt_posterior_in_left = sum(sum(posterior_left_sub)) / (sum(sum(posterior_left_sub))+sum(sum(posterior_right_sub)));
pcnt_posterior_in_right = sum(sum(posterior_right_sub)) / (sum(sum(posterior_left_sub))+sum(sum(posterior_right_sub)));


% figure()
% subplot(2,1,1)
% imagesc(posterior_left_sub');
% set(gca,'YDir','normal')
% colormap hot
% caxis([ 0 0.1])
% 
% subplot(2,1,2)
% imagesc(posterior_right_sub');
% set(gca,'YDir','normal')
% colormap hot
% caxis([ 0 0.1])



% If, in the selected sequence, there is actually more posterior in the
% wrong map, this sequence needs to be thrown out.
if best_map == 1 && pcnt_posterior_in_right > pcnt_posterior_in_left || ...
        best_map == 2 && pcnt_posterior_in_left > pcnt_posterior_in_right
    ind_remove = 1:length(indSub);
    boundaries = nan;
else 
    ind_remove = [];
end

x_NaN(ind_remove,2) = NaN;

candidateEvent.timeBins = timeBins;
candidateEvent.boundaries = boundaries;
candidateEvent.x_NaN = x_NaN;
candidateEvent.x = x;
candidateEvent.best_map = best_map;
candidateEvent.boundaries = boundaries;
candidateEvent.posterior_diff = posterior_diff;
candidate_event_posterior(isnan(x_NaN(:,2)),:) = NaN;

if plot_figure
figure()

subplot(2,2,1)
imagesc(posterior_left');
set(gca,'YDir','normal');
colormap hot;
caxis([0 0.1])

subplot(2,2,2)
imagesc(posterior_right');
set(gca,'YDir','normal');
colormap hot;
caxis([0 0.1])

if best_map == 1
    subplot(2,2,3)
    imagesc(candidate_event_posterior');
    set(gca,'YDir','normal');
    colormap hot;
    caxis([0 0.1])
else
    subplot(2,2,4)
    imagesc(candidate_event_posterior');
    set(gca,'YDir','normal');
    colormap hot;
    caxis([0 0.1])
end
drawnow

end 
end
