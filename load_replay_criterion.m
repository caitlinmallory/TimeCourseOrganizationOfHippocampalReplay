
function [criterion_str, duration_thr, max_jump_distance_thr, coverage_thr_percent_of_track, coverage_thr, weighted_r_thr, posterior_in_map_thr] = load_replay_criterion(long_replays_only,override_coverage_thr,override_weighted_r_thr,override_posterior_diff_thr)

normalized_range = 1;

 % Criterion for "good" replays:

if long_replays_only == 1
    if normalized_range == 0
        coverage_thr_percent_of_track = -inf; %0 to 1 (set to 0 to include everything; set to 0.5 for just 'good replays')
        coverage_thr = 50; % cm
    else
         coverage_thr_percent_of_track = 0.20; %0 to 1 (set to -inf to include everything; set to 0.5 for just 'good replays')'
         coverage_thr = 0; % cm  
    end
else
    coverage_thr_percent_of_track = -inf; %0 to 1 (set to 0 to include everything; set to 0.5 for just 'good replays')
    coverage_thr = -inf; % cm
end
    
max_jump_distance_thr = 0.4;
duration_thr = 0.05;
weighted_r_thr = 0.6;
posterior_in_map_thr = 0.33; %Ambrose: 0.33. Difference in left map/right map posterior

if ~isempty(override_coverage_thr)
    coverage_thr_percent_of_track = override_coverage_thr;
end
if ~isempty(override_weighted_r_thr)
    weighted_r_thr = override_weighted_r_thr;
end
if ~isempty(override_posterior_diff_thr)
    posterior_in_map_thr = override_posterior_diff_thr;
end


criterion_str = ['maxjump' num2str(max_jump_distance_thr*100) 'weightedR' num2str(weighted_r_thr*100) 'PercCov' num2str(coverage_thr_percent_of_track*100) 'Cov' num2str(coverage_thr)];
