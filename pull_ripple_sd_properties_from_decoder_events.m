
%% Pull our ripple properties:
% Default is 125 (12.5 ms smoothing)
t.replay_ripple_power = t.ripple_power_sd_metrics.zscored_ripple_power.peak_power_over_fraction_of_event(:,4);
t.replay_ripple_power_first_third = t.ripple_power_sd_metrics.zscored_ripple_power.peak_power_over_fraction_of_event(:,1);
t.replay_ripple_power_second_third = t.ripple_power_sd_metrics.zscored_ripple_power.peak_power_over_fraction_of_event(:,2);
t.replay_ripple_power_third_third = t.ripple_power_sd_metrics.zscored_ripple_power.peak_power_over_fraction_of_event(:,3);
t.mean_replay_ripple_power = t.ripple_power_sd_metrics.zscored_ripple_power.mean_power_over_fraction_of_event(:,4);
t.mean_replay_ripple_power_first_third = t.ripple_power_sd_metrics.zscored_ripple_power.mean_power_over_fraction_of_event(:,1);
t.mean_replay_ripple_power_second_third = t.ripple_power_sd_metrics.zscored_ripple_power.mean_power_over_fraction_of_event(:,2);
t.mean_replay_ripple_power_third_third = t.ripple_power_sd_metrics.zscored_ripple_power.mean_power_over_fraction_of_event(:,3);

% 1.5 ms smoothing:
t.replay_ripple_power_15 = t.ripple_power_sd_metrics.zscored_ripple_power_15.peak_power_over_fraction_of_event(:,4);
t.replay_ripple_power_15_first_third = t.ripple_power_sd_metrics.zscored_ripple_power_15.peak_power_over_fraction_of_event(:,1);
t.replay_ripple_power_15_second_third = t.ripple_power_sd_metrics.zscored_ripple_power_15.peak_power_over_fraction_of_event(:,2);
t.replay_ripple_power_15_third_third = t.ripple_power_sd_metrics.zscored_ripple_power_15.peak_power_over_fraction_of_event(:,3);
t.mean_replay_ripple_power_15 = t.ripple_power_sd_metrics.zscored_ripple_power_15.mean_power_over_fraction_of_event(:,4);
t.mean_replay_ripple_power_15_first_third = t.ripple_power_sd_metrics.zscored_ripple_power_15.mean_power_over_fraction_of_event(:,1);
t.mean_replay_ripple_power_15_second_third = t.ripple_power_sd_metrics.zscored_ripple_power_15.mean_power_over_fraction_of_event(:,2);
t.mean_replay_ripple_power_15_third_third = t.ripple_power_sd_metrics.zscored_ripple_power_15.mean_power_over_fraction_of_event(:,3);
%% Pull out spike density properties:
% Default is 125 (12.5 ms smoothing)
t.replay_spikeDensity_power = t.ripple_power_sd_metrics.zscored_sd_pyr_125.peak_power_over_fraction_of_event(:,4);
t.replay_spikeDensity_power_first_third = t.ripple_power_sd_metrics.zscored_sd_pyr_125.peak_power_over_fraction_of_event(:,1);
t.replay_spikeDensity_power_second_third = t.ripple_power_sd_metrics.zscored_sd_pyr_125.peak_power_over_fraction_of_event(:,2);
t.replay_spikeDensity_power_third_third = t.ripple_power_sd_metrics.zscored_sd_pyr_125.peak_power_over_fraction_of_event(:,3);
t.mean_replay_spikeDensity_power = t.ripple_power_sd_metrics.zscored_sd_pyr_125.mean_power_over_fraction_of_event(:,4);
t.mean_replay_spikeDensity_power_first_third = t.ripple_power_sd_metrics.zscored_sd_pyr_125.mean_power_over_fraction_of_event(:,1);
t.mean_replay_spikeDensity_power_second_third = t.ripple_power_sd_metrics.zscored_sd_pyr_125.mean_power_over_fraction_of_event(:,2);
t.mean_replay_spikeDensity_power_third_third = t.ripple_power_sd_metrics.zscored_sd_pyr_125.mean_power_over_fraction_of_event(:,3);

% 1.5 ms smoothing
t.replay_spikeDensity_power = t.ripple_power_sd_metrics.zscored_sd_pyr_15.peak_power_over_fraction_of_event(:,4);
t.replay_spikeDensity_power_first_third = t.ripple_power_sd_metrics.zscored_sd_pyr_15.peak_power_over_fraction_of_event(:,1);
t.replay_spikeDensity_power_second_third = t.ripple_power_sd_metrics.zscored_sd_pyr_15.peak_power_over_fraction_of_event(:,2);
t.replay_spikeDensity_power_third_third = t.ripple_power_sd_metrics.zscored_sd_pyr_15.peak_power_over_fraction_of_event(:,3);
t.mean_replay_spikeDensity_power = t.ripple_power_sd_metrics.zscored_sd_pyr_15.mean_power_over_fraction_of_event(:,4);
t.mean_replay_spikeDensity_power_first_third = t.ripple_power_sd_metrics.zscored_sd_pyr_15.mean_power_over_fraction_of_event(:,1);
t.mean_replay_spikeDensity_power_second_third = t.ripple_power_sd_metrics.zscored_sd_pyr_15.mean_power_over_fraction_of_event(:,2);
t.mean_replay_spikeDensity_power_third_third = t.ripple_power_sd_metrics.zscored_sd_pyr_15.mean_power_over_fraction_of_event(:,3);

%% Pull out inhibitory spike density properties:
% Just used 1.5 ms smoothing
t.replay_spikeDensity_power_inhib_15 = t.ripple_power_sd_metrics.zscored_sd_inhib_15.peak_power_over_fraction_of_event(:,4);
t.replay_spikeDensity_power_inhib_15_first_third = t.ripple_power_sd_metrics.zscored_sd_inhib_15.peak_power_over_fraction_of_event(:,1);
t.replay_spikeDensity_power_inhib_15_second_third = t.ripple_power_sd_metrics.zscored_sd_inhib_15.peak_power_over_fraction_of_event(:,2);
t.replay_spikeDensity_power_inhib_15_third_third = t.ripple_power_sd_metrics.zscored_sd_inhib_15.peak_power_over_fraction_of_event(:,3);
t.mean_replay_spikeDensity_power_inhib_15 = t.ripple_power_sd_metrics.zscored_sd_inhib_15.mean_power_over_fraction_of_event(:,4);
t.mean_spikeDensity_power_inhib_15_first_third = t.ripple_power_sd_metrics.zscored_sd_inhib_15.mean_power_over_fraction_of_event(:,1);
t.mean_spikeDensity_power_inhib_15_second_third = t.ripple_power_sd_metrics.zscored_sd_inhib_15.mean_power_over_fraction_of_event(:,2);
t.mean_spikeDensity_power_inhib_15_third_third = t.ripple_power_sd_metrics.zscored_sd_inhib_15.mean_power_over_fraction_of_event(:,3);
