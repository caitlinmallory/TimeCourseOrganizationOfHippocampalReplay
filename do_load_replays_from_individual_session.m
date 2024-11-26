


min_time_into_stopping_period = 0;
max_time_into_stopping_period = inf;
coverage_thr = 0; % for replays
weighted_r_thr = 0.6; % for replays
posterior_diff_thr = 0.33;

% Find your ripples
replay_selection = 'ripple_filtered';
sde_thr = -inf;
ripple_thr = 3;
use_duration_og = 1; % set to 1 if looking at ripples or sde's. set to 0 if looking at replays.
t_ripples = load_replays_from_individual_session(replay_selection,use_duration_og,-inf,-inf,posterior_diff_thr,sde_thr,ripple_thr,min_time_into_stopping_period,max_time_into_stopping_period);

% Find your replays
replay_selection = 'spike_filtered';
sde_thr = 3;
ripple_thr = -inf;
use_duration_og = 0; % set to 1 if looking at ripples or sde's. set to 0 if looking at replays.
t = load_replays_from_individual_session(replay_selection,use_duration_og,coverage_thr,weighted_r_thr,posterior_diff_thr,sde_thr,ripple_thr,min_time_into_stopping_period,max_time_into_stopping_period);
pull_ripple_sd_properties_from_decoder_events;
t.control_experimental_flag = zeros(height(t),1);
t.rat_label = ones(height(t),1);
t.day_label = ones(height(t),1);
t.weighted_r = abs(t.weighted_r);
t_replay = t(t.replay==1,:);

