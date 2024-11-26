clearvars -except dayFiles day directory rat

day = 1;
pre_period = 3; % seconds prior to stopping time to analyze
post_period = 11; % seconds post stopping time to analyze
range_wavelet_cycles = [4 7];
num_frequency_bins = 250;
freq_range = [1 251];
linear = 0;
Fs = 1500;
plot_fig = 1;
spikeSampRate = 30000;

signal_length = (pre_period + post_period)*30000/(1/Fs*30000);
dt = 1/Fs;
time_vec = (0:signal_length-1)*dt - pre_period;

load Experiment_Information.mat
run_segments = [];
for i = 1:length(Experiment_Information.Segments)
    if ismember(11,Experiment_Information.Segments(i).Flags)
        run_segments = [run_segments; i];
    end
end

load Position_Data;
LFP_file_theta = ['LFP_Data' num2str(Experiment_Information.theta_refs(1))];

load(LFP_file_theta)
LFP_Data_theta = LFP_Data;

LFP_file_ripple = ['LFP_Data' num2str(Experiment_Information.ripple_refs)];
load(LFP_file_ripple)
LFP_Data_ripple = LFP_Data;

Position_Data = compute_dataInterpolation(Position_Data,LFP_Data(:,1),4);
load laser_state

number_segments = length(run_segments);
%%
for run_segment_ind = 1:number_segments

    if Experiment_Information.spatialDim == 1
        load Behavior_Data.mat
        start_times = Reward_Epoch_Time_Boundaries_speed_thresholded(:,1) - pre_period*30000;
        stopping_periods = Reward_Epoch_Time_Boundaries_speed_thresholded;
    elseif Experiment_Information.spatialDim == 2
        load true_drink_periods.mat
        start_times = true_drink_periods_summary(run_segment_ind).true_drink_periods(:,1) - pre_period*spikeSampRate;
        stopping_periods = true_drink_periods_summary(run_segment_ind).true_drink_periods;
    end
    stopping_period_duration = (stopping_periods(:,2)-stopping_periods(:,1))./spikeSampRate;

    all_epoch_laser_state = zeros(size(start_times,1),1);
    end_times = start_times + (pre_period + post_period)*30000; % look at first ~10 seconds after stopping
    for i = 1:length(start_times)
        laser_state_sub = compute_dataTemporalConcatenation(laser_state,[start_times(i) end_times(i)]);
        if any(laser_state_sub(:,2)==1)
            all_epoch_laser_state(i) = 1;
        end
    end

    LFP_theta_during_stopping = nan(size(start_times,1),signal_length);
    LFP_ripple_during_stopping = nan(size(start_times,1),signal_length);
    all_speed_during_stopping = nan(size(start_times,1),signal_length);

    for i = 1:size(start_times,1)
        if isnan(start_times) | isnan(end_times)
            continue
        end

        LFP_theta_sub = compute_dataTemporalConcatenation(LFP_Data_theta,[start_times(i) end_times(i)]);
        LFP_ripple_sub = compute_dataTemporalConcatenation(LFP_Data_ripple,[start_times(i) end_times(i)]);
        if length(LFP_theta_sub)~= signal_length
            continue
        end
        Position_Data_sub = compute_dataTemporalConcatenation(Position_Data,[start_times(i) end_times(i)]);
        all_speed_during_stopping(i,:) = Position_Data_sub(:,5);
        LFP_theta_during_stopping(i,:) = LFP_theta_sub(:,2);
        LFP_ripple_during_stopping(i,:) = LFP_ripple_sub(:,2);
    end


    if Experiment_Information.spatialDim == 2
        if isfield(true_drink_periods_summary,'lfp_contaminated_by_artifacts')==0 || ...
                (isfield(true_drink_periods_summary,'lfp_contaminated_by_artifacts')==1 && isempty(true_drink_periods_summary(run_segment_ind).lfp_contaminated_by_artifacts))
            inspect_stopping_periods_for_lfp_artifacts
        end
        bad_stopping_periods_lfp = find(true_drink_periods_summary(run_segment_ind).lfp_contaminated_by_artifacts ==1);
    else
        if exist("Reward_Epoch_Time_Boundaries_speed_thresholded_lfp_contaminated")==0
            inspect_stopping_periods_for_lfp_artifacts_linear_track
        end
        bad_stopping_periods_lfp = find(Reward_Epoch_Time_Boundaries_speed_thresholded_lfp_contaminated==1);
    end

    rows_to_keep = setdiff(1:size(LFP_theta_during_stopping,1),[bad_stopping_periods_lfp; find(isnan(sum(LFP_theta_during_stopping,2)))]);


    low_stopping_period_duration_thr = 6;
    high_stopping_period_duration_thr = 10;
    all_trial_inds = rows_to_keep;
    short_trial_inds = intersect(rows_to_keep,find(stopping_period_duration >= low_stopping_period_duration_thr));
    long_trial_inds = intersect(rows_to_keep,find(stopping_period_duration >= high_stopping_period_duration_thr));

    % ALL trials, regardless of trial duration:
    epoch_laser_state = all_epoch_laser_state(all_trial_inds);
    [tf, tf_off, tf_on, tf_zscored_within_trial, tf_off_zscored_within_trial, tf_on_zscored_within_trial, ...
        tf_zscored_within_session, tf_off_zscored_within_session, tf_on_zscored_within_session,times2keep,inds2keep] = plot_time_frequency_spectrograms_cm(LFP_theta_during_stopping(all_trial_inds,:)',time_vec,ones(length(all_trial_inds),1), all_epoch_laser_state(all_trial_inds), range_wavelet_cycles,num_frequency_bins,freq_range,linear, Fs, plot_fig, 'none');
    speed_during_stopping = all_speed_during_stopping(all_trial_inds,inds2keep);
    save(['tf_spectrograms_segment' num2str(run_segments(run_segment_ind))],'tf','tf_off','tf_off_zscored_within_session','tf_off_zscored_within_trial',...
        'tf_on','tf_on_zscored_within_session','tf_on_zscored_within_trial','epoch_laser_state','time_vec','times2keep','inds2keep','speed_during_stopping')
    close all

    % Trials lasting at least low_stopping_period_duration_thr seconds long:
    epoch_laser_state = all_epoch_laser_state(short_trial_inds);
    [tf, tf_off, tf_on, tf_zscored_within_trial, tf_off_zscored_within_trial, tf_on_zscored_within_trial, ...
        tf_zscored_within_session, tf_off_zscored_within_session, tf_on_zscored_within_session,times2keep,inds2keep] = plot_time_frequency_spectrograms_cm(LFP_theta_during_stopping(short_trial_inds,:)',time_vec,ones(length(short_trial_inds),1), all_epoch_laser_state(short_trial_inds), range_wavelet_cycles,num_frequency_bins,freq_range,linear, Fs, plot_fig ,'6sec');
    speed_during_stopping = all_speed_during_stopping(short_trial_inds,inds2keep);
    save(['tf_spectrograms_low_duration_thr_segment' num2str(run_segments(run_segment_ind))],'tf','tf_off','tf_off_zscored_within_session','tf_off_zscored_within_trial',...
        'tf_on','tf_on_zscored_within_session','tf_on_zscored_within_trial','epoch_laser_state','time_vec','times2keep','inds2keep','speed_during_stopping')
    close all

    % Trials lasting at least high_stopping_period_duration_thr seconds long:
    epoch_laser_state = all_epoch_laser_state(long_trial_inds);
    [tf, tf_off, tf_on, tf_zscored_within_trial, tf_off_zscored_within_trial, tf_on_zscored_within_trial, ...
        tf_zscored_within_session, tf_off_zscored_within_session, tf_on_zscored_within_session,times2keep,inds2keep] = plot_time_frequency_spectrograms_cm(LFP_theta_during_stopping(long_trial_inds,:)',time_vec,ones(length(long_trial_inds),1), all_epoch_laser_state(long_trial_inds), range_wavelet_cycles,num_frequency_bins,freq_range,linear, Fs, plot_fig,'10sec');
    speed_during_stopping = all_speed_during_stopping(long_trial_inds,inds2keep);
    save(['tf_spectrograms_high_duration_thr_segment' num2str(run_segments(run_segment_ind))],'tf','tf_off','tf_off_zscored_within_session','tf_off_zscored_within_trial',...
        'tf_on','tf_on_zscored_within_session','tf_on_zscored_within_trial','epoch_laser_state','time_vec','times2keep','inds2keep','speed_during_stopping')
    close all
end


