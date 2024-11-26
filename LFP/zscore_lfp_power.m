clearvars -except dayFiles day directory rat windows hand_clustered_only

plot_figure = 0;
load Analysis_Information;
load Experiment_Information;
load Position_Data;

Fs = 1500; % LFP sampling rate
spikeSampRate = 30000; % spike sampling rate

t = table();
t.frequency_name = {'Delta';'Theta';'Beta';'Gamma';'Ripple'};
t.frequency_range = [1 6; 6 12; 15 20; 25 50; 150 250];
t.smoothing_sigma = [0.5; 0.3; 0.15; 0.085; 0.0125];
% Note: I am using the ripple reference tetrode for everything now, because
% I know it is in the pyramidal layer.
t.LFP_tetrode = [Experiment_Information.ripple_refs(1); Experiment_Information.ripple_refs(1); Experiment_Information.ripple_refs(1); Experiment_Information.ripple_refs(1); Experiment_Information.ripple_refs(1)];

Run_Times = Experiment_Information.Run_Times;
Sleep_Times = Experiment_Information.Sleep_Times;
unique_Session_Times = vertcat(Experiment_Information.Segments.Times);
Times_day = [min(unique_Session_Times) max(unique_Session_Times)];

%positions
times = load_timeBins_cm(Times_day,0.005*spikeSampRate,0.02*spikeSampRate);
Position_Data = load_positions_full(Run_Times,Sleep_Times,times,Position_Data,4);
Position_Data(isnan(Position_Data(:,5)),5) = 0;
zscored_lfp_power = table;

for i = 5

    lfp_tetrode = t.LFP_tetrode(i);
    lfp_file_denoised = ['LFP_Data' num2str(lfp_tetrode) '_denoised.mat'];
    lfp_file = ['LFP_Data' num2str(lfp_tetrode) '.mat'];
    
    load(lfp_file_denoised);
    load(lfp_file);

    % filter LFP for the desired frequency range and smooth
    LFP = LFP_Data;
    [LFP_phase,LFP_filt_Amp,LFP_filt] = compute_filteredLFP_cm(t.frequency_range(i,:),LFP(:,2),Fs);
    LFP = [LFP,LFP_filt,LFP_phase,LFP_filt_Amp];
    % Column 1 = time
    % Column 2 = raw lfp
    % Column 3 = x-band filtered lfp
    % Column 4 = x-band phase
    % Column 5 = x-band amplitude, smoothed
    % Column 6 = x-band zscored
    % Column 7 = x-band zscored, stationary

    smoothing_sigma = t.smoothing_sigma(i); % desired standard deviation of the Gaussian used for smoothing
    stepSize = 1/Fs;
    w = setUp_gaussFilt_sigma(smoothing_sigma,stepSize);
    smoothed_power = conv(LFP(:,5),w,'same');
    LFP(:,5) = smoothed_power;

    % Columns 8 and 9 are for theta only:
    % 8 = theta-cycle duration (trough to trough)
    % 9 = monotonically increasing flag
   
    if i == 1
        LFP(:,6:9) = nan(length(LFP),4);
        LFP_filt_theta_phase_zero_trough = mod(LFP(:,4) + pi,2*pi);
        LFP_filt_theta_phase_zero_trough(LFP_filt_theta_phase_zero_trough == 0) = 2*pi;
        % For Brad's code below to work, the trough must be 0/2*pi. I will
        % add pi to makes this the case, but save the original phases
        phase_difference=diff(LFP_filt_theta_phase_zero_trough);
        troughs=find(phase_difference<deg2rad(-345));
        troughs=troughs+1;
        monotonic_increasing=zeros(length(troughs),1);
        for n=1:(length(troughs)-1)
            if all(diff(LFP_filt_theta_phase_zero_trough(troughs(n):(troughs(n+1)-1)))>0)
                monotonic_increasing(n,1)=1;
            end
        end

        cycle_durations=zeros(length(troughs),1);
        for n=1:(length(troughs)-1)
            cycle_durations(n,1)=(LFP((troughs(n+1)-1),1)-LFP(troughs(n),1))./spikeSampRate;
        end

        for n=1:(length(troughs)-1)
            LFP(troughs(n):(troughs(n+1)-1),8:9)=ones(length(LFP(troughs(n):(troughs(n+1)-1),1)),2).*[cycle_durations(n),monotonic_increasing(n)];
        end

    end

    zscored_power = [];
    for j = 1:size(unique_Session_Times,1)

        LFP_sub = compute_dataTemporalConcatenation(LFP,unique_Session_Times(j,:));
        LFP_Data_denoised_sub = compute_dataTemporalConcatenation(LFP_Data_denoised,unique_Session_Times(j,:));

        % Do the following separately for each session:
        %load LFP ripple amp within session and zscore

        % position:
        Position_Data_sub = compute_dataTemporalConcatenation(Position_Data,unique_Session_Times(j,:));
        Position_Data_sub = compute_dataInterpolation(Position_Data_sub, LFP_sub(:,1), []);

        %mask out artifacts prior to z-scoring:
        artifact_inds = find(isnan(LFP_Data_denoised_sub(:,2)));
        moving_inds = find(Position_Data_sub(:,5)>speedThr);
        bad_inds = unique([artifact_inds; moving_inds]);

        % save the good inds (stationary) that will be used for computing mean
        % and variance.
        good_inds = setdiff(1:length(LFP_sub),bad_inds)';

        % compute mean and variance of the smoothed ripple amplitude, stationary periods only:
        mean_power = nanmean(LFP_sub(good_inds,5));
        std_power = nanstd(LFP_sub(good_inds,5));

        LFP_sub(:,6) = (LFP_sub(:,5)-mean_power)./std_power;
        LFP_sub(:,7) = LFP_sub(:,6);
        LFP_sub(bad_inds,7) = nan;

        %save a copy of the z-scored ripple power, without removing any data
        %points.

        zscored_power = [zscored_power; LFP_sub];

    end

    zscored_lfp_power.(t.frequency_name{i}) = zscored_power;

end

dt=whos('zscored_lfp_power');
if dt.bytes < 2e+09
    save('zscored_lfp_power','zscored_lfp_power')
else
    save('zscored_lfp_power','zscored_lfp_power','-v7.3')
end

% zscored_delta_power = zscored_lfp_power.Delta;
% zscored_theta_power = zscored_lfp_power.Theta;
% zscored_beta_power = zscored_lfp_power.Beta;
% zscored_gamma_power = zscored_lfp_power.Gamma;
% zscored_ripple_power = zscored_lfp_power.Ripple;

if plot_figure == 1
    figure()
    plot(zscored_lfp_power.Theta(:,1)./30000,zscored_lfp_power.Theta(:,6));
    hold on
    plot(zscored_lfp_power.Theta(:,1)./30000,zscored_lfp_power.Theta(:,7));
    hold on
    plot(zscored_lfp_power.Ripple(:,1)./30000,zscored_lfp_power.Ripple(:,7));
    legend({'Theta power','Theta power stationary','Ripple power stationary'})
end

