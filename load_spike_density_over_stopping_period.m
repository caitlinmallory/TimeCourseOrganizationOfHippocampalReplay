plot_fig = 0;

load zscored_sd_pyr
load Experiment_Information
spikeSampRate = 30000;

if Experiment_Information.spatialDim == 2
load true_drink_periods.mat
else
load Behavior_Data.mat
end

run_segments = [];
for i = 1:length(Experiment_Information.Segments)
    if ismember(11,Experiment_Information.Segments(i).Flags)
        run_segments = [run_segments; i];
    end
end

load Position_Data
number_segments = length(run_segments);
%%
for run_segment_ind = 1:number_segments

if Experiment_Information.spatialDim == 2
    stopping_period_times = true_drink_periods_summary(run_segment_ind).true_drink_periods;
else
    stopping_period_times = Reward_Epoch_Time_Boundaries_speed_thresholded;
end


pre_period = 10; %s;
post_period = 10; %s;

desired_length = (pre_period+post_period)./(0.1) + 1;

spiking_run_to_stop = nan(length(stopping_period_times),desired_length);
running_speed_run_to_stop = nan(length(stopping_period_times),desired_length);
for i = 1:length(stopping_period_times)
     trial_duration = (stopping_period_times(i,2)-stopping_period_times(i,1))./spikeSampRate;
    spiking_sub = compute_dataTemporalConcatenation(zscored_sd,[stopping_period_times(i,1)-pre_period*spikeSampRate stopping_period_times(i,1)+post_period*spikeSampRate]);
    running_speed_sub = compute_dataTemporalConcatenation(Position_Data,[stopping_period_times(i,1)-pre_period*spikeSampRate stopping_period_times(i,1)+post_period*spikeSampRate]);
    
    desired_times = (stopping_period_times(i,1)-pre_period*spikeSampRate:0.1*spikeSampRate:stopping_period_times(i,1)+post_period*spikeSampRate)';
    spiking_sub_interp = compute_dataInterpolation(spiking_sub,desired_times,[]);
    spiking_sub_interp(spiking_sub_interp(:,1)>stopping_period_times(i,2),2:3) = nan;
    spiking_run_to_stop(i,:) = spiking_sub_interp(:,2);
    running_speed_sub_interp = compute_dataInterpolation(Position_Data,desired_times,[]);
    running_speed_sub_interp(running_speed_sub_interp(:,1)>stopping_period_times(i,2),2:end) = nan;
    running_speed_run_to_stop(i,:) = running_speed_sub_interp(:,5);
end



if plot_fig == 1
figure()
subplot(3,3,1)
imagesc(spiking_run_to_stop)
xticks(linspace(1,desired_length,pre_period+post_period+1))
xticklabels(num2cell(-1*pre_period:1:post_period))
hold on
h = xline(pre_period/0.1+1); h.Color = [1 1 1]; h.LineWidth = 2;
caxis([0 2])
xlabel('Time since arrival (s)')
ylabel('Trial')

subplot(3,3,4)
plot(linspace(-1*pre_period,post_period,desired_length),nanmean(spiking_run_to_stop));
hold on
xline(0)
subplot(3,3,7)
plot(linspace(-1*pre_period,post_period,desired_length),nanmean(running_speed_run_to_stop));
hold on
xline(0)


subplot(3,3,2)
imagesc(spiking_run_to_stop(true_drink_periods_summary(run_segment_ind).home_period==1,:))
xticks(linspace(1,desired_length,pre_period+post_period+1))
xticklabels(num2cell(-1*pre_period:1:post_period))
hold on
h = xline(pre_period/0.1+1); h.Color = [1 1 1]; h.LineWidth = 2;
caxis([0 2])
xlabel('Time since arrival (s)')
ylabel('Home trials')

subplot(3,3,5)
plot(linspace(-1*pre_period,post_period,desired_length),nanmean(spiking_run_to_stop(true_drink_periods_summary(run_segment_ind).home_period==1,:)));
hold on
xline(0)
subplot(3,3,8)
plot(linspace(-1*pre_period,post_period,desired_length),nanmean(running_speed_run_to_stop(true_drink_periods_summary(run_segment_ind).home_period==1,:)));
hold on
xline(0)


subplot(3,3,3)
imagesc(spiking_run_to_stop(true_drink_periods_summary(run_segment_ind).away_period==1,:))
xticks(linspace(1,desired_length,pre_period+post_period+1))
xticklabels(num2cell(-1*pre_period:1:post_period))
hold on
h = xline(pre_period/0.1+1); h.Color = [1 1 1]; h.LineWidth = 2;
caxis([0 2])
h = xline(5*1500); h.Color = [1 1 1]; h.LineWidth = 2;
caxis([0 2])
xlabel('Time since arrival (s)')
ylabel('Away trials')

subplot(3,3,6)
plot(linspace(-1*pre_period,post_period,desired_length),nanmean(spiking_run_to_stop(true_drink_periods_summary(run_segment_ind).away_period==1,:)));
hold on
xline(0)
subplot(3,3,9)
plot(linspace(-1*pre_period,post_period,desired_length),nanmean(running_speed_run_to_stop(true_drink_periods_summary(run_segment_ind).away_period==1,:)));
hold on
xline(0)

end


save(['spiking_run_to_stop_segment_' num2str(run_segments(run_segment_ind))],'spiking_run_to_stop','running_speed_run_to_stop')

end

