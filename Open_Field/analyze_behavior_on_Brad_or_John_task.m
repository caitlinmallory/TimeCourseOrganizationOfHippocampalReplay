load Behavior_Analysis
load Position_Data
load laser_state
load Experiment_Information.mat


run_segments = [];
for i = 1:length(Experiment_Information.Segments)
if ismember(11,[Experiment_Information.Segments(i).Flags])
    run_segments = [run_segments; i];
end
end

for run = 1:max(size(Behavior_Analysis))
   
goal_well_list = Behavior_Analysis(run).goal_well_list;
entry_exit_times = Behavior_Analysis(run).entry_exit_times;

trial_times = [[Experiment_Information.Segments(run_segments(run)).Times(1); entry_exit_times(1:end-1,2)] entry_exit_times(1:end,1)];
latencies = (trial_times(:,2)-trial_times(:,1))./30000;
home_latencies = latencies(goal_well_list(:,1) == mode(goal_well_list(:,1)));
away_latencies = latencies(goal_well_list(:,1) ~= mode(goal_well_list(:,1)));

latencies(latencies<0) = nan;

ranksum(home_latencies,away_latencies)
mean(home_latencies)
mean(away_latencies)

path_length = nan(length(trial_times),1);
path_velocity = nan(length(trial_times),1);
for trial = 1:length(trial_times)
path_sub = compute_dataTemporalConcatenation(Position_Data,[trial_times(trial,1) trial_times(trial,2)]);
laser_state_sub = compute_dataTemporalConcatenation(laser_state,[trial_times(trial,1) trial_times(trial,2)]);
laser_state_sub = compute_dataInterpolation(laser_state_sub,path_sub(:,1),[]);
path_sub(laser_state_sub(:,1)==1,:)= nan;

path_length(trial) = sum(sqrt((diff(path_sub(:,2))).^2 + (diff(path_sub(:,3))).^2));
path_velocity(trial) = nanmean(path_sub(:,5));
end
path_length(isnan(latencies)) = nan;
path_velocity(isnan(latencies)) = nan;

home_path_lengths = path_length(goal_well_list(:,1) == mode(goal_well_list(:,1)));
away_path_lengths = path_length(goal_well_list(:,1) ~= mode(goal_well_list(:,1)));
ranksum(home_path_lengths,away_path_lengths)
nanmean(home_path_lengths)
nanmean(away_path_lengths)
Behavior_Analysis(run).latencies = latencies;
Behavior_Analysis(run).path_length = path_length;
Behavior_Analysis(run).path_velocity = path_velocity;
end

save('Behavior_Analysis','Behavior_Analysis')