fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures';

colors = [.4660 0.6740 0.1880; 0.4940 0.1840 0.5560;];
transparency_pcnt = 0.5;
colors_2 = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
    [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];

new_dispThr = 0; % set to 6 to match john
new_sdThr = -inf; % default = -inf. set to 3 to match brad
new_rippleThr = -inf; % default = -inf.
new_startStopThr = 0; %set to 20 for 40 cm like brad
drink_times = [0 10];
fig_ylim = [60 100];
replay_sub = replay(replay.time_since_real_drink_onset>=drink_times(1) & replay.time_since_real_drink_onset <= drink_times(2) & ...
    replay.ripple_power >= new_rippleThr & ...
    replay.dispersion >= new_dispThr & ...
    replay.start_to_end_distance>= new_startStopThr & ...
    replay.spike_density_power>=new_sdThr,:);

% group1 = 'home';
% group2 = 'away';
% replay_sub_group1 = replay_sub(replay_sub.home_event == 1,:);
% replay_sub_group2 = replay_sub(replay_sub.home_event == 0,:);

group1 = 'off';
group2 = 'on';
replay_sub_group1 = replay_sub(replay_sub.laser_state_binary==0,:);
replay_sub_group2 = replay_sub(replay_sub.laser_state_binary==1,:);

%%
data1_fut = abs(cell2mat(replay_sub_group1.angDisplacement_futPath));
data1_past = abs(cell2mat(replay_sub_group1.angDisplacement_pastPath));
data1_fut(isnan(data1_fut))=[];
data1_past(isnan(data1_past))=[];


% inds2rmv = unique([find(isnan(data1_fut)); find(isnan(data1_past))]);
% data1_fut(inds2rmv) = [];
% data1_past(inds2rmv) = [];

data2_fut = abs(cell2mat(replay_sub_group2.angDisplacement_futPath));
data2_past = abs(cell2mat(replay_sub_group2.angDisplacement_pastPath));
data2_fut(isnan(data2_fut))=[];
data2_past(isnan(data2_past))=[];

% inds2rmv = unique([find(isnan(data2_fut)); find(isnan(data2_past))]);
% data2_fut(inds2rmv) = [];
% data2_past(inds2rmv) = [];

% group1_n = length(data1_fut)
% group2_n = length(data2_fut)

group1_n_future = [sum(~isnan(data1_fut))]
group2_n_future = [sum(~isnan(data2_fut))]
group1_n_past = [sum(~isnan(data1_past))]
group2_n_past = [sum(~isnan(data2_past))]

mean_abs_angle_future_group1 = nanmean(data1_fut);
sem_abs_angle_future_group1 = nanstd(data1_fut)./sqrt(length(data1_fut));
mean_abs_angle_past_group1 = nanmean(data1_past);
sem_abs_angle_past_group1 = nanstd(data1_past)./sqrt(sum(~isnan(data1_past)));
mean_abs_angle_future_group2 = nanmean(data2_fut);
sem_abs_angle_future_group2 = nanstd(data2_fut)./sqrt(length(data2_fut));
mean_abs_angle_past_group2 = nanmean(data2_past);
sem_abs_angle_past_group2 = nanstd(data2_past)./sqrt(sum(~isnan(data2_past)));

data = [mean_abs_angle_future_group1, mean_abs_angle_future_group2, mean_abs_angle_past_group1, mean_abs_angle_past_group2];
err = [sem_abs_angle_future_group1, sem_abs_angle_future_group2, sem_abs_angle_past_group1, sem_abs_angle_past_group2];


[p,h,z] = ranksum(data1_fut,data2_fut)
[p,h,z] = ranksum(data1_past,data2_past)

%%
figure('Position',[1986 1051 100 100])
colors = [.4660 0.6740 0.1880; 0.4940 0.1840 0.5560;];
transparency_pcnt = 0.5;
colors_2 = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
    [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];
b1 = bar(1,data(1)); hold on;
e1 = errorbar(1,data(1),err(1),'k','linestyle','none');
e1.CapSize = 4;
hold on
b2 = bar(2,data(2)); hold on;
e2 = errorbar(2,data(2),err(2),'k','linestyle','none');
e2.CapSize = 4;
b3 = bar(4,data(3)); hold on;
e3 = errorbar(4,data(3),err(3),'k','linestyle','none');
e3.CapSize = 4;
b4 = bar(5,data(4)); hold on;
e4 = errorbar(5,data(4),err(4),'k','linestyle','none');
e4.CapSize = 4;
b1.FaceColor = colors(1,:);
b1.EdgeColor = 'none';
b2.FaceColor = colors_2(1,:);
b2.EdgeColor = 'none';
b3.FaceColor = colors(2,:);
b3.EdgeColor = 'none';
b4.FaceColor = colors_2(2,:);
b4.EdgeColor = 'none';


box off
xticks([1.5 4.5])
xticklabels({})
ylim(fig_ylim)
ylabel('|Displacement|')



set(gcf, 'Color', 'white','Renderer','painters', 'PaperPositionMode', 'auto');
fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures';
saveas(gcf,fullfile(fig_path,['summary_crossing_angles_without_time_' group1 '_' group2 'rat' num2str(rats)]),'pdf')
figure()

data_table_sub = table();
data_table_sub.angles = [data1_fut; data2_fut; data1_past; data2_past];
data_table_sub.group = [zeros(length(data1_fut),1); ones(length(data2_fut),1); zeros(length(data1_past),1); ones(length(data2_past),1)];
data_table_sub.fut_past = [ones(length(data1_fut),1); ones(length(data2_fut),1); 2*ones(length(data1_past),1); 2*ones(length(data2_past),1)];
groupOrder = ["1","2"];
namedGroup = categorical(data_table_sub.group,0:1,groupOrder);
directionOrder = ["Future","Past"];
namedDirection = categorical(data_table_sub.fut_past,[1,2],directionOrder);
data_sub = data_table_sub.angles;

% Anova effects:
[p,tbl,stats] = anovan(data_sub',{namedGroup,namedDirection},'model','interaction','varnames',{'namedGroup','namedDirection'});
[results,~,~,gnames] = multcompare(stats,"Dimension",[1 2])

%%
num_shuffles=10000;
all_future = [data1_fut; data2_fut];
all_past = [data1_past; data2_past];

shuffled_future_diff = nan(num_shuffles,1);
shuffled_past_diff = nan(num_shuffles,1);

for shuffle = 1:num_shuffles
    data1_fut_shuffle_inds = randperm(height(all_future),height(data1_fut));
    data2_fut_shuffle_inds = setdiff([1:height(all_future)],data1_fut_shuffle_inds);

    data1_past_shuffle_inds = randperm(height(all_past),height(data1_past));
    data2_past_shuffle_inds = setdiff([1:height(all_past)],data1_past_shuffle_inds);

    shuffled_future_diff(shuffle) = nanmean(all_future(data1_fut_shuffle_inds)) - nanmean(all_future(data2_fut_shuffle_inds));
    shuffled_past_diff(shuffle) = nanmean(all_past(data1_past_shuffle_inds)) - nanmean(all_past(data2_past_shuffle_inds));
end
shuffled_interaction_differences = shuffled_past_diff - shuffled_future_diff;

shuffled_interaction_pval = (sum(shuffled_interaction_differences > abs(real_interaction_diff)) + sum(shuffled_interaction_differences < -1*abs(real_interaction_diff)) + 1 )/(num_shuffles+1);
shuffled_future_pval = (sum(shuffled_future_diff > abs(real_future_diff)) + sum(shuffled_future_diff < -1*abs(real_future_diff)) + 1 )/(num_shuffles+1);
shuffled_past_pval = (sum(shuffled_past_diff > abs(real_past_diff)) + sum(shuffled_past_diff < -1*abs(real_past_diff)) + 1 )/(num_shuffles+1);




