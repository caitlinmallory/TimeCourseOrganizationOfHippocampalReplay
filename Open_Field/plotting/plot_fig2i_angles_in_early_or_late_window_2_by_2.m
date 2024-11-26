time_windows = [0 3; 3 10];

future_early = replay.meanAngDisplacement_futPath(replay.time_since_real_drink_onset>=time_windows(1,1) & ...
    replay.time_since_real_drink_onset<time_windows(1,2));
past_early = replay.meanAngDisplacement_pastPath(replay.time_since_real_drink_onset>=time_windows(1,1) & ...
    replay.time_since_real_drink_onset<time_windows(1,2));

future_late = replay.meanAngDisplacement_futPath(replay.time_since_real_drink_onset>=time_windows(2,1) & ...
    replay.time_since_real_drink_onset<time_windows(2,2));
past_late = replay.meanAngDisplacement_pastPath(replay.time_since_real_drink_onset>=time_windows(2,1) & ...
    replay.time_since_real_drink_onset<time_windows(2,2));

n_early = height(future_early);
n_late = height(future_late);

[p,h,z] = signrank(future_early,past_early);
[p,h,z] = signrank(future_late,past_late);
[p,h,z] = ranksum(future_early,future_late);
[p,h,z] = ranksum(past_early,past_late);

sum(~isnan(future_early)) 
sum(~isnan(past_early))
sum(~isnan(future_late)) 
sum(~isnan(past_late))

data = [nanmean(future_early),nanmean(past_early),nanmean(future_late) nanmean(past_late)];
err = [nanstd(future_early)/sqrt(nansum(~isnan(future_early))),...
    nanstd(past_early)/sqrt(nansum(~isnan(past_early))),...
    nanstd(future_late)/sqrt(nansum(~isnan(future_late))),...
     nanstd(past_late)/sqrt(nansum(~isnan(past_late)))]  ; 

colors = [.4660 0.6740 0.1880; 0.4940 0.1840 0.5560;];
transparency_pcnt = 0.5;
colors_2 = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
    [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];

figure('Position',[1986 1051 100 100])
b1 = bar(1,data(1)); hold on;
e1 = errorbar(1,data(1),err(1),'k','linestyle','none')
e1.CapSize = 4;
hold on
b2 = bar(2,data(2)); hold on;
e2 = errorbar(2,data(2),err(2),'k','linestyle','none')
e2.CapSize = 4;
b3 = bar(4,data(3)); hold on;
e3 = errorbar(4,data(3),err(3),'k','linestyle','none')
e3.CapSize = 4;
b4 = bar(5,data(4)); hold on;
e4 = errorbar(5,data(4),err(4),'k','linestyle','none')
e4.CapSize = 4;
b1.FaceColor = colors(1,:);
b1.EdgeColor = 'none';
b2.FaceColor = colors(2,:);
b2.EdgeColor = 'none';
b3.FaceColor = colors(1,:);
b3.EdgeColor = 'none';
b4.FaceColor = colors(2,:);
b4.EdgeColor = 'none';
ylim([50 100])


set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
saveas(gcf,fullfile(fig_path,'ang_displacements_forward_v_reverse_window'),'pdf')


early_diff = mean(past_early)-mean(future_early);
late_diff = mean(past_late)-mean(future_late);
real_interaction = early_diff-late_diff;

all_early_angles = [past_early; future_early];
all_late_angles = [past_late; future_late];
%shuffle the angles?

early_diff_shuffled = nan(10000,1);
late_diff_shuffled = nan(10000,1);
for n=1:10000
    shuffled_past_early_inds = randperm(length(all_early_angles),length(past_early));
    shuffled_future_early_inds = setdiff(1:length(all_early_angles),shuffled_past_early_inds);

    shuffled_past_late_inds = randperm(length(all_late_angles),length(past_late));
    shuffled_future_late_inds = setdiff(1:length(all_late_angles),shuffled_past_late_inds);

    early_diff_shuffled(n) = mean(all_early_angles(shuffled_past_early_inds))-mean(all_early_angles(shuffled_future_early_inds));
    late_diff_shuffled(n) = mean(all_late_angles(shuffled_past_late_inds))-mean(all_late_angles(shuffled_future_late_inds));
end
shuffled_interaction = early_diff_shuffled - late_diff_shuffled;

early_pval = (sum(early_diff_shuffled>abs(early_diff)) + sum(early_diff_shuffled < -1*abs(early_diff)) + 1)/(10001);
late_pval = (sum(late_diff_shuffled>abs(late_diff)) + sum(late_diff_shuffled < -1*abs(late_diff)) + 1)/(10001);
interaction_pval = (sum(shuffled_interaction>abs(real_interaction)) + sum(shuffled_interaction < -1*abs(real_interaction)) + 1)/(10001);

%% Alternative: 2 x 2 anova
% 
% data_table_sub = table();
% data_table_sub.rates = [future_early; past_early; future_late; past_late];
% data_table_sub.group = [zeros(size(future_early)); zeros(size(past_early)); ones(size(future_late)); ones(size(past_late))];
% data_table_sub.replay_class = [ones(size(future_early)); 2*ones(size(past_early)); ones(size(future_late)); 2*ones(size(past_late))];
% groupOrder = ["1","2"];
% namedGroup = categorical(data_table_sub.group,0:1,groupOrder);
% directionOrder = ["Reverse","Forward"];
% namedDirection = categorical(data_table_sub.replay_class,[2,1],directionOrder);
% data_sub = data_table_sub.rates;
% 
% % Anova effects:
% [p,tbl,stats] = anovan(data_sub',{namedGroup,namedDirection},'model','interaction','varnames',{'namedGroup','namedDirection'});
% [results,~,~,gnames] = multcompare(stats,"Dimension",[1 2]);

