function custom_boxplot_paired(data,laser_state,fig_ylabel,fig_ylim,fig_title,plot_text)


data_table = table(); data_table.data=data; data_table.laser_state = laser_state;
laserStateOrder = ["Laser-off","Laser-on"];
namedLaserState = categorical(data_table.laser_state,0:1,laserStateOrder);
b = boxchart(data_table.data,'GroupByColor',namedLaserState,'MarkerStyle','none');
b(1).BoxFaceColor = [0 0 0]; b(1).MarkerColor = [0 0 0];
b(2).BoxFaceColor = [1 0 0]; b(2).MarkerColor = [1 0 0];

differences = data(laser_state==0)-data(laser_state==1);
if sum(~isnan(differences)) > 3
h = lillietest(differences);
else
h = 0;
end
if h == 1
[p,~,z]=signrank(data(laser_state==0),data(laser_state==1));
    stat_test = 'WSR';
else
    [~,p,~,test_stat_val] = ttest(data(laser_state==0),data(laser_state==1));
    stat_test = 'ttest';
end

ylabel(fig_ylabel)
% legend
title(fig_title,'Interpreter','none')
xticklabels({})
xticks([])

if ~isempty(fig_ylim)
    ylim(fig_ylim);
else
    % Loop through each boxchart object
upperbound = []; 
lowerbound = []; 
for i = 1:numel(b)
    % Compute outlier bounds: box edges +/- (1.5 * IQR)
    groups = findgroups(b(i).XData);
    qtile.lower = splitapply(@(x)quantile(x,0.25),b(i).YData,groups);
    qtile.upper = splitapply(@(x)quantile(x,0.75),b(i).YData,groups);
    iqr = qtile.upper - qtile.lower; 
    upperbound = [upperbound; qtile.upper + 1.5*iqr]; %#ok<*AGROW>
    lowerbound = [lowerbound; qtile.lower - 1.5*iqr]; 
end
ybound = [min(lowerbound), max(upperbound)];
% Set y axis limit
ylim(ybound)
end
fig_ylim = gca().YLim;

fig_xlim = 0.05;
n = sum(~isnan(data_table.data(data_table.laser_state==0)) & ~isnan(data_table.data(data_table.laser_state==1)));

if plot_text==1
if strcmp(stat_test, 'WSR') && isfield(z,'zval')
text(fig_xlim,fig_ylim(2)-0.1*fig_ylim(2),[stat_test ' p=' num2str(p,2) ', z=' num2str(z.zval,2)])
else
 text(fig_xlim,fig_ylim(2)-0.1*fig_ylim(2),[stat_test ' p=' num2str(p,2) ', t=' num2str(test_stat_val.tstat,2)])
end

text(fig_xlim,fig_ylim(2)-0.2*fig_ylim(2),['n=' num2str(n,2)])
end