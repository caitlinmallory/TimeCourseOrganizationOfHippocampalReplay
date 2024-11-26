function custom_boxplot_unpaired(data,groupingVariableOrder,grouping_variable,colors,fig_ylabel,fig_ylim,fig_title,plot_text)


namedVariable = categorical(grouping_variable,[0,1],groupingVariableOrder);
b = boxchart(data,'GroupByColor',namedVariable,'MarkerStyle','none');
b(1).BoxFaceColor = colors(1,:); b(1).MarkerColor = colors(1,:);
b(2).BoxFaceColor = colors(2,:); b(2).MarkerColor = colors(2,:);

if length(data(grouping_variable==0))>3 && length(data(grouping_variable==1))>3
h1 = lillietest(data(grouping_variable==0));
h2 = lillietest(data(grouping_variable==1));
else
    h1 = 1;
    h2 = 1;
end


if h1 == 0 && h2 == 0
[~,p,~,test_stat]=ttest2(data(grouping_variable==0),data(grouping_variable==1));
test = 'ttest';
z=test_stat.tstat;
else
[p,h,z]=ranksum(data(grouping_variable==0),data(grouping_variable==1));
test = 'WRS';
z=z.zval;
end

ylabel(fig_ylabel,'interpreter','none')
% legend
title(fig_title,'interpreter','none')
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

laser_off_n = nansum(grouping_variable==0);
laser_on_n = nansum(grouping_variable==1);


fig_ylim = gca().YLim; fig_yrange = fig_ylim(2) - fig_ylim(1);
if plot_text
text(0.01,fig_ylim(1)+0.95*fig_yrange,[test ' p=' num2str(p,2), 'stat=' num2str(z,2)])
text(0.01,fig_ylim(1)+0.85*fig_yrange,[groupingVariableOrder{1} ' n=' num2str(laser_off_n)]);
text(0.01,fig_ylim(1)+0.75*fig_yrange,[groupingVariableOrder{2} ' n=' num2str(laser_on_n)]);
end
