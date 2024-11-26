function fig = paired_boxplot(measures, colors, lineThicknesses, nCats, ylimit, figure_xticklabels, figure_title, pval)
nDatas = size(measures,1);
%% Plot

fig = figure();
boxplot(measures(1:nDatas, 1:nCats), 'Symbol', 'k.'); hold on;

set(findobj(gca,'type','line'),'linew',2,'color','k')

for i = 1:length(measures)
    line([1 2],measures(i,:),'Color', colors(i,:),  'Marker', '.', 'MarkerSize', 10, 'LineWidth', lineThicknesses(i));
    hold on
end

ylim(ylimit)
title([figure_title ' p=' num2str(pval,2)])
xticklabels(figure_xticklabels)
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12)