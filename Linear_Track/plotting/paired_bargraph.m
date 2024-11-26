function fig = paired_bargraph(measures, colors, nCats, ylimit, figure_xticklabels, figure_title, pval)
nDatas = size(measures,1);
%% Plot


x = [1;2];
mean_data = nanmean(measures,1);
sem_data = nanstd(measures,1)./sqrt(length(measures));

fig = figure();
b = bar(x,mean_data);
b.FaceColor = [0.5 0.5 0.5]
hold on
er = errorbar(x,mean_data,sem_data,sem_data);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 1;

for i = 1:length(measures)
    line([1 2],measures(i,:),'Color', colors(i,:),  'Marker', '.', 'MarkerSize', 10, 'LineWidth', 2);
    hold on
end

ylim(ylimit)
title([figure_title ' p=' num2str(pval,2)])
xticklabels(figure_xticklabels)
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12)