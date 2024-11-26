function customBoxchart(data_for_boxplot,namedLaserState,namedDirection,fig_title,fig_ylim,fig_ylabel,figure_stats)
b=boxchart(namedLaserState,data_for_boxplot,'GroupByColor',namedDirection,'MarkerStyle','none');

b(1).BoxFaceColor = [0 0 0];
b(1).MarkerColor = [0 0 0];
b(2).BoxFaceColor = [1 0 0];
b(2).MarkerColor =  [1 0 0];
title(fig_title,'Interpreter','none')
ylabel(fig_ylabel)
if ~ isempty(fig_ylim)
    ylim(fig_ylim);
end
set(gca,'FontSize',12)

text_box_x = 0.2;
% text_box_x = gca().XLim(1) + 0.1*(gca().XLim(2) - gca().XLim(1));
for i = 1:length(figure_stats)
    text_box_y = gca().YLim(2) - i*0.1*(gca().YLim(2) - gca().YLim(1));
    text(text_box_x,text_box_y,num2str(figure_stats(i),2))
hold on
end