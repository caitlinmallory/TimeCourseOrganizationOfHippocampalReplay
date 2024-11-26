function simple_barplot(data, group_colors, figure_xticklabels, figure_ylim, figure_ylabel, figure_stats)


y_00 = data{1};
y_01 = data{2};


% remove nans
y_00 = y_00(~isnan(y_00));
y_01 = y_01(~isnan(y_01));


y_00_mean = nanmean(y_00);
y_01_mean = nanmean(y_01);


y_00_sem = nanstd(y_00)/sqrt(length(y_00));
y_01_sem  = nanstd(y_01)/sqrt(length(y_01));


y = [y_00_mean; y_01_mean];
err = [y_00_sem; y_01_sem; ];


hb = bar(y,'FaceColor','flat'); % get the bar handles
hold on;
for k = 1:size(y,2)
    % get x positions per group
    xpos = hb(k).XData + hb(k).XOffset;
    % draw errorbar
    errorbar(xpos, y(:,k), err(:,k), 'LineStyle', 'none', ...
        'Color', 'k', 'LineWidth', 1);
end
hb.CData(1,:) = group_colors(1,:);
hb.CData(2,:) = group_colors(2,:);

% Set Axis properties
set(gca,'xticklabel',figure_xticklabels);
if ~isempty(figure_ylim)
ylim(figure_ylim);
end
ylabel(figure_ylabel,'Interpreter','none')
xticklabels(figure_xticklabels);
set(gca,'FontSize',12)
title(figure_stats,'Interpreter','none')


