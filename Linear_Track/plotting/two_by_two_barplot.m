function two_by_two_barplot(data, group_colors, figure_xticklabels, figure_ylim, figure_ylabel, figure_legend, figure_stats, figure_title)


y_00 = data{1};
y_01 = data{2};
y_10 = data{3};
y_11 = data{4};

% remove nans
y_00 = y_00(~isnan(y_00));
y_01 = y_01(~isnan(y_01));
y_10 = y_10(~isnan(y_10));
y_11 = y_11(~isnan(y_11));

y_00_mean = nanmean(y_00);
y_01_mean = nanmean(y_01);
y_10_mean = nanmean(y_10);
y_11_mean = nanmean(y_11);

y_00_sem = nanstd(y_00)/sqrt(length(y_00));
y_01_sem  = nanstd(y_01)/sqrt(length(y_01));
y_10_sem  = nanstd(y_10)/sqrt(length(y_10));
y_11_sem  = nanstd(y_11)/sqrt(length(y_11));


y = [y_00_mean y_01_mean; y_10_mean y_11_mean];
err = [y_00_sem y_01_sem; y_10_sem y_11_sem];


hb = bar(y); % get the bar handles
hold on;
for k = 1:size(y,2)
    % get x positions per group
    xpos = hb(k).XData + hb(k).XOffset;
    % draw errorbar
    errorbar(xpos, y(:,k), err(:,k), 'LineStyle', 'none', ...
        'Color', 'k', 'LineWidth', 1);
end


hb(1).FaceColor = 'flat';
hb(2).FaceColor = 'flat';
hb(1).EdgeColor = 'none';
hb(2).EdgeColor = 'none';
hb(1).CData(1,:) = group_colors{1,1};
hb(1).CData(2,:) = group_colors{2,1};
hb(2).CData(1,:) = group_colors{1,2};
hb(2).CData(2,:) = group_colors{2,2};


% Set Axis properties
set(gca,'xticklabel',figure_xticklabels);
if ~isempty(figure_ylim)
ylim(figure_ylim);
end

if ~isempty(figure_ylabel)
ylabel(figure_ylabel,'Interpreter','none')
end

if ~isempty(figure_legend)
legend(figure_legend)
end
set(gca,'FontSize',12)

if ~isempty(figure_title)
title(figure_title,'Interpreter','none')
end

if ~isempty(figure_stats)
end

text_box_x = gca().XLim(1) + 0.1*(gca().XLim(2) - gca().XLim(1));
for i = 1:length(figure_stats)
    text_box_y = gca().YLim(2) - i*0.1*(gca().YLim(2) - gca().YLim(1));
    text(text_box_x,text_box_y,num2str(figure_stats(i),2))
hold on
end