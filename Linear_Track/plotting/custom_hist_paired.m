function custom_hist_paired(data,laser_state,hist_bins,fig_title)

differences = data(laser_state==0)-data(laser_state==1);
median_difference = nanmedian(differences);

if ~isempty(hist_bins)
[values,hist_bins] = hist(differences,hist_bins);
else
    [values,hist_bins] = hist(differences,20);
end
b = bar(hist_bins,values);
b.FaceColor = [0 0 0];
b.EdgeColor = [0 0 0];

if length(~isnan(differences)) >= 4 % lillietest must have at least 4 observations
h = lillietest(differences);
else 
    h = 0;
end
if h == 1
[p, ~, z] =signrank(data(laser_state==0),data(laser_state==1));
test = 'WSR';
else
    [~,p] = ttest(data(laser_state==0),data(laser_state==1));
    test = 'ttest';
end


% legend
title(fig_title,'Interpreter','none')

fig_xlim = gca().XLim;
fig_ylim = gca().YLim;
n = sum(~isnan(data(laser_state==0)) & ~isnan(data(laser_state==1)));
if strcmp(test, 'WSR') & isfield(z,'zval')
text((fig_xlim(2)-fig_xlim(1))*0.05 + fig_xlim(1),fig_ylim(2)-0.1*fig_ylim(2),[test ' p=' num2str(p,2) ', z=' num2str(z.zval,2)])
else
 text((fig_xlim(2)-fig_xlim(1))*0.05 + fig_xlim(1),fig_ylim(2)-0.1*fig_ylim(2),[test ' p=' num2str(p,2)])
end

text((fig_xlim(2)-fig_xlim(1))*0.05 + fig_xlim(1),fig_ylim(2)-0.2*fig_ylim(2),['n=' num2str(n,2)])
text((fig_xlim(2)-fig_xlim(1))*0.05 + fig_xlim(1),fig_ylim(2)-0.3*fig_ylim(2),['med diff=' num2str(median_difference,4)])

xline(0,'r','LineWidth',2)
xline(median_difference,'g','LineWidth',2)
