

%%
% setup_properties_table % contains a list of all properties you might want to look at/plot

control_experimental_animals_to_include = [1];
log_transform_data = 0;


upper_limit_time_since_reward_zone_entry = 10; % inf is most permissive
lower_limit_time_since_reward_zone_entry = 0; % -inf is most permissive
events = replay(replay.time_since_real_drink_onset >= lower_limit_time_since_reward_zone_entry & replay.time_since_real_drink_onset <= upper_limit_time_since_reward_zone_entry,:);
% Forward versus reverse:
events.replay_class(events.past==1) = 2;
events.replay_class(events.future==1) = 1;


events.group(events.laser_state_binary==0 & events.replay_class==2) = 1; % reverse or past
events.group(events.laser_state_binary==1 & events.replay_class==2) = 2; % reverse or past
events.group(events.laser_state_binary==0 & events.replay_class==1) = 3; % forward or future
events.group(events.laser_state_binary==1 & events.replay_class==1) = 4; % forward or future

% events = events(ismember(events.control_experimental_flag,control_experimental_animals_to_include) & events.time_since_reward_zone_entry > lower_limit_time_since_reward_zone_entry & events.time_since_reward_zone_entry < upper_limit_time_since_reward_zone_entry,:);

plot_individual_graphs = 0;
plot_boxplots=1;
plot_barplots=0;
plot_pdfs=0;
plot_cdfs=0;

fig_path = '/home/caitlin/Insync/caitlinmallory@berkeley.edu/Google Drive/Processed_Data/';


properties_to_plot = {'slope_distance','spike_density_power','ripple_power','duration','distance','dispersion','meanAngDisplacement_futPath','meanAngDisplacement_pastPath'};

%plot replay mosiaic
numPanels = 36;
subplot_xDim = 8;
subplot_yDim = (ceil((numPanels)/subplot_xDim));
panel_count = 1;

if plot_boxplots==1
    fig_boxplot = figure(); fig_boxplot.Position = [50 50 1800 900]; ha_fig_boxplot = tight_subplot(subplot_yDim,subplot_xDim,[.05 .05],[.05 .05],[.05 .05]);
end
if plot_barplots == 1
    fig_barplot  = figure(); fig_barplot.Position = [50 50 1800 900]; ha_fig_barplot = tight_subplot(subplot_yDim,subplot_xDim,[.05 .05],[.05 .05],[.05 .05]);
end
if plot_pdfs == 1
    fig_pdf = figure(); fig_pdf.Position = [50 50 1800 900]; ha_fig_pdf = tight_subplot(subplot_yDim,subplot_xDim,[.05 .05],[.05 .05],[.05 .05]);
end
if plot_cdfs == 1
    fig_cdf = figure();  fig_cdf.Position = [50 50 1800 900]; ha_fig_cdf = tight_subplot(subplot_yDim,subplot_xDim,[.05 .05],[.05 .05],[.05 .05]);
end



for i = 1:length(properties_to_plot)

   % properties_table_row = find(strcmp(properties_to_plot{i},properties.names)==1);
    num_panel = i;

    data = [...
        {events.(properties_to_plot{i})(events.laser_state_binary==0 & events.replay_class==2)};
        {events.(properties_to_plot{i})(events.laser_state_binary==0 & events.replay_class==1)};
        {events.(properties_to_plot{i})(events.laser_state_binary==1 & events.replay_class==2)};
        {events.(properties_to_plot{i})(events.laser_state_binary==1 & events.replay_class==1)}];

    data{1}(isinf(data{1}) | isnan(data{1})) = [];
    data{2}(isinf(data{2}) | isnan(data{2})) = [];
    data{3}(isinf(data{3}) | isnan(data{3})) = [];
    data{4}(isinf(data{4}) | isnan(data{4})) = [];

    if log_transform_data==1
        data{1} = log10(data{1});
        data{2} = log10(data{2});
        data{3} = log10(data{3});
        data{4} = log10(data{4});
    end

    data_table_sub = events;
    data_table_sub(isinf(data_table_sub.(properties_to_plot{i})),:) = [];
    data_sub = data_table_sub.(properties_to_plot{i});
    if log_transform_data==1
        data_sub = log10(data_sub);
    end

    laserStateOrder = ["Off","On"];
    namedLaserState = categorical(data_table_sub.laser_state_binary,0:1,laserStateOrder);
    directionOrder = ["Reverse","Forward"];
    namedDirection = categorical(data_table_sub.replay_class,[2,1],directionOrder);
    %groupOrder = ["Off_rev","On_rev","Off_for","On_for"];


    % Anova effects:
    p_anova = anovan(data_sub,{namedLaserState,namedDirection},'model','interaction','varnames',{'namedLaserState','namedDirection'},'display','off');

    if ~isempty(data{3})
    WRS.p1 = ranksum(data{1},data{3});
    WRS.p2 = ranksum(data{2},data{4});
    else
        WRS.p1 = nan;
        WRS.p2 = nan;
    end
    WRS.p3 = ranksum(data{1},data{2});


    if plot_boxplots == 1
        %BoxChart:
        sub_axes = ha_fig_boxplot(num_panel);
        axes(sub_axes);
        b=customBoxchart(data_sub,namedLaserState,namedDirection,sub_axes, p_anova, WRS);
        if plot_individual_graphs==1
            figure(); sub_axes = gca;
            b=customBoxchart(data_sub,namedLaserState,namedDirection,sub_axes, p_anova, WRS);
            saveas(gcf,fullfile(fig_path,[properties_to_plot{i} 'boxplot']),'jpeg');
        end
    end

    if plot_barplots==1
        % BarPlot
        sub_axes = ha_fig_barplot(num_panel);
        axes(sub_axes);
        two_by_two_barplot(data, [0.4940 0.1840 0.5560; .4660 0.6740 0.1880], {'Laser-Off'; 'Laser-On'}, properties.ylims{properties_table_row}(1,1:2), properties.ylabels{properties_table_row}, {}, [],properties.titles{properties_table_row});
        if plot_individual_graphs==1
            figure();
            two_by_two_barplot(data, [0.4940 0.1840 0.5560; .4660 0.6740 0.1880], {'Laser-Off'; 'Laser-On'}, properties.ylims{properties_table_row}(1,1:2), properties.ylabels{properties_table_row}, {}, [],properties.titles{properties_table_row});
            saveas(gcf,fullfile(fig_path,[properties_to_plot{i} 'barplot']),'jpeg');
        end
    end

    if plot_pdfs==1
        % PDF/Histogram plot:
        sub_axes = ha_fig_pdf(num_panel);
        axes(sub_axes);
        customPDF(data,properties(properties_table_row,:))
        if plot_individual_graphs==1
            figure()
            customPDF(data,properties(properties_table_row,:))
            saveas(gcf,fullfile(fig_path,[properties_to_plot{i} 'pdfplot']),'jpeg');
        end
    end

    if plot_cdfs==1

        % CDF plot:
        sub_axes = ha_fig_cdf(num_panel);
        axes(sub_axes);
        customCDF(data,properties(properties_table_row,:))
        if plot_individual_graphs==1
            figure()
            customCDF(data,properties(properties_table_row,:))
            saveas(gcf,fullfile(fig_path,[properties_to_plot{i} 'cdffplot']),'jpeg');
        end
    end
end

%% Make pretty figure for publication

properties_to_plot = {...
    'ripple_power'; ...
    'dispersion';...
    'distance';...
    'slope_distance'};

% fig_xlims = [...
%     [-2 4];...
%     [-2 4]
%     [-2 4]];

fig_xlims = [...
    [0 15];...
    [0 15]
    [0 50]
    [0 500]];

fig_titles = {...
    'Ripple power (zscored)';...
    'Dispersion';...
    'Distance';...
    'Slope using distance'};


figure('Position',[653 502 800 200])
tiledlayout(1,4,'TileSpacing','tight','Padding','tight')
for i = 1:length(properties_to_plot)
    nexttile
    num_panel = i;
    data_table_sub = events;
    data_table_sub(isinf(data_table_sub.(properties_to_plot{i})),:) = [];

    data = [...
        {data_table_sub.(properties_to_plot{i})(data_table_sub.laser_state_binary==0 & data_table_sub.replay_class==2)};
        {data_table_sub.(properties_to_plot{i})(data_table_sub.laser_state_binary==1 & data_table_sub.replay_class==2)};
        {data_table_sub.(properties_to_plot{i})(data_table_sub.laser_state_binary==0 & data_table_sub.replay_class==1)};
        {data_table_sub.(properties_to_plot{i})(data_table_sub.laser_state_binary==1 & data_table_sub.replay_class==1)}];
    customCDF(data,fig_titles{i},fig_xlims(i,:))
box off
end
set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
export_fig(fullfile(fig_path,'pretty_summary_cdf'),'-jpeg')
saveas(gcf,fullfile(fig_path,'pretty_summary_cdf'),'pdf')
%% Make pretty figure for publication-- barplots

figure('Position',[653 502 800 200])
tiledlayout(1,length(properties_to_plot))
for i = 1:length(properties_to_plot)
    sub_axes = nexttile;

    data_table_sub = events;
    data_table_sub(isinf(data_table_sub.(properties_to_plot{i})),:) = [];

    data = [...
        {data_table_sub.(properties_to_plot{i})(data_table_sub.laser_state_binary==0 & data_table_sub.replay_class==1)};
        {data_table_sub.(properties_to_plot{i})(data_table_sub.laser_state_binary==1 & data_table_sub.replay_class==1)};
        {data_table_sub.(properties_to_plot{i})(data_table_sub.laser_state_binary==0 & data_table_sub.replay_class==2)};
        {data_table_sub.(properties_to_plot{i})(data_table_sub.laser_state_binary==1 & data_table_sub.replay_class==2)}];

    laserStateOrder = ["Off","On"];
    namedLaserState = categorical(data_table_sub.laser_state_binary,0:1,laserStateOrder);
    directionOrder = ["Forward","Reverse"];
    namedDirection = categorical(data_table_sub.replay_class,[1,2],directionOrder);
    groupOrder = ["Off_for","On_for","Off_rev","On_rev"];
    namedGroup = categorical(data_table_sub.group,1:4,groupOrder);
    % Anova effects:
    p_anova = anovan(data_table_sub.(properties_to_plot{i}),{namedLaserState,namedDirection},'model','interaction','varnames',{'namedLaserState','namedDirection'},'display','off');
   
    if ~isempty(data{2})
    WRS.p1 = ranksum(data{1},data{2});
    else
        WRS.p1 = nan;
    end
    if ~isempty(data{4})
    WRS.p2 = ranksum(data{3},data{4});
    else
        WRS.p2 = nan;
    end
    WRS.p3 = ranksum(data{1},data{3});


    rev_color_off = [0.4940 0.1840 0.5560];
    for_color_off = [0.4660 0.6740 0.1880];
    transparency_pcnt = 0.5;
    rev_color_on = 1-transparency_pcnt.*(1-rev_color_off);
    for_color_on = 1-transparency_pcnt.*(1-for_color_off);


    two_by_two_barplot(data, [{for_color_off} {for_color_on}; {rev_color_off} {rev_color_on}], {'For'; 'Rev'}, [],[],[],[],[]);
    y_range = gca().YLim(2)-gca().YLim(1);
    text(0.55, gca().YLim(1) + 0.95*y_range,['laser: ' num2str(p_anova(1),2)])
    text(0.55, gca().YLim(1) + 0.88*y_range,['dir: ' num2str(p_anova(2),2)])
    text(0.55, gca().YLim(1) + 0.82*y_range,['int: ' num2str(p_anova(3),2)])
    box off
end
set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
export_fig(fullfile(fig_path,'pretty_summary_barplot'),'-jpeg')
saveas(gcf,fullfile(fig_path,'pretty_summary_barplot'),'pdf')


%% Make pretty 2 by 2 boxplots
figure('Position',[653 502 861 237])
tiledlayout(1,3)
for i = 1:length(properties_to_plot)
    sub_axes = nexttile;
       data_table_sub = events;
    data_table_sub(isinf(data_table_sub.(properties_to_plot{i})),:) = [];

    data = [...
        {data_table_sub.(properties_to_plot{i})(data_table_sub.laser_state_binary==0 & data_table_sub.replay_class==2)};
        {data_table_sub.(properties_to_plot{i})(data_table_sub.laser_state_binary==1 & data_table_sub.replay_class==2)};
        {data_table_sub.(properties_to_plot{i})(data_table_sub.laser_state_binary==0 & data_table_sub.replay_class==1)};
        {data_table_sub.(properties_to_plot{i})(data_table_sub.laser_state_binary==1 & data_table_sub.replay_class==1)}];
    laserStateOrder = ["Off","On"];
    namedLaserState = categorical(data_table_sub.laser_state_binary,0:1,laserStateOrder);
    directionOrder = ["Forward","Reverse"];
    namedDirection = categorical(data_table_sub.replay_class,[1,2],directionOrder);
    groupOrder = ["Off_for","On_for","Off_rev","On_rev"];
    namedGroup = categorical(data_table_sub.group,[3 4 1 2],groupOrder);

    % Anova effects:
    p_anova = anovan(data_table_sub.(properties_to_plot{i}),{namedLaserState,namedDirection},'model','interaction','varnames',{'namedLaserState','namedDirection'},'display','off');
   if ~isempty(data{2})
    WRS.p1 = ranksum(data{1},data{2});
   else
       WRS.p1 = nan;
   end
   if ~isempty(data{4})
    WRS.p2 = ranksum(data{3},data{4});
   else
       WRS.p2 = nan;
   end
    WRS.p3 = ranksum(data{1},data{3});

    % b=customBoxchart2(data_sub.(properties_to_plot{i}),namedGroup,properties(properties_table_row,:),1,sub_axes, p_anova, WRS);
    b=customBoxchart2(data_table_sub.(properties_to_plot{i}),namedGroup,sub_axes, p_anova, WRS);
   
    groupOrder = ["Off_for","On_for","Off_rev","On_rev"];
    namedGroup = categorical(data_table_sub.group,[3 4 1 2],groupOrder);
 
    % custom_boxplot_unpaired(data,[1 2],{'Off_for','Off_rev'},[{for_color_off}; {rev_color_off}],'',[],'',1);
    %     y_range = gca().YLim(2)-gca().YLim(1);
    %     text(0.55, gca().YLim(1) + 0.95*y_range,['laser: ' num2str(p_anova(1),2)])
    %     text(0.55, gca().YLim(1) + 0.88*y_range,['dir: ' num2str(p_anova(2),2)])
    %     text(0.55, gca().YLim(1) + 0.82*y_range,['int: ' num2str(p_anova(3),2)])
end
set(gcf, 'Color', 'white','Renderer','painters');
export_fig(fullfile(fig_path,'pretty_summary_boxplot'),'-jpeg')
export_fig(fullfile(fig_path,'pretty_summary_boxplot'),'-pdf')

%% Histograms with boxplots on top

figure('Position',[653 502 861 510])
tiledlayout(2,3)
for i = 1:length(properties_to_plot)

  
    data_table_sub = events;
    data_table_sub(isinf(data_table_sub.(properties_to_plot{i})),:) = [];

    data_01 = data_table_sub.(properties_to_plot{i})(data_table_sub.laser_state_binary==0 & data_table_sub.replay_class==1);
    data_02 = data_table_sub.(properties_to_plot{i})(data_table_sub.laser_state_binary==0 & data_table_sub.replay_class==2);
    data_11 = data_table_sub.(properties_to_plot{i})(data_table_sub.laser_state_binary==1 & data_table_sub.replay_class==1);
    data_12 = data_table_sub.(properties_to_plot{i})(data_table_sub.laser_state_binary==1 & data_table_sub.replay_class==2);
    if log_transform_data == 1
    data_01 = log10(data_01);
    data_02 = log10(data_02);
    data_11 = log10(data_11);
    data_12 = log10(data_12);
    end


    nexttile(i)

    clear A
    A.For_Off = data_01;
    A.For_On = data_11;

    rev_color_off = [0.4940 0.1840 0.5560];
    for_color_off = [0.4660 0.6740 0.1880];
    transparency_pcnt = 0.5;
    rev_color_on = 1-transparency_pcnt.*(1-rev_color_off);
    for_color_on = 1-transparency_pcnt.*(1-for_color_off);

    colors_forward_on_off = [for_color_off; for_color_on];
    colormap(colors_forward_on_off)

    if i == length(properties_to_plot)
        nhist(A,'samebins',1,'noerror',1,'pdf','box',1,'color','colormap')
    else
        nhist(A,'samebins',1,'noerror',1,'pdf','box',1,'color','colormap','nolegend')
    end

    % Plot reverse separately:
    nexttile(i+length(properties_to_plot))
    clear A
    A.Rev_Off = data_02;
    A.Rev_On = data_12;
    colors_reverse_on_off = [rev_color_off; rev_color_on];
    colormap(colors_reverse_on_off)

    if i == length(properties_to_plot)
        nhist(A,'samebins',1,'noerror',1,'pdf','box',1,'color','colormap')
    else
        nhist(A,'samebins',1,'noerror',1,'pdf','box',1,'color','colormap','nolegend')
    end

end
set(gcf, 'Color', 'white','Renderer','painters');
export_fig(fullfile(fig_path,'pretty_summary_hist_with_boxplot'),'-jpeg')
export_fig(fullfile(fig_path,'pretty_summary_hist_with_boxplot'),'-pdf')






%%


if plot_boxplots==1
    saveas(fig_boxplot,fullfile(fig_path,'Replay properties boxplot'),'jpeg')
    saveas(fig_boxplot,fullfile(fig_path,'Replay properties boxplot'),'pdf')
end
if plot_barplots==1
    saveas(fig_barplot,fullfile(fig_path,'Replay properties barplot'),'jpeg')
    saveas(fig_barplot,fullfile(fig_path,'Replay properties barplot'),'pdf')
end
if plot_pdfs==1
    saveas(fig_pdf,fullfile(fig_path,'Replay properties hist'),'jpeg')
    saveas(fig_pdf,fullfile(fig_path,'Replay properties hist'),'pdf')    
end
if plot_cdfs==1
    saveas(fig_cdf,fullfile(fig_path,'Replay properties cdf'),'jpeg')
    saveas(fig_cdf,fullfile(fig_path,'Replay properties cdf'),'pdf')    
end

function b =  customBoxchart(data_for_boxplot,groupingVar1,groupingVar2,sub_axes,p_anova,WRS)
b=boxchart(groupingVar1,data_for_boxplot,'GroupByColor',groupingVar2);

b(1).BoxFaceColor = [0.4940 0.1840 0.5560];
b(1).BoxFaceAlpha = 0.5;
b(1).MarkerColor = [0.4940 0.1840 0.5560];

b(2).BoxFaceColor = [.4660 0.6740 0.1880];
b(2).BoxFaceAlpha = 0.5;
b(2).MarkerColor =  [.4660 0.6740 0.1880];
set(gca,'FontSize',12)

if ~isempty(p_anova)

    text(0.05,sub_axes.YLim(2)-(sub_axes.YLim(2)- sub_axes.YLim(1))/10 ,num2str(p_anova(1),2));
    text(0.05,sub_axes.YLim(2)-2*((sub_axes.YLim(2)- sub_axes.YLim(1))/10),num2str(p_anova(2),2));
    text(0.05,sub_axes.YLim(2)-3*((sub_axes.YLim(2)- sub_axes.YLim(1))/10),num2str(p_anova(3),2));
end

if ~isempty(WRS)
    text(0.05,sub_axes.YLim(2)-7*(sub_axes.YLim(2)- sub_axes.YLim(1))/10 ,num2str(WRS.p1,2));
    text(0.05,sub_axes.YLim(2)-8*((sub_axes.YLim(2)- sub_axes.YLim(1))/10),num2str(WRS.p2,2));
    text(0.05,sub_axes.YLim(2)-9*((sub_axes.YLim(2)- sub_axes.YLim(1))/10),num2str(WRS.p3,2));
end
end

function b =  customBoxchart2(data_for_boxplot,groups,sub_axes,p_anova,WRS)
keyboard
b=boxchart(data_for_boxplot,'GroupByColor',groups,'MarkerStyle','none');

if size(b,1)>2
b(3).BoxFaceColor = [0.4940 0.1840 0.5560];
b(3).BoxFaceAlpha = 0.7;
b(3).MarkerColor = [0.4940 0.1840 0.5560];

b(4).BoxFaceColor = [0.4940 0.1840 0.5560];
b(4).BoxFaceAlpha = 0.2;
b(4).MarkerColor = [0.4940 0.1840 0.5560];


b(1).BoxFaceColor = [.4660 0.6740 0.1880];
b(1).BoxFaceAlpha = 0.7;
b(1).MarkerColor =  [.4660 0.6740 0.1880];

b(2).BoxFaceColor = [.4660 0.6740 0.1880];
b(2).BoxFaceAlpha = 0.2;
b(2).MarkerColor =  [.4660 0.6740 0.1880];
set(gca,'FontSize',12)
end

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


if ~isempty(p_anova)

    text(0.05,sub_axes.YLim(2)-(sub_axes.YLim(2)- sub_axes.YLim(1))/10 ,num2str(p_anova(1),2));
    text(0.05,sub_axes.YLim(2)-2*((sub_axes.YLim(2)- sub_axes.YLim(1))/10),num2str(p_anova(2),2));
    text(0.05,sub_axes.YLim(2)-3*((sub_axes.YLim(2)- sub_axes.YLim(1))/10),num2str(p_anova(3),2));
end

if ~isempty(WRS)
    text(0.05,sub_axes.YLim(2)-7*(sub_axes.YLim(2)- sub_axes.YLim(1))/10 ,num2str(WRS.p1,2));
    text(0.05,sub_axes.YLim(2)-8*((sub_axes.YLim(2)- sub_axes.YLim(1))/10),num2str(WRS.p2,2));
    text(0.05,sub_axes.YLim(2)-9*((sub_axes.YLim(2)- sub_axes.YLim(1))/10),num2str(WRS.p3,2));
end
end

function customPDF(data)
r_0 = data{1}; f_0 = data{2}; r_1 = data{3}; f_1 = data{4};

p = hist(f_0);
plot(hist_axis,p./sum(p),'color',"#77AC30",'LineWidth',2);
hold on
p = hist(f_1);
plot(hist_axis,p./sum(p),'color',"#77AC30",'LineWidth',2,'LineStyle',":");
hold on
p = hist(r_0);
plot(hist_axis,p./sum(p),'color',"#7E2F8E",'LineWidth',2);
hold on
p = hist(r_1);
plot(hist_axis,p./sum(p),'color',"#7E2F8E",'LineWidth',2,'LineStyle',":");
title(properties_sub.titles{:},'Interpreter','none')
xlabel(properties_sub.ylabels{:})
ylabel('Frequency')
set(gca,'FontSize',12)


end

function customCDF(data,fig_title,fig_xlim)
r_0 = data{1}; r_1 = data{2}; f_0 = data{3}; f_1 = data{4};
if ~isempty(f_0)
[p,x] = ecdf(f_0);
plot(x,p,'color',"#77AC30",'LineWidth',2);
hold on
end
if ~isempty(f_1)
[p,x] = ecdf(f_1);
plot(x,p,'color',"#77AC30",'LineWidth',2,'LineStyle',":");
hold on
end
if ~isempty(r_0)
[p,x] = ecdf(r_0);
plot(x,p,'color',"#7E2F8E",'LineWidth',2);
hold on
end
if ~isempty(r_1)
[p,x] = ecdf(r_1);
plot(x,p,'color',"#7E2F8E",'LineWidth',2,'LineStyle',":");
end
xlabel(fig_title)
ylabel('Cumulative frequency')
set(gca,'FontSize',10)

if ~isempty(r_1)
[h,kstest_p_reverse] = kstest2(r_0,r_1);
else
    kstest_p_reverse = nan;
end
if ~isempty(f_1)
[h,kstest_p_forward] = kstest2(f_0,f_1);
else
    kstest_p_forward = nan;
end

if ~isempty(fig_xlim)
    xlim(fig_xlim)
end


x_range = gca().XLim(2)-gca().XLim(1);
text(gca().XLim(1) + 0.7*x_range, 0.2,num2str(kstest_p_forward,2),'color',"#77AC30")
text(gca().XLim(1) + 0.7*x_range, 0.1,num2str(kstest_p_reverse,2),'color',"#7E2F8E")
end



