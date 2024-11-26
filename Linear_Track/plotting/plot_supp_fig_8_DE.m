%% Make figure for publication-- barplots

% NOTE: make sure to only analyze data from rats 1 2 3 6 for supp fig. 8D
% (rats expressing Jaws, GFP).
% for supp figure 8E, only analyze rats 4 and 5.

properties_to_plot = {...
    'replay_ripple_power'; ...
    'fraction_of_excitatory_cells_participating';...
    'range'};
titles = {'Ripple power','% Participation','Range (cm)'};
if control_experimental_animals_to_include == 1
    ylims = [2 6; 20 40; 120 160];
    xlims = [0 15; 0 60; 50 250];
    xticks_fig = [{0:5:15},{0:20:60},{50:100:250}];
elseif control_experimental_animals_to_include == 0
    ylims = [2 6; 20 40; 80 140];
    xlims = [0 15; 0 60; 50 250];
    xticks_fig = [{0:5:15},{0:20:60},{50:100:250}];
end
figure('Position',[653 502 450 250])
tiledlayout(2,length(properties_to_plot),'TileSpacing','tight')

colors = [.4660 0.6740 0.1880; 0.4940 0.1840 0.5560;];
transparency_pcnt = 0.5;
colors_2 = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
    [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];

for i = 1:length(properties_to_plot)
    sub_axes = nexttile(i);

    data_table_sub = events;
    if i == 2
        data_table_sub.(properties_to_plot{i}) = (data_table_sub.(properties_to_plot{i})).*100;
    end
    data_table_sub(isinf(data_table_sub.(properties_to_plot{i})),:) = [];

    data = [...
        {data_table_sub.(properties_to_plot{i})(data_table_sub.laser_state==0 & data_table_sub.replay_class==1)};
        {data_table_sub.(properties_to_plot{i})(data_table_sub.laser_state==1 & data_table_sub.replay_class==1)};
        {data_table_sub.(properties_to_plot{i})(data_table_sub.laser_state==0 & data_table_sub.replay_class==2)};
        {data_table_sub.(properties_to_plot{i})(data_table_sub.laser_state==1 & data_table_sub.replay_class==2)}];

    length(data{1})
    length(data{2})
    length(data{3})
    length(data{4})

    [p,h,z]=ranksum(data{1},data{2})
    [p,h,z]=ranksum(data{3},data{4})

    [h,p,z] =  kstest2(data{1},data{2})
    [h,p,z] =  kstest2(data{3},data{4})

    laserStateOrder = ["Off","On"];
    namedLaserState = categorical(data_table_sub.laser_state,0:1,laserStateOrder);
    directionOrder = ["Forward","Reverse"];
    namedDirection = categorical(data_table_sub.replay_class,[1,2],directionOrder);
    groupOrder = ["Off_for","On_for","Off_rev","On_rev"];
    namedGroup = categorical(data_table_sub.group,1:4,groupOrder);

%     p_anova = anovan( data_table_sub.(properties_to_plot{i}),{namedLaserState,namedDirection},'model','interaction','varnames',{'namedLaserState','namedDirection'},'display','on');
    real_diff = (nanmean(data_table_sub.(properties_to_plot{i})(data_table_sub.laser_state==0 & data_table_sub.direction == 2)) - ...
        nanmean(data_table_sub.(properties_to_plot{i})(data_table_sub.laser_state==1 & data_table_sub.direction == 2))) - ...
        (nanmean(data_table_sub.(properties_to_plot{i})(data_table_sub.laser_state==0 & data_table_sub.direction == 1)) - ...
        nanmean(data_table_sub.(properties_to_plot{i})(data_table_sub.laser_state==1 & data_table_sub.direction == 1))) ;

    %shuffle_labels
    n_laser_0_dir2 = height(data_table_sub.(properties_to_plot{i})(data_table_sub.laser_state==0 & data_table_sub.direction == 2));
    n_laser_0_dir1 = height(data_table_sub.(properties_to_plot{i})(data_table_sub.laser_state==0 & data_table_sub.direction == 1));

    reverse = data_table_sub(data_table_sub.direction==2,:);
    forward = data_table_sub(data_table_sub.direction==1,:);

    num_shuffles = 5000;
    group_label_shuffle_diffs = nan(num_shuffles,1);
    for j = 1:num_shuffles
        inds1 = randperm(height(reverse),n_laser_0_dir2);
        inds2 = setdiff(1:height(reverse),inds1);
        inds3 = randperm(height(forward),n_laser_0_dir1);
        inds4 = setdiff(1:height(forward),inds3);

        group_label_shuffle_diffs(j) = (nanmean(reverse.(properties_to_plot{i})(inds1))-nanmean(reverse.(properties_to_plot{i})(inds2))) - ...
            (nanmean(forward.(properties_to_plot{i})(inds3))-nanmean(forward.(properties_to_plot{i})(inds4)));
    end

    shuffle_pval_1_tailed = (sum(group_label_shuffle_diffs > real_diff) + 1)/(num_shuffles + 1)
    shuffle_pval_2_tailed = (sum(group_label_shuffle_diffs > real_diff | group_label_shuffle_diffs < (-1)*real_diff) + 1)/(num_shuffles + 1)

    % Anova effects:
    properties_table_row = find(strcmp(properties_to_plot{i},properties.names)==1);
    rev_color_off = [0.4940 0.1840 0.5560];
    for_color_off = [0.4660 0.6740 0.1880];
    transparency_pcnt = 0.5;
    rev_color_on = 1-transparency_pcnt.*(1-rev_color_off);
    for_color_on = 1-transparency_pcnt.*(1-for_color_off);
    box off

    x_mean = [nanmean(data{1}),nanmean(data{2}),nanmean(data{3}) nanmean(data{4})];
    x_err = [nanstd(data{1})/sqrt(nansum(~isnan(data{1}))),...
        nanstd(data{2})/sqrt(nansum(~isnan(data{2}))),...
        nanstd(data{3})/sqrt(nansum(~isnan(data{3}))),...
        nanstd(data{4})/sqrt(nansum(~isnan(data{4})))]

    b1 = bar(1,x_mean(1)); hold on;
    e1 = errorbar(1,x_mean(1),x_err(1),'k','linestyle','none');
    e1.CapSize = 4;
    hold on
    b2 = bar(2,x_mean(2)); hold on;
    e2 = errorbar(2,x_mean(2),x_err(2),'k','linestyle','none');
    e2.CapSize = 4;
    b3 = bar(4,x_mean(3)); hold on;
    e3 = errorbar(4,x_mean(3),x_err(3),'k','linestyle','none');
    e3.CapSize = 4;
    b4 = bar(5,x_mean(4)); hold on;
    e4 = errorbar(5,x_mean(4),x_err(4),'k','linestyle','none');
    e4.CapSize = 4;
    b1.FaceColor = colors(1,:);
    b1.EdgeColor = 'none';
    b2.FaceColor = colors_2(1,:);
    b2.EdgeColor = 'none';
    b3.FaceColor = colors(2,:);
    b3.EdgeColor = 'none';
    b4.FaceColor = colors_2(2,:);
    b4.EdgeColor = 'none';

    box off
    xticks([1.5 4.5])
    xticklabels({'For.','Rev.'})
    ylabel(titles{i})
    ylim(ylims(i,:))

    nexttile(i+3)
    f_0 = data{1}; f_1 = data{2}; r_0 = data{3}; r_1 = data{4};
    [p,x] = ecdf(f_0);
    plot(x,p,'color',"#77AC30",'LineWidth',2);
    hold on
    [p,x] = ecdf(f_1);
    plot(x,p,'color',"#77AC30",'LineWidth',2,'LineStyle',":");
    hold on
    [p,x] = ecdf(r_0);
    plot(x,p,'color',"#7E2F8E",'LineWidth',2);
    hold on
    [p,x] = ecdf(r_1);
    plot(x,p,'color',"#7E2F8E",'LineWidth',2,'LineStyle',":");
    %xlabel(fig_titles{i},'Interpreter','none')
    box off
    if i == 1
        ylabel('Cum. frequency')
    end
    xlim(xlims(i,:))
    xticks(xticks_fig{i})
    xlabel(titles{i})
    xtickangle(0)
end
set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
export_fig(fullfile(fig_path,'pretty_summary_barplot'),'-jpeg')
% export_fig(fullfile(fig_path,'pretty_summary_barplot'),'-pdf')
saveas(gcf,fullfile(fig_path,'pretty_summary_barplot'),'pdf')