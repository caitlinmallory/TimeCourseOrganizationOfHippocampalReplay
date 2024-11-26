load true_drink_periods.mat 
load Position_Data
load Experiment_Information.mat;

Position_Data(:,2) = Position_Data(:,2)./2;
Position_Data(:,3) = Position_Data(:,3)./2;

for sessionNum = 1:size(true_drink_periods_summary,2)
        
    Position_Data_sub = compute_dataTemporalConcatenation(Position_Data,[true_drink_periods_summary(sessionNum).true_drink_periods(1,1) true_drink_periods_summary(sessionNum).true_drink_periods(end,2)]);

    well_list = true_drink_periods_summary(sessionNum).goal_well_list;
    home_well = mode(well_list);
    Experiment_Information.homeLoc{sessionNum} = Experiment_Information.rewardLocs{sessionNum}(home_well,:);
    figure('Position',[600 600 100 100])
    plot(Position_Data_sub(:,2),Position_Data_sub(:,3),'color',[0.5 0.5 0.5],'LineWidth',0.25); hold on
    wells = (Experiment_Information.rewardLocs{sessionNum})./2;
    % mark the location of all wells with open circles
    for i = 1:length(wells)
        plot(wells(i,1),wells(i,2),'ok','markerSize',3,'LineWidth',1)
        hold on
    end

    % mark the location of the home well with a filled circle
    plot(wells(home_well,1),wells(home_well,2),'ok','markerSize',3,'MarkerFaceColor','k','LineWidth',1)

    barriers = Experiment_Information.barrierLocs{sessionNum};
    % draw in the barriers for this session
    for i = 1:6
        if ~isempty(barriers{i})
        line([barriers{i}(1)./2 barriers{i}(2)./2], [barriers{i}(3)./2 barriers{i}(4)./2],'LineWidth',1,'Color','k'); hold on;
        end
    end
   
axis square
xticks([])
yticks([])
box on
xlim([0 (Experiment_Information.maze_size)/2]);
ylim([0 (Experiment_Information.maze_size)/2]);
ax = gca;
% ax.LineWidth = 1;
rectangle('position',[0 0 (Experiment_Information.maze_size)/2 (Experiment_Information.maze_size)/2])
set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
saveas(gcf,['Segment_' num2str(sessionNum)],'jpg')
saveas(gcf,['Segment_' num2str(sessionNum)],'pdf')
end
% save('Experiment_Information.mat','Experiment_Information')