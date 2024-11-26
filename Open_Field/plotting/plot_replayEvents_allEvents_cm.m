load Experiment_Information
load Analysis_Information
load Position_Data
load replayEvents.mat

colormap_choice = cbrewer2('Blues');

figure('Position',[600 600 100 100])
for sessionNum_decoder = 1
    sessionNum = sessionNum_decoder;

    %load replays
    replay_dispersionThr = 3;
    replay_durationThr = 0;
    load_ind_replayEvents
    [replays,replays_NaN,replays_NaNremoved,replays_interp,~,~,~,~,replays_smoothed] = compute_replaysConcatenation(decoder_replay(sessionNum),ind_replay);

    set(gcf,'color','w')
    cplot(replays_NaNremoved(:,2),replays_NaNremoved(:,3),replays_NaNremoved(:,1),'linewidth',0.5)
%    plot_mazeProperties_cm
    set(gca,'xtick',[],'ytick',[]), axis square, box on, axis([1 numSpatialBins(2) 1 numSpatialBins(1)])
    colormap(colormap_choice) 
    hold on

    well_list = true_drink_periods_summary(sessionNum).goal_well_list;
    home_well = mode(well_list);
    Experiment_Information.homeLoc{sessionNum} = Experiment_Information.rewardLocs{sessionNum}(home_well,:);
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
saveas(gcf,['All_replays_disp' num2str(replay_dispersionThr,2)],'jpg');
saveas(gcf,['All_replays_disp' num2str(replay_dispersionThr,2)],'pdf');
end
