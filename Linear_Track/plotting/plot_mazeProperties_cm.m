

% if isfield(Experiment_Information,'rewardLocs')
%     if ~isempty(barrierConfig_IDs{sessionNum}) && ~isempty(rewardLocs{sessionNum})
%         rewardLocs_scaled = compute_locsToBins(rewardLocs{sessionNum},numSpatialBins,x_edges,y_edges);
%         [~,barriers_mat] = load_idealizedBarrierLocations(rewardLocs_scaled,barrierConfig_IDs{sessionNum});
%         hold on, plot(barriers_mat(:,1),barriers_mat(:,2),'k','linewidth',4), hold off
%     end
% else
%    hold on, plot_barrierLocs(barrierLocs{sessionNum},1,1,numSpatialBins,4), hold off
% end


if isfield(Experiment_Information,'homeLoc') && ~isempty(Experiment_Information.homeLoc{sessionNum})
    homeLoc_scaled = compute_locsToBins(Experiment_Information.homeLoc{sessionNum},numSpatialBins,x_edges,y_edges);                 
   for home_location = 1:size(homeLoc_scaled,1)
    hold on, plot(homeLoc_scaled(home_location,1),homeLoc_scaled(home_location,2),'k.','markersize',20), hold off
   end
end
if isfield(Experiment_Information,'rewardLocs') && ~isempty(Experiment_Information.rewardLocs{sessionNum})
    rewardLocs_scaled = compute_locsToBins(Experiment_Information.rewardLocs{sessionNum},numSpatialBins,x_edges,y_edges);                 
%     hold on, plot(rewardLocs_scaled(:,1),rewardLocs_scaled(:,2),'o','linewidth',0.5,'markersize',6,'color',[0.5 0.5 0.5]), hold off
    hold on, plot(rewardLocs_scaled(:,1),rewardLocs_scaled(:,2),'o','linewidth',0.5,'markersize',6,'color','k'), hold off
end

axis([1 numSpatialBins(2) 1 numSpatialBins(1)]), axis square, set(gca,'xtick',[],'ytick',[],'ydir','normal','fontsize',14)
set(gcf,'color','w'), box on, axis square, 
hold on

if exist('walls.mat')==2
    load('walls.mat');
    for wall = 1:length(walls)
    x_y = compute_locsToBins([walls(wall).x1 walls(wall).y1],numSpatialBins,x_edges,y_edges);
    walls_scaled(wall).x1 = x_y(1); walls_scaled(wall).y1 = x_y(2);
    x_y = compute_locsToBins([walls(wall).x2 walls(wall).y2],numSpatialBins,x_edges,y_edges);
    walls_scaled(wall).x2 = x_y(1); walls_scaled(wall).y2 = x_y(2);
    plot([walls_scaled(wall).x1,walls_scaled(wall).x2],[walls_scaled(wall).y1,walls_scaled(wall).y2],'LineWidth', 3,'color','k');
    hold on
    end
end
    
    
    
    
% if exist('rewardLocs')
%     if ~isempty(barrierConfig_IDs{sessionNum}) && ~isempty(rewardLocs{sessionNum})
%         rewardLocs_scaled = compute_locsToBins(rewardLocs{sessionNum},numSpatialBins,x_edges,y_edges);
%         [~,barriers_mat] = load_idealizedBarrierLocations(rewardLocs_scaled,barrierConfig_IDs{sessionNum});
%         hold on, plot(barriers_mat(:,1),barriers_mat(:,2),'k','linewidth',2), hold off
%     end
% else
%    hold on, plot_barrierLocs(barrierLocs{sessionNum},1,1,numSpatialBins,2), hold off
% end
% if exist('homeLoc') && ~isempty(homeLoc{sessionNum})
%     homeLoc_scaled = compute_locsToBins(homeLoc{sessionNum},numSpatialBins,x_edges,y_edges);                 
%     hold on, plot(homeLoc_scaled(1),homeLoc_scaled(2),'k.','markersize',11), hold off
% end
% if exist('rewardLocs') && ~isempty(rewardLocs{sessionNum})
%     rewardLocs_scaled = compute_locsToBins(rewardLocs{sessionNum},numSpatialBins,x_edges,y_edges);                 
%     hold on, plot(rewardLocs_scaled(:,1),rewardLocs_scaled(:,2),'ko','linewidth',1,'markersize',4), hold off
% end
% 
% axis([1 numSpatialBins(2) 1 numSpatialBins(1)]), axis square, set(gca,'xtick',[],'ytick',[],'ydir','normal','fontsize',14)
% set(gcf,'color','w'), box on, axis square, 


% if exist('rewardLocs')
%     if ~isempty(barrierConfig_IDs{sessionNum}) && ~isempty(rewardLocs{sessionNum})
%         rewardLocs_scaled = compute_locsToBins(rewardLocs{sessionNum},numSpatialBins,x_edges,y_edges);
%         [~,barriers_mat] = load_idealizedBarrierLocations(rewardLocs_scaled,barrierConfig_IDs{sessionNum});
%         hold on, plot(barriers_mat(:,1),barriers_mat(:,2),'k','linewidth',6), hold off
%     end
% else
%    hold on, plot_barrierLocs(barrierLocs{sessionNum},1,1,numSpatialBins,6), hold off
% end
% if exist('homeLoc') && ~isempty(homeLoc{sessionNum})
%     homeLoc_scaled = compute_locsToBins(homeLoc{sessionNum},numSpatialBins,x_edges,y_edges);                 
%     hold on, plot(homeLoc_scaled(1),homeLoc_scaled(2),'k.','markersize',30), hold off
% end
% if exist('rewardLocs') && ~isempty(rewardLocs{sessionNum})
%     rewardLocs_scaled = compute_locsToBins(rewardLocs{sessionNum},numSpatialBins,x_edges,y_edges);                 
%     hold on, plot(rewardLocs_scaled(:,1),rewardLocs_scaled(:,2),'ko','linewidth',3,'markersize',8), hold off
% end
% 
% axis([1 numSpatialBins(2) 1 numSpatialBins(1)]), axis square, set(gca,'xtick',[],'ytick',[],'ydir','normal','fontsize',14)
% set(gcf,'color','w'), box on, axis square, 
