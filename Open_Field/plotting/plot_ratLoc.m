% hold on, plot(ratPos(1),ratPos(2),'k*','markersize',10,'linewidth',1), hold off
% hold on, plot(ratPos(1),ratPos(2),'ko','markersize',12,'linewidth',2), hold off
% hold on, quiver(ratPos(1),ratPos(2),cos(ratHD),sin(ratHD),3,'k','linewidth',3), hold off


% if exist('replay_laser_state') &&  sum(replay_laser_state) >= 1
%     hold on, plot(ratPos(1),ratPos(2),'r*','markersize',6,'linewidth',1), hold off
%     hold on, plot(ratPos(1),ratPos(2),'ro','markersize',8,'linewidth',2), hold off
%     hold on, quiver(ratPos(1),ratPos(2),cos(ratHD),sin(ratHD),5,'r','linewidth',2), hold off
% else
%     hold on, plot(ratPos(1),ratPos(2),'k*','markersize',6,'linewidth',1), hold off
%     hold on, plot(ratPos(1),ratPos(2),'ko','markersize',8,'linewidth',2), hold off
%     hold on, quiver(ratPos(1),ratPos(2),cos(ratHD),sin(ratHD),5,'k','linewidth',2), hold off
% end

if exist('replay_laser_state') &&  sum(replay_laser_state) >= 1
    hold on, plot(ratPos(1),ratPos(2),'r*','markersize',5,'linewidth',1), hold off
    hold on, plot(ratPos(1),ratPos(2),'ro','markersize',6,'linewidth',1), hold off
    hold on, quiver(ratPos(1),ratPos(2),cos(ratHD),sin(ratHD),5,'r','linewidth',2), hold off
else
    hold on, plot(ratPos(1),ratPos(2),'k*','markersize',5,'linewidth',1), hold off
    hold on, plot(ratPos(1),ratPos(2),'ko','markersize',6,'linewidth',1), hold off
    hold on, quiver(ratPos(1),ratPos(2),cos(ratHD),sin(ratHD),5,'k','linewidth',1), hold off
end