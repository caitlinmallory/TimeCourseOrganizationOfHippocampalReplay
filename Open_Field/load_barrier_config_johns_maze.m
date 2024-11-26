well_well_distance = 25;
well_wall_distance = (90 - 2*well_well_distance)/2;

wells = [well_wall_distance well_wall_distance; 45 well_wall_distance; 90-well_wall_distance well_wall_distance;...
    well_wall_distance 45; 45 45; 90-well_wall_distance 45;
    well_wall_distance  90-well_wall_distance; 45 90-well_wall_distance; 90-well_wall_distance 90-well_wall_distance];

% for i = 1:length(wells)
% well = circle(wells(i,1),wells(i,2),2);
% plot(well(:,1),well(:,2),'k','LineWidth',2)
% hold on
% end

half_barrier = 10;
% x1 x2 y1 y2
% Horizontal Barriers:
barrier1 = [wells(1,1)-half_barrier, wells(1,1)+half_barrier+2.5, wells(1,2)+well_well_distance/2 ,wells(1,2)+well_well_distance/2];
barrier2 = [wells(2,1)-half_barrier-2.5, wells(2,1)+half_barrier+2.5, wells(2,2)+well_well_distance/2 ,wells(2,2)+well_well_distance/2];
barrier3 = [wells(3,1)-half_barrier-2.5, wells(3,1)+half_barrier, wells(3,2)+well_well_distance/2 ,wells(3,2)+well_well_distance/2];
barrier4 = [wells(4,1)-half_barrier, wells(4,1)+half_barrier+2.5, wells(4,2)+well_well_distance/2 ,wells(4,2)+well_well_distance/2];
barrier5 = [wells(5,1)-half_barrier-2.5, wells(5,1)+half_barrier+2.5, wells(5,2)+well_well_distance/2 ,wells(5,2)+well_well_distance/2];
barrier6 = [wells(6,1)-half_barrier-2.5, wells(6,1)+half_barrier, wells(6,2)+well_well_distance/2 ,wells(6,2)+well_well_distance/2];
% Vertical Barriers;
barrier7 = [wells(1,1)+well_well_distance/2 wells(1,1)+well_well_distance/2 wells(1,2)-half_barrier wells(1,2)+half_barrier+2.5];
barrier8 = [wells(4,1)+well_well_distance/2 wells(4,1)+well_well_distance/2 wells(4,2)-half_barrier-2.5 wells(4,2)+half_barrier+2.5];
barrier9 = [wells(7,1)+well_well_distance/2 wells(7,1)+well_well_distance/2 wells(7,2)-half_barrier-2.5 wells(7,2)+half_barrier];
barrier10 = [wells(2,1)+well_well_distance/2 wells(2,1)+well_well_distance/2 wells(2,2)-half_barrier wells(2,2)+half_barrier+2.5];
barrier11 = [wells(5,1)+well_well_distance/2 wells(5,1)+well_well_distance/2 wells(5,2)-half_barrier-2.5 wells(5,2)+half_barrier+2.5];
barrier12 = [wells(8,1)+well_well_distance/2 wells(8,1)+well_well_distance/2 wells(8,2)-half_barrier-2.5 wells(8,2)+half_barrier];

if rat == 12 || rat == 13 || rat == 14
    wells = wells([7 8 9 4 5 6 1 2 3],:);
end


barriers = [barrier1; barrier2; barrier3; barrier4; barrier5; barrier6; barrier7; barrier8; barrier9; barrier10; barrier11; barrier12];
% hold on
for i = 1:length(barriers)
line([barriers(i,1) barriers(i,2)], [barriers(i,3) barriers(i,4)],'LineWidth',4,'Color','k'); hold on;
end


