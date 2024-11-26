
% mark the location of all wells with open circles
for i = 1:length(wells)
well = circle(wells(i,1),wells(i,2),2);
plot(well(:,1),well(:,2),'k','LineWidth',2)
hold on
end

% mark the location of the home well with a filled circle
plot(wells(home_well,1),wells(home_well,2),'.k','markerSize',40)

% draw in the barriers for this session
if sum(Experiment_Information.barriers_session(segment,:) > 0)
    for i = Experiment_Information.barriers_session(segment,:)
        line([barriers(i,1) barriers(i,2)], [barriers(i,3) barriers(i,4)],'LineWidth',6,'Color','k'); hold on;
    end
end