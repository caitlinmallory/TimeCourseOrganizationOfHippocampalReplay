function [points_of_replay_crossing] = computePointsOfIntersectionWithConcentricCircles(replay_path,ratPos)

% Brad's method draws circles centered at the rat's current location,
% starting with a radius of 15 cm (7.5 bins) and increasing by 2 cm (1 bin).

ring_radius = 1:1:90;


points_of_replay_crossing = nan(size(ring_radius,2),2);
distance_rat_to_replay = sqrt((replay_path(:,1)-ratPos(1)).^2+(replay_path(:,2)-ratPos(2)).^2);
% this is the actual distance between the rat and the replay at all sampling points. We will use this
% to consider when this radius matches the one we are looking for.

if sum(~isnan(distance_rat_to_replay))>0
    for i = 1:length(ring_radius)

        desired_radius = ring_radius(i);

        % Find the first time the replay path crosses the circle
        distances = [distance_rat_to_replay(1:end-1) distance_rat_to_replay(2:end)];

        points_of_crossing = replay_path((distances(:,1) < desired_radius & distances(:,2) > desired_radius) | ...
            (distances(:,1) > desired_radius & distances(:,2)  < desired_radius) | ...
            (distances(:,1) == desired_radius),:);

        if ~isempty(points_of_crossing)
            points_of_replay_crossing(i,:) = points_of_crossing(1,1:2);
        end
    end
end