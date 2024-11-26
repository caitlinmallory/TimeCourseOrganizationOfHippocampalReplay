function [angles, points_of_path1_crossing, points_of_path2_crossing] = computeAngleBetweenTwoPaths(path1, path2, ratPos)
% computeAngleBetweenTwoPaths computes the angle between two paths at varying distances from the rat's current position.
%
% Inputs:
%   - path1: An N x 2 matrix representing the first path (each row is a point in 2D space).
%   - path2: An N x 2 matrix representing the second path (each row is a point in 2D space).
%   - ratPos: A 1 x 2 vector representing the current position of the rat in 2D space [x, y].
%
% Outputs:
%   - angles: A column vector containing the angles between the two paths at various radii from the rat.
%   - points_of_path1_crossing: A matrix of crossing points on path1 at each radius from the rat.
%   - points_of_path2_crossing: A matrix of crossing points on path2 at each radius from the rat.

% Define the range of ring radii, from 1 to 90, increasing by 1 unit at each step.
ring_radius = 1:1:90;

% Initialize output variables.
angles = nan(size(ring_radius, 2), 1);  % Preallocate the angles vector to store the computed angles.
points_of_path1_crossing = nan(size(ring_radius, 2), 2);  % Preallocate for crossing points of path1.
points_of_path2_crossing = nan(size(ring_radius, 2), 2);  % Preallocate for crossing points of path2.

% Calculate the distance from the rat's position to each point on path1 and path2.
distance_rat_to_path1 = sqrt((path1(:, 1) - ratPos(1)).^2 + (path1(:, 2) - ratPos(2)).^2);
distance_rat_to_path2 = sqrt((path2(:, 1) - ratPos(1)).^2 + (path2(:, 2) - ratPos(2)).^2);

% Loop through each desired radius to find the crossing points and calculate angles.
for i = 1:length(ring_radius)
    % Set the desired radius for the current loop iteration.
    desired_radius = ring_radius(i);

    % --- Find the first crossing point on path1 ---
    % Compute the distances between consecutive points along path1.
    distances = [distance_rat_to_path1(1:end-1), distance_rat_to_path1(2:end)];

    % Find points on path1 where the distance crosses the desired radius (either entering or leaving the circle).
    points_of_crossing = path1((distances(:, 1) < desired_radius & distances(:, 2) > desired_radius) | ...
        (distances(:, 1) > desired_radius & distances(:, 2) < desired_radius) | ...
        (distances(:, 1) == desired_radius), :);

    % If a crossing point is found, store it in points_of_path1_crossing.
    if ~isempty(points_of_crossing)
        points_of_path1_crossing(i, :) = points_of_crossing(1, 1:2);  % Store the first crossing point.
    end

    % --- Find the first crossing point on path2 ---
    % Compute the distances between consecutive points along path2.
    distances = [distance_rat_to_path2(1:end-1), distance_rat_to_path2(2:end)];

    % Find points on path2 where the distance crosses the desired radius.
    points_of_crossing = path2((distances(:, 1) < desired_radius & distances(:, 2) > desired_radius) | ...
        (distances(:, 1) > desired_radius & distances(:, 2) < desired_radius) | ...
        (distances(:, 1) == desired_radius), :);

    % If a crossing point is found, store it in points_of_path2_crossing.
    if ~isempty(points_of_crossing)
        points_of_path2_crossing(i, :) = points_of_crossing(1, 1:2);  % Store the first crossing point.
    end

    % --- Calculate the angle between the two vectors ---
    % Compute the vectors u and v from the rat's position to the crossing points on path1 and path2, respectively.
    u = [points_of_path1_crossing(i, 1) - ratPos(1), points_of_path1_crossing(i, 2) - ratPos(2)];
    v = [points_of_path2_crossing(i, 1) - ratPos(1), points_of_path2_crossing(i, 2) - ratPos(2)];

    % Calculate the angle between the two vectors using the 2D cross product and dot product.
    % The atan2 function computes the signed angle between the vectors u and v.
    angles(i) = rad2deg(atan2(u(1) * v(2) - u(2) * v(1), u(1) * v(1) + u(2) * v(2)));
end

end
