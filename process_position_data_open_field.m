function [Position_Data, maze_coordinates, new_maze_processing_params] = process_position_data_open_field4(ts, posx, posy, head_direction, time_clip, maze_coordinates, params, maze_processing_params)

% NOTE: In this version the edges of the arena are defined by the user
% selected coordinates. Everything outside of these bounds will be
% eliminated. This is NOT what we want to do on the linear track. There,
% the edges should be the rewards.


use_2d_speed = 1;
posx0 = posx; 
posy0 = posy;

clear posx; clear posy;

if ~isempty(time_clip)
    posx0 = posx0(ts >= time_clip(1) & ts < time_clip(2));
    posy0 = posy0(ts >= time_clip(1) & ts < time_clip(2));
    head_direction = head_direction(ts >= time_clip(1) & ts < time_clip(2));
    ts = ts(ts >= time_clip(1) & ts < time_clip(2));
end

% Interpolate the position data so that its sampled at a constant rate.
% Although the spike data is sampled at 30,000 Hz, the position data is sampled around 20 Hz
% However it's not consistent, so interpolate so that there is position data at 20 Hz (Should increase by 1500 samples)
% linearally interpolate x and y position so that they are sampled at the new timestamps: vq = interp1(x,v,xq
% Interpolate the position data so that its sampled at a constant rate.
% Although the spike data is sampled at 30,000 Hz, the position data is sampled around 20 Hz
% However it's not consistent, so interpolate so that there is position data at 20 Hz (Should increase by 1500 samples)
% linearally interpolate x and y position so that they are sampled at the new timestamps: vq = interp1(x,v,xq)

% Before interpolation, you need to remove duplicate timestamps
[ts_in_clip_clean, ia, ic] = unique(ts,'stable');
posx_in_clip_clean = posx0(ia);
posy_in_clip_clean = posy0(ia);
head_direction_clean = head_direction(ia);


ts_interp = min(ts_in_clip_clean):params.spike_sampling_rate/params.desired_position_sampling_rate:max(ts_in_clip_clean);
posx_in_clip_interp = interp1(ts_in_clip_clean,posx_in_clip_clean,ts_interp);
posy_in_clip_interp = interp1(ts_in_clip_clean,posy_in_clip_clean,ts_interp);
head_direction_in_clip_interp = interp1(ts_in_clip_clean,head_direction_clean,ts_interp);

% Write over the original time and position values:
posx0 = posx_in_clip_interp';
posy0 = posy_in_clip_interp';
head_direction = head_direction_in_clip_interp';
ts = ts_interp';


if isempty(maze_processing_params)
    maze_processing_params.straighten_image = [];
    maze_processing_params.flip_vertically = [];
    maze_processing_params.flip_horizontally = [];
    maze_processing_params.rotate_clockwise = [];
    maze_processing_params.rotate_counterclockwise = [];
    maze_processing_params.boundaries = [];
end



if ~isempty(maze_processing_params)
    new_maze_processing_params = [];
end

figure()
plot(posx0,posy0)

% Inputs: timestamps (ts), x-position (posx) and
% y-position (posy)
% time_clip: the start and stop timestamp for the segment you want to
% analyze
% params: constants

% Outputs: Position_Data is a matrix containing the interpolated timestamps
% (column 1), smoothed x-position (column 2), smoothed y-position (column
% 3), head direction (column 4), smoothed velocity (column 5), x-movement
% direction (column 6) and y-movement direction (column 7)

if maze_processing_params.straighten_image == 1
    position_data = [posx0';posy0'];
    [rotated_position_data,~, ~, maze_coordinates] = rotate_position_data_by_specified_angle(position_data,position_data, position_data, maze_coordinates, maze_processing_params.theta);
    posx0 = rotated_position_data(1,:);
    posy0 = rotated_position_data(2,:);
elseif isempty(maze_processing_params.straighten_image)
    % Ask user if they want to straighten the image
    % Also save the response to "new_maze_processing_params" for future use
    new_maze_processing_params.straighten_image = 0;
    fprintf('Do you want to straighten the image? Click Y or N\n')
    [~,~,response] = ginput(1);

    if response == 89
        position_data = [posx0';posy0'];
        [rotated_position_data, ~,~, ~, ~, ~, ~, maze_coordinates, theta] = rotate_position_data(position_data,position_data, position_data,maze_coordinates);

        posx0 = rotated_position_data(1,:);
        posy0 = rotated_position_data(2,:);
        new_maze_processing_params.straighten_image = 1;
        new_maze_processing_params.theta = theta;
        head_direction = wrapTo2Pi(head_direction+theta);
  
    end
end
clf
plot(posx0,posy0)
%%
% flip vertically if requested
if maze_processing_params.flip_vertically == 1
    posy0 = -1.*(posy0);
    maze_coordinates(:,2) = -1.*(maze_coordinates(:,2));
elseif isempty(maze_processing_params.flip_vertically)
    % Ask user if they want to flip the image vertically:
    % Also save the response to "new_maze_processing_params" for future use
    new_maze_processing_params.flip_vertically = 0;
    fprintf('Do you want to flip the image vertically? Click Y or N\n')

    [~,~,response] = ginput(1);

    if response == 89
        new_maze_processing_params.flip_vertically = 1;
        posy0 = -1.*(posy0);
        maze_coordinates(:,2) = -1.*(maze_coordinates(:,2));
        % head_direction = head_direction+pi; % Not correct - would need to figure out the right transform
    end
end

% plot the new version:
clf
plot(posx0,posy0);

%% flip horizontally if requested
if maze_processing_params.flip_horizontally == 1
    posx0 = -1.*(posx0);
    maze_coordinates(:,1) = -1.*(maze_coordinates(:,1));
elseif isempty(maze_processing_params.flip_horizontally)
    % Ask user if they want to flip the image horizontally:
    % Also save the response to "new_maze_processing_params" for future use
    new_maze_processing_params.flip_horizontally = 0;
    fprintf('Do you want to flip the image horizontally? Click Y or N\n')
    [~,~,response] = ginput(1);

    if response == 89
        new_maze_processing_params.flip_horizontally = 1;
        posx0 = -1.*(posx0);
        maze_coordinates(:,1) = -1.*(maze_coordinates(:,1));
         % head_direction = head_direction+pi; % Not correct - would need to figure out the right transform
    end
end

clf
plot(posx0,posy0);
%% Rotate 90 degrees clockwise if requested
if maze_processing_params.rotate_clockwise == 1
    position_matrix = [posx0'; posy0'];
    maze_coordinate_matrix = [maze_coordinates(:,1)'; maze_coordinates(:,2)'];
    theta = -1*pi/2;
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    rot_position_matrix = R*position_matrix;
    rot_maze_coordinate_matrix = R*maze_coordinate_matrix;
    maze_coordinates(:,1) = rot_maze_coordinate_matrix(1,:)';
    maze_coordinates(:,2) = rot_maze_coordinate_matrix(2,:)';
    posx0 = (rot_position_matrix(1,:))';
    posy0 = (rot_position_matrix(2,:))';
    posx1 = (rot_position_matrix(1,:))';
    posy1 = (rot_position_matrix(2,:))';
    posx2 = (rot_position_matrix(1,:))';
    posy2 = (rot_position_matrix(2,:))';


elseif isempty(maze_processing_params.rotate_clockwise)
    fprintf('Do you want to rotate the image 90 degrees clockwise? Click Y or N\n')
    % Also save the response to "new_maze_processing_params" for future use
    new_maze_processing_params.rotate_clockwise = 0;

    [~,~,response] = ginput(1);
    if response == 89
        new_maze_processing_params.rotate_clockwise = 1;
        position_matrix = [posx0'; posy0'];
        maze_coordinate_matrix = [maze_coordinates(:,1)'; maze_coordinates(:,2)'];
        theta = -1*pi/2;
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        rot_position_matrix = R*position_matrix;
        rot_maze_coordinate_matrix = R*maze_coordinate_matrix;
        maze_coordinates(:,1) = rot_maze_coordinate_matrix(1,:)';
        maze_coordinates(:,2) = rot_maze_coordinate_matrix(2,:)';
        posx0 = (rot_position_matrix(1,:))';
        posy0 = (rot_position_matrix(2,:))';
        posx1 = (rot_position_matrix(1,:))';
        posy1 = (rot_position_matrix(2,:))';
        posx2 = (rot_position_matrix(1,:))';
        posy2 = (rot_position_matrix(2,:))';
    end
end
clf
plot(posx0,posy0);
%% Rotate 90 degrees counterclockwise if requested
if maze_processing_params.rotate_counterclockwise == 1
    position_matrix = [posx0'; posy0'];
    maze_coordinate_matrix = [maze_coordinates(:,1)'; maze_coordinates(:,2)'];
    theta = pi/2;
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    rot_position_matrix = R*position_matrix;
    rot_maze_coordinate_matrix = R*maze_coordinate_matrix;
    maze_coordinates(:,1) = rot_maze_coordinate_matrix(1,:)';
    maze_coordinates(:,2) = rot_maze_coordinate_matrix(2,:)';
    posx0 = (rot_position_matrix(1,:))';
    posy0 = (rot_position_matrix(2,:))';
    posx1 = (rot_position_matrix(1,:))';
    posy1 = (rot_position_matrix(2,:))';
    posx2 = (rot_position_matrix(1,:))';
    posy2 = (rot_position_matrix(2,:))';
elseif isempty(maze_processing_params.rotate_counterclockwise)
    %Ask user if they want to rotate the image 90 degrees counterclockwise
    % Also save the response to "new_maze_processing_params" for future use
    new_maze_processing_params.rotate_counterclockwise = 0;
    fprintf('Do you want to rotate the image 90 degrees counterclockwise? Click Y or N\n')

    [~,~,response] = ginput(1);
    if response == 89
        new_maze_processing_params.rotate_counterclockwise = 1;
        position_matrix = [posx0'; posy0'];
        maze_coordinate_matrix = [maze_coordinates(:,1)'; maze_coordinates(:,2)'];
        theta = pi/2;
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        rot_position_matrix = R*position_matrix;
        rot_maze_coordinate_matrix = R*maze_coordinate_matrix;
        maze_coordinates(:,1) = rot_maze_coordinate_matrix(1,:)';
        maze_coordinates(:,2) = rot_maze_coordinate_matrix(2,:)';
        posx0 = (rot_position_matrix(1,:))';
        posy0 = (rot_position_matrix(2,:))';
        posx1 = (rot_position_matrix(1,:))';
        posy1 = (rot_position_matrix(2,:))';
        posx2 = (rot_position_matrix(1,:))';
        posy2 = (rot_position_matrix(2,:))';
    end
end

posx0 = posx0';
posy0 = posy0';
plot(posx0,posy0);

% % Calculate the current head direction from the 2 leds:
% head_direction = (calculate_headDirection(posx1,posy1,posx2,posy2));

%%
if isempty(maze_processing_params.boundaries)
    % Re-select the boundaries:
    fprintf('Click on the new left boundary of the maze\n')
    coordinates_end_1 = zoomginput(1);
    fprintf('Click on the new right boundary of the maze\n')
    coordinates_end_2 = zoomginput(1);
    fprintf('Click on the new bottom of the maze\n')
    coordinates_bottom = zoomginput(1);
    fprintf('Click on the new top of the maze\n')
    coordinates_top = zoomginput(1);

    new_maze_processing_params.boundaries.left = coordinates_end_1;
    new_maze_processing_params.boundaries.right = coordinates_end_2;
    new_maze_processing_params.boundaries.bottom = coordinates_bottom;
    new_maze_processing_params.boundaries.top = coordinates_top;

else
    coordinates_end_1 = maze_processing_params.boundaries.left;
    coordinates_end_2 = maze_processing_params.boundaries.right;
    coordinates_top = maze_processing_params.boundaries.top;
    coordinates_bottom = maze_processing_params.boundaries.bottom;
end

% smooth:
if params.smooth_position_data

    % BRAD's METHOD: seems to be over smoothing.
    %         [filterA,filterB]=butter(2,0.025);
%     [filterA,filterB]=butter(2,0.05); % Caitlin likes this better.
%     posx0=filtfilt(filterA,filterB,posx0);
%     posy0=filtfilt(filterA,filterB,posy0);

    % ARCHT's METHOD: looks much closer to the real data.
    %     filter = ones(3,1);
    %     filter = filter/length(filter);
    %     posx = conv(posx,filter,'same');
    %     posy= conv(posy,filter,'same');

    % Caitlin's Method: Used for deeplabcut-tracked position as well.
         filter = ones(13,1)/13;
         posx0 = filtfilt(filter,1,posx0);
         posy0 = filtfilt(filter,1,posy0);
   % remove edge effects:
        posx0(1) = posx0(2);
        posx0(end) = posx0(end-1);
        posy0(1)= posy0(2);
        posy0(end) = posy0(end-1);
end

% remove points outside the boundaries
posx0(posx0<coordinates_end_1(1)) = deal(coordinates_end_1(1));
posx0(posx0>coordinates_end_2(1)) = deal(coordinates_end_2(1));
posy0(posy0<coordinates_bottom(2)) = deal(coordinates_bottom(2));
posy0(posy0>coordinates_top(2)) = deal(coordinates_top(2));
clf;
plot(posx0,posy0)

%scale the position data by the arena size
% maze_length_x = max(abs(coordinates_end_2(1)-coordinates_end_1(1)));
% maze_length_y = max(abs(coordinates_top(2)-coordinates_bottom(2)));
maze_length_x = max(posx0)-min(posx0);
maze_length_y = max(posy0)-min(posy0);
scaling_factor_x = params.maze_size/maze_length_x;

if params.spatialDim == 1
    disp('Linear track session! Only using X-dim for scaling')
    scaling_factor_y = scaling_factor_x;
else
    scaling_factor_y = params.maze_size/maze_length_y;
end
posx0 = posx0.*scaling_factor_x;
posy0 = posy0.*scaling_factor_y;


maze_coordinates(:,1) = maze_coordinates(:,1).*scaling_factor_x;
maze_coordinates(:,2) = maze_coordinates(:,2).*scaling_factor_y;

% min_x = min([coordinates_end_1(1),coordinates_end_2(1)].*scaling_factor_x);
% min_y = min([coordinates_top(2),coordinates_bottom(2)].*scaling_factor_y);
min_x = min(posx0); 
min_y = min(posy0);

posx0 = posx0-min_x+0.001;
posy0 = posy0-min_y+0.001;
maze_coordinates(:,1) = maze_coordinates(:,1)-min_x+0.001;
maze_coordinates(:,2) = maze_coordinates(:,2)-min_y+0.001;

close all

posx = posx0';
posy = posy0';

% determine speed of the animal at each time point (Euclidean difference in
% position/time between position points, where time between position points
% is equal to the interpolated sampling rate (20 Hz)).
if use_2d_speed == 1
    speed = sqrt((diff(posx)).^2 + (diff(posy)).^2)/(1/params.desired_position_sampling_rate);
else
    speed = sqrt(diff(posx).^2)/(1/params.desired_position_sampling_rate);
end
speed = [speed; speed(end)];

% remove stray points where the speed is way too high to be real.
bad_inds = find(speed >= 150);

all_bad_inds = [];
for i = 1:length(bad_inds)
    all_bad_inds = [all_bad_inds bad_inds(i)-3:bad_inds(i)+3];
end
all_bad_inds = unique(all_bad_inds);

good_inds = setdiff(1:length(speed),all_bad_inds);


posx_interp = (interp1(good_inds,posx(good_inds),1:length(posx)))';
posy_interp = (interp1(good_inds,posy(good_inds),1:length(posy)))';
speed_interp = (interp1(good_inds,speed(good_inds),1:length(speed)))';

if ~isempty(find(isnan(posx_interp)))
    if abs(1-find(isnan(posx_interp),1,'last')) < abs(length(posx_interp)-find(isnan(posx_interp),1,'last'))
        posx_interp(isnan(posx_interp)) = posx_interp(find(~isnan(posx_interp),1,'first'));
        posy_interp(isnan(posy_interp)) = posy_interp(find(~isnan(posy_interp),1,'first'));
        speed_interp(isnan(speed_interp)) = speed_interp(find(~isnan(speed_interp),1,'first'));
    elseif abs(1-find(isnan(posx_interp),1,'last')) > abs(length(posx_interp)-find(isnan(posx_interp),1,'last'))
        posx_interp(isnan(posx_interp)) = posx_interp(find(~isnan(posx_interp),1,'last'));
        posy_interp(isnan(posy_interp)) = posy_interp(find(~isnan(posy_interp),1,'last'));
        speed_interp(isnan(speed_interp)) = speed_interp(find(~isnan(speed_interp),1,'last'));
    end
end

posx = posx_interp;
posy = posy_interp;
speed = speed_interp;

smoothing_weights = ones(params.speed_smoothing_window_length,1)/params.speed_smoothing_window_length;
speed_smoothed = filtfilt(smoothing_weights,1,speed);

figure(99)
plot(ts,speed_smoothed);


% For consistnecy with Brad's code:
% create an array called Position Data
% containing each position timestamp and the associated xpos, ypos, hd
% (nan's for now), and speed
% create an array called Spike_Information containing each spike timestamp
% and the corresponding position data for that timestamp.
x_movement = diff(posx);
x_movement = [x_movement; x_movement(end)];
y_movement = diff(posy);
y_movement = [y_movement; y_movement(end)];


% Sanity check:
figure()
subplot(1,3,1);
hist(head_direction)
subplot(1,3,2)
hist(head_direction(posx>params.maze_size-15));
title('HD when rat is at right end');
subplot(1,3,3)
hist(head_direction(posx<15));
title('HD when rat is at left end');

% calculate angular velocity
change_in_hd = circ_dist(head_direction(1:end-1),head_direction(2:end));
% not currently being used


Position_Data = [ts posx posy head_direction speed_smoothed x_movement y_movement];

figure()
plot(posx,posy); hold on;
plot(maze_coordinates(:,1),maze_coordinates(:,2),'ok')