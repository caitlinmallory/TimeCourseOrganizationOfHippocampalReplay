function [Position_Data, maze_coordinates, new_maze_processing_params] = process_position_data_open_field3(ts, posx, posy, posx1, posy1, posx2, posy2, time_clip, maze_coordinates, params, maze_processing_params, laser_state)

align_laser_on_laser_off_positions = 0;
use_2d_speed = 1;
led_choice = 3; % 0 = use the 'middle point' detected by deeplabcut; 1 = use green led; 2 = use red led; 3 = use midpoint of green and red

posx0 = posx; 
posy0 = posy;

clear posx; clear posy;

if ~isempty(time_clip)
    posx0 = posx0(ts >= time_clip(1) & ts < time_clip(2));
    posy0 = posy0(ts >= time_clip(1) & ts < time_clip(2));
    posx1 = posx1(ts >= time_clip(1) & ts < time_clip(2));
    posy1 = posy1(ts >= time_clip(1) & ts < time_clip(2));
    posx2 = posx2(ts >= time_clip(1) & ts < time_clip(2));
    posy2 = posy2(ts >= time_clip(1) & ts < time_clip(2));
    ts = ts(ts >= time_clip(1) & ts < time_clip(2));
end

% Interpolate the position data so that its sampled at a constant rate.
% Although the spike data is sampled at 30,000 Hz, the position data is sampled around 20 Hz
% However it's not consistent, so interpolate so that there is position data at 20 Hz (Should increase by 1500 samples)
% linearally interpolate x and y position so that they are sampled at the new timestamps: vq = interp1(x,v,xq)

% Before interpolation, you need to remove duplicate timestamps
[ts_no_duplicates, ia, ic] = unique(ts,'stable');
posx0_no_duplicates = posx0(ia);
posy0_no_duplicates = posy0(ia);
posx1_no_duplicates = posx1(ia);
posy1_no_duplicates = posy1(ia);
posx2_no_duplicates = posx2(ia);
posy2_no_duplicates = posy2(ia);

ts = (min(ts_no_duplicates):params.spike_sampling_rate/params.desired_position_sampling_rate:max(ts_no_duplicates))';
posx0 = interp1(ts_no_duplicates,posx0_no_duplicates,ts);
posy0 = interp1(ts_no_duplicates,posy0_no_duplicates,ts);
posx1 = interp1(ts_no_duplicates,posx1_no_duplicates,ts);
posy1 = interp1(ts_no_duplicates,posy1_no_duplicates,ts);
posx2 = interp1(ts_no_duplicates,posx2_no_duplicates,ts);
posy2 = interp1(ts_no_duplicates,posy2_no_duplicates,ts);



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



% Inputs: timestamps (ts), x-position (posx) and
% y-position (posy)
% time_clip: the start and stop timestamp for the segment you want to
% analyze
% params: constants

% Outputs: Position_Data is a matrix containing the interpolated timestamps
% (column 1), smoothed x-position (column 2), smoothed y-position (column
% 3), head direction (column 4), smoothed velocity (column 5), x-movement
% direction (column 6) and y-movement direction (column 7)


% ToDo: figure out how to remove any wildly crazy points:
diff_in_position = sqrt( (diff(posx0)).^2 + (diff(posy0)).^2 );
bad_points = find([0; diff_in_position] > 200);
bad_points = [bad_points; bad_points-1];
figure()
plot(posx0,posy0);
hold on
plot(maze_coordinates(:,1),maze_coordinates(:,2),'.k','MarkerSize',20)
hold on;
plot(maze_coordinates(1,1),maze_coordinates(1,2),'.g','MarkerSize',20);
hold on
plot(maze_coordinates(36,1),maze_coordinates(36,2),'.r','MarkerSize',20);
hold on
plot(posx0(bad_points),posy0(bad_points),'.m','MarkerSize',20);
% posx_in_clip(bad_points) = nan;
% posy_in_clip(bad_points) = nan;
% bad_points = sort(bad_points);
% diff_bad_points = [1; diff(bad_points)];

if maze_processing_params.straighten_image == 1

    position_data = [posx0';posy0'];
    position_data1 = [posx1';posy1'];
    position_data2 = [posx2';posy2'];
    [rotated_position_data,rotated_position_data1, rotated_position_data2, maze_coordinates] = rotate_position_data_by_specified_angle(position_data,position_data1, position_data2, maze_coordinates, maze_processing_params.theta);
    posx0 = rotated_position_data(1,:);
    posy0 = rotated_position_data(2,:);
    posx1 = rotated_position_data1(1,:);
    posy1 = rotated_position_data1(2,:);
    posx2 = rotated_position_data2(1,:);
    posy2 = rotated_position_data2(2,:);

elseif isempty(maze_processing_params.straighten_image)
    % Ask user if they want to straighten the image
    % Also save the response to "new_maze_processing_params" for future use
    new_maze_processing_params.straighten_image = 0;

    fprintf('Do you want to straighten the image? Click Y or N\n')
    [~,~,response] = ginput(1);

    if response == 89
        position_data = [posx0';posy0'];
        position_data1 = [posx1';posy1'];
        position_data2 = [posx2';posy2'];

        [rotated_position_data, rotated_position_data1, rotated_position_data2, ~, ~, ~, ~, maze_coordinates, theta] = rotate_position_data(position_data,position_data1, position_data2,maze_coordinates);

        posx0 = rotated_position_data(1,:);
        posy0 = rotated_position_data(2,:);
        posx1 = rotated_position_data1(1,:);
        posy1 = rotated_position_data1(2,:);
        posx2 = rotated_position_data2(1,:);
        posy2 = rotated_position_data2(2,:);

        new_maze_processing_params.straighten_image = 1;
        new_maze_processing_params.theta = theta;
    end
end
clf
plot(posx0,posy0)
hold on
plot(posx1,posy1,'.g','MarkerSize',1);
hold on
plot(posx2,posy2,'.r','MarkerSize',1);
hold on
plot(maze_coordinates(:,1),maze_coordinates(:,2),'.k','MarkerSize',10);
hold on
plot(maze_coordinates(1,1),maze_coordinates(1,2),'.g','MarkerSize',10);
hold on
plot(maze_coordinates(36,1),maze_coordinates(36,2),'.r','MarkerSize',10);
%%
% flip vertically if requested
if maze_processing_params.flip_vertically == 1
    posy0 = -1.*(posy0);
    posy1 = -1.*(posy1);
    posy2 = -1.*(posy2);
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
        posy1 = -1.*(posy1);
        posy2 = -1.*(posy2);
        maze_coordinates(:,2) = -1.*(maze_coordinates(:,2));
    end
end

% plot the new version:
clf
plot(posx0,posy0);
hold on
plot(posx1,posy1,'.g','MarkerSize',1);
hold on
plot(posx2,posy2,'.r','MarkerSize',1);
hold on
plot(maze_coordinates(:,1),maze_coordinates(:,2),'.k','MarkerSize',10);
hold on
plot(maze_coordinates(1,1),maze_coordinates(1,2),'.g','MarkerSize',10);
hold on
plot(maze_coordinates(36,1),maze_coordinates(36,2),'.r','MarkerSize',10);

%% flip horizontally if requested
if maze_processing_params.flip_horizontally == 1
    posx0 = -1.*(posx0);
    posx1 = -1.*(posx1);
    posx2 = -1.*(posx2);
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
        posx1 = -1.*(posx1);
        posx2 = -1.*(posx2);
        maze_coordinates(:,1) = -1.*(maze_coordinates(:,1));
    end
end

clf
plot(posx0,posy0);
hold on
plot(posx1,posy1,'.g','MarkerSize',1);
hold on
plot(posx2,posy2,'.r','MarkerSize',1);
hold on
plot(maze_coordinates(:,1),maze_coordinates(:,2),'.k','MarkerSize',20);
hold on
plot(maze_coordinates(1,1),maze_coordinates(1,2),'.g','MarkerSize',20);
hold on
plot(maze_coordinates(36,1),maze_coordinates(36,2),'.r','MarkerSize',20);

%% Rotate 90 degrees clockwise if requested
if maze_processing_params.rotate_clockwise == 1
    position_matrix = [posx0'; posy0'; posx1'; posy1'; posx2'; posy2'];
    maze_coordinate_matrix = [maze_coordinates(:,1)'; maze_coordinates(:,2)'];
    theta = -1*pi/2;
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    rot_position_matrix = R*position_matrix;
    rot_maze_coordinate_matrix = R*maze_coordinate_matrix;
    maze_coordinates(:,1) = rot_maze_coordinate_matrix(1,:)';
    maze_coordinates(:,2) = rot_maze_coordinate_matrix(2,:)';
    posx0 = (rot_position_matrix(1,:))';
    posy0 = (rot_position_matrix(2,:))';
    posx1 = (rot_position_matrix(3,:))';
    posy1 = (rot_position_matrix(4,:))';
    posx2 = (rot_position_matrix(5,:))';
    posy2 = (rot_position_matrix(6,:))';


elseif isempty(maze_processing_params.rotate_clockwise)
    fprintf('Do you want to rotate the image 90 degrees clockwise? Click Y or N\n')
    % Also save the response to "new_maze_processing_params" for future use
    new_maze_processing_params.rotate_clockwise = 0;

    [~,~,response] = ginput(1);
    if response == 89
        new_maze_processing_params.rotate_clockwise = 1;
        position_matrix = [posx0'; posy0'; posx1'; posy1'; posx2'; posy2'];
        maze_coordinate_matrix = [maze_coordinates(:,1)'; maze_coordinates(:,2)'];
        theta = -1*pi/2;
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        rot_position_matrix = R*position_matrix;
        rot_maze_coordinate_matrix = R*maze_coordinate_matrix;
        maze_coordinates(:,1) = rot_maze_coordinate_matrix(1,:)';
        maze_coordinates(:,2) = rot_maze_coordinate_matrix(2,:)';
        posx0 = (rot_position_matrix(1,:))';
        posy0 = (rot_position_matrix(2,:))';
        posx1 = (rot_position_matrix(3,:))';
        posy1 = (rot_position_matrix(4,:))';
        posx2 = (rot_position_matrix(5,:))';
        posy2 = (rot_position_matrix(6,:))';
    end
end
clf
plot(posx0,posy0);
hold on
plot(posx1,posy1,'.g','MarkerSize',1);
hold on
plot(posx2,posy2,'.r','MarkerSize',1);
hold on
plot(maze_coordinates(:,1),maze_coordinates(:,2),'.k','MarkerSize',20);
hold on
plot(maze_coordinates(1,1),maze_coordinates(1,2),'.g','MarkerSize',20);
hold on
plot(maze_coordinates(36,1),maze_coordinates(36,2),'.r','MarkerSize',20);
%% Rotate 90 degrees counterclockwise if requested
if maze_processing_params.rotate_counterclockwise == 1
    position_matrix = [posx0'; posy0'; posx1'; posy1'; posx2'; posy2'];
    maze_coordinate_matrix = [maze_coordinates(:,1)'; maze_coordinates(:,2)'];
    theta = pi/2;
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    rot_position_matrix = R*position_matrix;
    rot_maze_coordinate_matrix = R*maze_coordinate_matrix;
    maze_coordinates(:,1) = rot_maze_coordinate_matrix(1,:)';
    maze_coordinates(:,2) = rot_maze_coordinate_matrix(2,:)';
    posx0 = (rot_position_matrix(1,:))';
    posy0 = (rot_position_matrix(2,:))';
    posx1 = (rot_position_matrix(3,:))';
    posy1 = (rot_position_matrix(4,:))';
    posx2 = (rot_position_matrix(5,:))';
    posy2 = (rot_position_matrix(6,:))';
elseif isempty(maze_processing_params.rotate_counterclockwise)
    %Ask user if they want to rotate the image 90 degrees counterclockwise
    % Also save the response to "new_maze_processing_params" for future use
    new_maze_processing_params.rotate_counterclockwise = 0;
    fprintf('Do you want to rotate the image 90 degrees counterclockwise? Click Y or N\n')

    [~,~,response] = ginput(1);
    if response == 89
        new_maze_processing_params.rotate_counterclockwise = 1;
        position_matrix = [posx0'; posy0'; posx1'; posy1'; posx2'; posy2'];
        maze_coordinate_matrix = [maze_coordinates(:,1)'; maze_coordinates(:,2)'];
        theta = pi/2;
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        rot_position_matrix = R*position_matrix;
        rot_maze_coordinate_matrix = R*maze_coordinate_matrix;
        maze_coordinates(:,1) = rot_maze_coordinate_matrix(1,:)';
        maze_coordinates(:,2) = rot_maze_coordinate_matrix(2,:)';
        posx0 = (rot_position_matrix(1,:))';
        posy0 = (rot_position_matrix(2,:))';
        posx1 = (rot_position_matrix(3,:))';
        posy1 = (rot_position_matrix(4,:))';
        posx2 = (rot_position_matrix(5,:))';
        posy2 = (rot_position_matrix(6,:))';
    end
end

posx0 = posx0';
posy0 = posy0';
posx1 = posx1';
posy1 = posy1';
posx2 = posx2';
posy2 = posy2';

plot(posx0,posy0);
hold on
plot(posx1,posy1,'.g','MarkerSize',1);
hold on
plot(posx2,posy2,'.r','MarkerSize',1);
hold on
hold on
plot(maze_coordinates(:,1),maze_coordinates(:,2),'.k','MarkerSize',20);
hold on
plot(maze_coordinates(1,1),maze_coordinates(1,2),'.g','MarkerSize',20);
hold on
plot(maze_coordinates(36,1),maze_coordinates(36,2),'.r','MarkerSize',20);

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
    [filterA,filterB]=butter(2,0.05); % Caitlin likes this better.
    posx0=filtfilt(filterA,filterB,posx0);
    posy0=filtfilt(filterA,filterB,posy0);
    posx1=filtfilt(filterA,filterB,posx1);
    posy1=filtfilt(filterA,filterB,posy1);
    posx2=filtfilt(filterA,filterB,posx2);
    posy2=filtfilt(filterA,filterB,posy2);
    % ARCHT's METHOD: looks much closer to the real data.
    %     filter = ones(3,1);
    %     filter = filter/length(filter);
    %     posx = conv(posx,filter,'same');
    %     posy= conv(posy,filter,'same');

    % Caitlin's Method: Used for deeplabcut-tracked position as well.
    %      filter = ones(13,1)/13;
    %      posx = filtfilt(filter,1,posx);
    %      posy = filtfilt(filter,1,posy);
    %remove edge effects:
    %     posx(1) = posx(2);
    %     posx(end) = posx(end-1);
    %     posy(1)= posy(2);
    %     posy(end) = posy(end-1);
end

%scale the position data by the arena size
maze_length_x = max(abs(coordinates_end_2(1)-coordinates_end_1(1)));
maze_length_y = max(abs(coordinates_top(2)-coordinates_bottom(2)));
scaling_factor_x = params.maze_size/maze_length_x;

if params.spatialDim == 1
    disp('Linear track session! Only using X-dim for scaling')
    scaling_factor_y = scaling_factor_x;
else
    scaling_factor_y = params.maze_size/maze_length_y;
end
posx0 = posx0.*scaling_factor_x;
posy0 = posy0.*scaling_factor_y;
posx1 = posx1.*scaling_factor_x;
posy1 = posy1.*scaling_factor_y;
posx2 = posx2.*scaling_factor_x;
posy2 = posy2.*scaling_factor_y;

maze_coordinates(:,1) = maze_coordinates(:,1).*scaling_factor_x;
maze_coordinates(:,2) = maze_coordinates(:,2).*scaling_factor_y;


min_x = min([coordinates_end_1(1),coordinates_end_2(1)].*scaling_factor_x);
min_y = min([coordinates_top(2),coordinates_bottom(2)].*scaling_factor_y);

posx0 = posx0-min_x+0.001;
posy0 = posy0-min_y+0.001;
posx1 = posx1-min_x + 0.001;
posy1 = posy1- min_y+0.001;
posx2 = posx2-min_x + 0.001;
posy2 = posy2- min_y+0.001;

maze_coordinates(:,1) = maze_coordinates(:,1) - min_x + 0.001;
maze_coordinates(:,2) = maze_coordinates(:,2) - min_y + 0.001;

figure()
plot(posx1,posy1);
hold on
plot(maze_coordinates(:,1),maze_coordinates(:,2),'.k','MarkerSize',20);
hold on
plot(maze_coordinates(1,1),maze_coordinates(1,2),'.g','MarkerSize',20);
hold on
plot(maze_coordinates(36,1),maze_coordinates(36,2),'.r','MarkerSize',20);


laser_state = compute_dataTemporalConcatenation(laser_state,[ts(1) ts(end)]);
laser_state = compute_dataInterpolation(laser_state,ts,[]);
laser_state = laser_state(:,2);

% Now calculate the head direction from the 2 leds:
head_direction = calculate_headDirection(posx1,posy1,posx2,posy2);

close all
% ax1 = subplot(2,2,1);
% plot(posx1(laser_state == 0),posy1(laser_state == 0),'.g');
% hold on
% plot(posx1(laser_state == 1),posy1(laser_state == 1),'.b');
% ax2 = subplot(2,2,2);
% plot(posx2(laser_state == 0),posy2(laser_state == 0),'.r');
% hold on
% plot(posx2(laser_state == 1),posy2(laser_state == 1),'.m');

save('position_tracking_data','laser_state','ts','posx0','posy0','posx1','posy1','posx2','posy2')


if align_laser_on_laser_off_positions == 1

align_deeplabcut_tracking_laser_off_on_v3


% load position_tracking_data_BlockSession_transformed_by_Dave.mat
% posx1_dave = Expression1(:,2);
% posy1_dave = Expression1(:,3);
% posx2_dave = Expression1(:,4);
% posy2_dave = Expression1(:,5);
% posx_dave = (posx1_dave+posx2_dave)./2;
% posy_dave = (posy1_dave+posy2_dave)./2;
% % 
% ax1 = subplot(2,2,1);
% plot(posx1(laser_state == 0),posy1(laser_state == 0),'.g');
% hold on
% plot(posx1(laser_state == 1),posy1(laser_state == 1),'.b');
% ax2 = subplot(2,2,2);
% plot(posx2(laser_state == 0),posy2(laser_state == 0),'.r');
% hold on
% plot(posx2(laser_state == 1),posy2(laser_state == 1),'.m');
% ax3 = subplot(2,2,3);
% plot(posx1_dave(laser_state == 0),posy1_dave(laser_state == 0),'.g')
% hold on
% plot(posx1_dave(laser_state == 1),posy1_dave(laser_state == 1),'.b')
% ax4 = subplot(2,2,4);
% plot(posx2_dave(laser_state == 0),posy2_dave(laser_state == 0),'.r');
% hold on
% plot(posx2_dave(laser_state == 1),posy2_dave(laser_state == 1),'.m');
% linkaxes([ax1,ax2,ax3,ax4],'xy')
% 
% 
% posx = posx_dave; 
% posy = posy_dave;
% posx1 = posx1_dave;
% posy1 = posy1_dave;
% posx2 = posx2_dave;
% posy2 = posy2_dave;


end

if led_choice == 0
    posx = posx0;
    posy = posy0;
elseif led_choice == 1
    posx = posx1;
    posy = posy1;
elseif led_choice == 2
    posx = posx2;
    posy = posy2;
elseif led_choice == 3
    posx = (posx1 + posx2)./2;
    posy = (posy1 + posy2)./2;
end


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

% figure(100)
% plot(posx,posy); hold on;
% plot(posx_interp,posy_interp)
% legend({'data before odd point removal','data after odd point removal'})
% 
% figure(101)
% plot(ts,posx)
% hold on
% plot(ts,posx_interp)
% legend({'data before odd point removal','data after odd point removal'})

posx = posx_interp;
posy = posy_interp;
speed = speed_interp;

smoothing_weights = ones(params.speed_smoothing_window_length,1)/params.speed_smoothing_window_length;
speed_smoothed = filtfilt(smoothing_weights,1,speed);

figure(99)
plot(ts,speed_smoothed);

figure()
plot(ts(laser_state==0),speed_smoothed(laser_state==0),'k'); hold on
plot(ts(laser_state==1),speed_smoothed(laser_state==1),'r')

keyboard
% Alterinatively, instead of masking out the laser times we could low pass velocity?
order = 2;
cutoff = 150;
nyquistFrequency = params.spike_sampling_rate/2;
filterCutoff = cutoff/nyquistFrequency;

[b,a] = butter(order,filterCutoff); % define the filter
% [b,a] = besself(order,filterCutoff);

filtered_speed = filtfilt(b,a,speed);
figure
plot(speed)
hold on
plot(filtered_speed)

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


% Now calculate the head direction from the 2 leds:
head_direction = calculate_headDirection(posx1,posy1,posx2,posy2);

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

figure()
plot(posx,'k');
hold on
plot(posy1,'g'); hold on;
plot(posy2,'r');
hold on
plot(speed_smoothed,'m');

title('green-red led y position, green must be on top when rat running right')

% calculate angular velocity
change_in_hd = circ_dist(head_direction(1:end-1),head_direction(2:end));
% not currently being used


Position_Data = [ts posx posy head_direction speed_smoothed x_movement y_movement];