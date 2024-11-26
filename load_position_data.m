% load laser_state
load Experiment_Information.mat
deeplabcut_tracking = 0;
camera_module_tracking = 1;

load params
session_path = strsplit(pwd,'/');
session_dir = fullfile('/',session_path{1:end-1});
session_folder = session_path{end};

%% Pre-process the position data:
% If there are going to be multiple sessions in the day, it's very
% important that maze_processing_params stays the same for each behavioral
% session. The first time, you can run with maze_processing_params empty:
% this will generate 'maze_processing_params_new' which you can then save
% and apply to all subsequent sessions.

if strcmp(Experiment_Information.epoch,'run')
    params.smooth_position_data = 1;
    params.speed_smoothing_window_length = 3;
    params.maze_size = Experiment_Information.maze_size;
    try 
        params.track_type = Experiment_Information.track;
    catch
        params.track_type = Experiment_Information.Track_Type;
    end

    if camera_module_tracking == 1
        video_tracking_file = dir(fullfile(session_dir,session_folder,'*.videoPositionTracking')).name;
        position_data = readTrodesExtractedDataFile(video_tracking_file);
        
        ts = double(position_data.fields(1).data);
        posx = double(position_data.fields(2).data);
        posy = double(position_data.fields(3).data);
        
    elseif deeplabcut_tracking == 1
        if exist('preprocessed_position_data.mat') ~=2
        preprocess_deeplabcut_output
        else
        load(fullfile(session_dir,session_folder,'preprocessed_position_data'))
        end
    end

    try
        load Raw_Well_Coordinates
    catch
        Raw_Well_Coordinates = nan(36,2);
    end
%     
%     f = figure();
%     f.Position = [10 10 400 300];
%     plot(position_ts,posx)
    hold on
    if isfield(Experiment_Information,'time_clip') && ~isempty(Experiment_Information.time_clip)
        time_clip = Experiment_Information.time_clip;
        [val, closest_time_1] = min(abs(ts-time_clip(1)));
        time_clip(1) = ts(closest_time_1);
        [val, closest_time_2] = min(abs(ts-time_clip(2)));
        time_clip(2) = ts(closest_time_2);
%         plot(position_ts(position_ts==time_clip(1)),posx(position_ts==time_clip(1)),'ok','MarkerFaceColor', 'k')
%         hold on
%         plot(position_ts(position_ts==time_clip(2)),posx(position_ts==time_clip(2)),'ok', 'MarkerFaceColor', 'k')
    else
        disp('no time clip')
            % If deeplabcut was used to track position data, there may be a chunk
    % of NaN's at the beginning or end of the rat wasn't in the image.
%     time_clip = [min(position_ts(~isnan(posx)))+1 max(position_ts(~isnan(posx)))-1];

    time_clip = [min(ts) max(ts)];
    Experiment_Information.time_clip = time_clip;
    save('Experiment_Information','Experiment_Information')
    disp(['time clip: ' num2str(time_clip)])
        
    end
    % pause here to determine the actual start and stop times for the run and,
    % if necessary, set/save them.

    try
        load('maze_processing_params.mat')
    catch
        maze_processing_params = [];
    end

    if Experiment_Information.spatialDim == 1
        params.spatialDim = 1;
    else
        params.spatialDim = 2;
    end

    if deeplabcut_tracking == 1
        params.smooth_position_data = 0; % deeplabcut tracked position data was already smoothed!
            [Position_Data, Well_Coordinates, new_maze_processing_params] = process_position_data_open_field_deeplabcut(ts, posx_0a, posy_0a, ...
        posx_1a,posy_1a,posx_2a,posy_2a,time_clip, Raw_Well_Coordinates, params, maze_processing_params,laser_state);
    else
        params.smooth_position_data = 1;
            [Position_Data, Well_Coordinates, new_maze_processing_params] = process_position_data_open_field(ts, posx, posy, nan(size(posx)), ...
            time_clip, Raw_Well_Coordinates, params, maze_processing_params);
    end

    if isempty(maze_processing_params)
        maze_processing_params = new_maze_processing_params;
        save('maze_processing_params','maze_processing_params')
    end

    % save processed position data for use in matlab and python:
    timestamps = Position_Data(:,1);
    posx = Position_Data(:,2);
    posy = Position_Data(:,3);
    save('Position_Data.mat','Position_Data','timestamps','posx','posy');
    
    % Also the processed Well coordinate data (if it exists) in npy format
    if nansum(nansum(Well_Coordinates)) > 0
        Well_Coordinates_X = Well_Coordinates(:,1);
        Well_Coordinates_Y = Well_Coordinates(:,2);
        save('Well_Coordinates','Well_Coordinates','Well_Coordinates_X','Well_Coordinates_Y');
    end
end 


