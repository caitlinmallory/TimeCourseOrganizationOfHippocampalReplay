%returns the state of the laser (0 = off, 1 = on) at the times specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except dayFiles day directory rat windows hand_clustered_only
load Experiment_Information
load Analysis_Information


spikeSampRate = 30000;


times_list = [];
if isfield(Experiment_Information,'Run_Times')
    run_times_list = cat(2, Experiment_Information.Run_Times{:})';
    times_list = [times_list; run_times_list];
end
if isfield(Experiment_Information,'Sleep_Times')
    sleep_times_list = cat(2, Experiment_Information.Sleep_Times{:})';
    times_list = [times_list; sleep_times_list];
end

Times_day = [min(times_list) max(times_list)];


if isfile('LaserTimestampData.mat')
    load LaserTimestampData
    
    
    % get the timestamps of the laser pulses (within the above position
    % segment); *assumes that laser always started off and ended off. add in
    % some checks about this.
    % if strcmp(EPOCH,'run')
    if laser_timestamps_cat(1,2) == 0
        laser_timestamps_cat(1,:) = [];
    end


    % If multiple sessions were stiched together, the second or 3rd session
    % may have a zero for the first entry. If there are two zeros in a row,
    % delete the second one.
%     differences = diff(laser_timestamps_cat(:,2));
%     if any(differences==0)
%     keyboard
%     end


    laser_timestamps_cat = [laser_timestamps_cat(laser_timestamps_cat(:,2) == 1,1) laser_timestamps_cat(laser_timestamps_cat(:,2) == 0,1)];
    %laser_timestamps_cat = (laser_timestamps_cat/1000)*params.spike_sampling_rate;
    
    laser_timestamps_cat(laser_timestamps_cat(:,2) < Times_day(1),:) = [];
    laser_timestamps_cat(laser_timestamps_cat(:,1) > Times_day(2),:) = [];
    
%     laser_off_timestamps = [0,0];
%     for pulse = 1:size(laser_timestamps_cat,1)
%         laser_off_timestamps(pulse,2) = laser_timestamps_cat(pulse,1)-1;
%     end
%     laser_off_timestamps = [laser_off_timestamps; [0 Times_day(2)]];
%     
%     for pulse = 1:size(laser_timestamps_cat,1)
%         laser_off_timestamps(pulse + 1,1) = laser_timestamps_cat(pulse,2)+1;
%     end
%     laser_off_timestamps(1,1) = Times_day(1);
    
    
    laser_state_times = [Times_day(1):0.001*spikeSampRate:Times_day(2)]';
    start_inds = find(histc(laser_timestamps_cat(:,1),laser_state_times) == 1);
    stop_inds = find(histc(laser_timestamps_cat(:,2),laser_state_times) == 1);
    if length(start_inds) ~= length(stop_inds)
        disp('warning: laser was either on at the start or end of position data tracking')
    end
    if length(start_inds)<length(stop_inds) % the laser was already on when the position tracking started
        start_inds = [1; start_inds];
    end
    if length(start_inds)>length(stop_inds) % the laser was stil on when the position tracking started
        stop_inds = [stop_inds;  length(laser_state_times(:,1))];
    end

    laser_state = zeros(size(laser_state_times));
    
    for i = 1:length(start_inds)
        laser_state(start_inds(i):stop_inds(i)) = 1;
    end
    
    laser_state = [laser_state_times laser_state];
    save('laser_state','laser_state')
    
    % Super slow but maybe slightly more acccurate?
    % laser_state_times = [Times_day(1):0.001*spikeSampRate:Times_day(2)]';
    % laser_state = zeros(size(laser_state_times));
    %
    % for i = 1:length(laser_timestamps_cat)
    %     [val, start_ind] = min(abs(laser_state_times-laser_timestamps_cat(i,1)));
    %     [val, stop_ind] = min(abs(laser_state_times-laser_timestamps_cat(i,2)));
    %     laser_state(start_ind:stop_ind) = 1;
    % end
    %
    % laser_state = [laser_state_times laser_state];
    % save('laser_state','laser_state')
    
else
    laser_state_times = [Times_day(1):0.001*spikeSampRate:Times_day(2)]';
    laser_state = [laser_state_times zeros(size(laser_state_times))];
    save('laser_state','laser_state')
end