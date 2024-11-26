function [lfp_ts, lfp_data] = load_lfp_data(session_dir, session_folder, lfp_tetrodes, params)

lfp_data_path = dir(fullfile(session_dir,session_folder,'*.LFP')).name;
lfp_ts_file = dir(fullfile(lfp_data_path, '*.timestamps.dat')).name;
formatSpec = 'Loading lfp data from %s\n';
fprintf(formatSpec,lfp_data_path)
lfp_ts_info = readTrodesExtractedDataFile(fullfile(lfp_data_path,lfp_ts_file));

% load the lfp timestamp data:
lfp_ts = double(lfp_ts_info.fields.data);

% load the lfp data:
lfp_data = [];
for t_idx = 1:numel(lfp_tetrodes)
    lfp_file = dir(fullfile(lfp_data_path,  ['*.LFP_nt', num2str(lfp_tetrodes(t_idx)), 'ch*' ])).name;
    lfp_info = readTrodesExtractedDataFile(fullfile(lfp_data_path,lfp_file));
    lfp_data(t_idx,:) = double(lfp_info.fields.data).*lfp_info.voltage_scaling;
end

if size(lfp_ts,1) ~= size(lfp_data,2)
    disp('Warning! Wrong number of LFP timestamps')
end

% Sometimes there are huge weird jumps in the lfp timestamps. Pause in
% this case to investigate what's going on.
if any(diff(lfp_ts) > 0.01*params.spike_sampling_rate)
    disp('Warning! Large jumps in LFP timestamps')
    ts_diff = diff(lfp_ts);
    [max_jump, ~] = max(abs(ts_diff));
    disp(['Largest jump was ' num2str(max_jump) ' timestamps'])
    
end

% Sometimes there are jumps in the lfp timestamps. If so, interpolate
% the lfp at a consistent sampling rate:

if size(unique(diff(lfp_ts)),1) > 1 || any(unique(diff(lfp_ts)) ~= params.spike_sampling_rate/params.lfp_sampling_rate)
    disp('Warning! Inconsistent LFP timestamps')
    
    % Before interpolation, you need to remove data with duplicate timestamps
    [lfp_ts_clean, ia, ~] = unique(lfp_ts,'stable');
    
    lfp_data_clean = lfp_data(:,ia);
    
    lfp_ts_interp = min(lfp_ts_clean):(params.spike_sampling_rate/params.lfp_sampling_rate):max(lfp_ts_clean);
  
    lfp_data_interp = nan(size(lfp_data,1),size(lfp_ts_interp,2));
    for i = 1:size(lfp_data,1)
        lfp_data_interp(i,:) = interp1(lfp_ts_clean,lfp_data_clean(i,:),lfp_ts_interp);
    end
    
    % Write over the original lfp timestamps and values:
    lfp_data = lfp_data_interp;
    lfp_ts = lfp_ts_interp';
    
end




