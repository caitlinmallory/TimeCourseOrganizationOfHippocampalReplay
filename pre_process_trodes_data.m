deeplabcut_tracking = 0;
camera_module_tracking = 1;

load params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First Pre-Process Individual Sessions. In the next section will we combined multiple sessions for analysis.
%
% session_dir = 'L:/';
% session_folder = '20220608_162649_novel_lin_track_1';
% combined_session_folder = 'L:/20220608_162649_novel_lin_track_1/processed';
% sessions_to_concatenate = {'L:/20220608_162649_novel_lin_track_1'};


session_dir = '/media/caitlin/Caitlin_Drive_6/raw/';
session_folder = 'bradTask__20210419_182343';
combined_session_folder = '/home/caitlin/Insync/caitlinmallory@berkeley.edu/Google Drive/Processed_Data/Clover/open_field/20210420';


cd(fullfile(session_dir,session_folder));

load Experiment_Information.mat
%NOTE: if the session was wireless, you need to put the environment rec file in a
%separate folder for the remaining steps:
merged_rec_file = dir(fullfile(session_dir,session_folder, '*merged.rec'));


if ~isempty(merged_rec_file)
    wireless_recording = true;
    merged_rec_file = merged_rec_file.name;
    session_prefix = strsplit(merged_rec_file,'.rec'); session_prefix = session_prefix{1};

else
    wireless_recording = false;
    env_rec_file = dir(fullfile(session_dir,session_folder, '*.rec')).name;
    session_prefix = strsplit(env_rec_file,'.rec'); session_prefix = session_prefix{1};
end

% The following should be run from the command line or from Trodes2
% application:
%% Extract Continuous Time:
% If there are are multiple rec files, you need to do a pre-processing step
% that resets time to zero for the second one.
% extractTimeBinaryFile(session_prefix);
%% Extract Spikes:
% extractSpikeBinaryFiles(session_prefix)
%save_mclust_hand_clustered_data(session_dir,session_folder)
%% Extract LFP:
% extractLFPBinaryFiles(session_prefix)
%% Save selected LFP channels as a matlab file:
LFP_channels_to_save = unique([Experiment_Information.theta_refs Experiment_Information.gamma_refs Experiment_Information.ripple_refs]);
LFP_channels_to_save = 1:64;
load_LFPconversionToMatFiles_spikeGadgets_selectTetrodes_cm(session_dir, session_folder, LFP_channels_to_save);
%% Pre-process the position data output by deeplabcut:
preprocess_deeplabcut_output
%% Process the position data
load_position_data

%% If there was no video for this session, you still may want to set the clip to use using the LFP data.
LFP_file = dir('LFP*').name;
load(LFP_file)
plot(LFP_Data(:,1),LFP_Data(:,2));
LFP_ts = LFP_Data(:,1);
hold on
if isfield(Experiment_Information,'time_clip') && ~isempty(Experiment_Information.time_clip)
    time_clip = Experiment_Information.time_clip;

    [val, closest_time_1] = min(abs(LFP_ts-time_clip(1)));
    time_clip(1) = LFP_ts(closest_time_1);
    [val, closest_time_2] = min(abs(LFP_ts-time_clip(2)));
    time_clip(2) = LFP_ts(closest_time_2);
    plot(LFP_ts(LFP_ts==time_clip(1)),LFP_Data(LFP_ts==time_clip(1),2),'ok','MarkerFaceColor', 'k')
    hold on
    plot(LFP_ts(LFP_ts==time_clip(2)),LFP_Data(LFP_ts==time_clip(2),2),'ok','MarkerFaceColor', 'k')
else
    disp('no time clip')
    time_clip = [min(LFP_ts) max(LFP_ts)];
    keyboard
end

Experiment_Information.time_clip = time_clip;
save('Experiment_Information','Experiment_Information')
disp(['time clip: ' num2str(time_clip)])

%%
% Experiment_Information.rewardLocs{1} = Well_Coordinates;
% Experiment_Information.rewardLocs{2} = Well_Coordinates;
% Experiment_Information.homeLoc{1} = Well_Coordinates(21,:);
% Experiment_Information.homeLoc{2} = Well_Coordinates(29,:);

%% Extract DIO's and get laser ON/OFF times (only works if laser TTL pulses were sent to MCU)
%extractDioBinaryFiles_cm(fullfile(session_dir,session_prefix,[session_prefix '.rec']));
load_laser_times_from_mcu(session_dir,session_folder);

%% Extract laser ON/OFF times from statescript file
% generateTimestampData_from_statescript_log(session_dir,session_folder);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Combining sessions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct a file containing a list of the last timestamps for each session, and
% and a stucture containing the time offset for that session (in other
% words, you will add that value to the timestamps for that session).

% Another way to get time offsets is to put the rec files (for merged
% sessions, use the environment recording) into one folder, all with the
% same prefix (i.e., 'mysession_01.rc','mysession_02.rec') and run
% extractTimeFile('mysession'); The time offset that is generated using
% this method adds together the last timestamp from the previous session +
% the actual time (in timestamps) between the end of the previous session
% and the start of the next). This code would need to be changed to put
% those values into Session_Information.time_offset.


% Specify the folder containing the combined data:
%combined_session_folder = 'F:/20210420';

cd(combined_session_folder)

load Experiment_Information.mat
% List the sessions you want to combine IN ORDER:
load session_list;

% Create the Analysis_Information file (check that all entries are correct
% before running this!)
load_AnalysisInformation_cm


for i = 1:size(sessions_to_concatenate,1)

    lfp_data_path = dir(fullfile(sessions_to_concatenate{i},'*.LFP')).name;
    lfp_ts_file = dir(fullfile(sessions_to_concatenate{i},lfp_data_path, '*.timestamps.dat')).name;
    formatSpec = 'Loading lfp data from %s\n';
    fprintf(formatSpec,lfp_data_path)
    lfp_ts_info = readTrodesExtractedDataFile(fullfile(sessions_to_concatenate{i},lfp_data_path,lfp_ts_file));

    Session_Information(i).name = sessions_to_concatenate{i};
    Session_Information(i).first_timestamp = double(min(lfp_ts_info.fields.data));
    Session_Information(i).last_timestamp = double(max(lfp_ts_info.fields.data));
    Session_Information(i).last_timestamp_plus_decimation = double(max(lfp_ts_info.fields.data) + lfp_ts_info.decimation);
end

for i = 2:length(sessions_to_concatenate)
    Session_Information(i).time_offset = sum([Session_Information(1:i-1).last_timestamp_plus_decimation]);
end
Session_Information(1).time_offset = 0;

if ~exist(fullfile(combined_session_folder,'Session_Information.mat'))
    save(fullfile(combined_session_folder,'Session_Information.mat'), 'Session_Information')
else
    save(fullfile(combined_session_folder,'Session_Information'),'Session_Information','-append');
end

%% Load the experiment information for each session. To the structure Session_Information add the run timestamp information for each session (without adding the time offsets).

for i = 1:length(sessions_to_concatenate)
    load(fullfile(sessions_to_concatenate{i},'Experiment_Information.mat'));
    if strcmp(Experiment_Information.epoch,'run')
        session_run_times = Experiment_Information.time_clip;
        Session_Information(i).Run_Times = session_run_times;
    else
        Session_Information(i).Run_Times = [];
    end

    if strcmp(Experiment_Information.epoch,'sleep')
        session_sleep_times = Experiment_Information.time_clip;
        Session_Information(i).Sleep_Times = session_sleep_times;
    else
        Session_Information(i).Sleep_Times = [];
    end

end

if ~exist(fullfile(combined_session_folder,'Session_Information.mat'))
    save(fullfile(combined_session_folder,'Session_Information','Session_Information'))
else
    save(fullfile(combined_session_folder,'Session_Information'),'Session_Information','-append');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Combine position data, LFP data, spike, and laser timestamp data across multiple sessions
% Here we will concatenate the sessions. The key is to update the
% timestamps by adding the time offset associated with the session- this
% will ensure that timestamps are strictly increasing across sessions.


% 1. Combine LFP matfiles into one long file with updated timestamps (i.e., with offsets added)
combined_session_folder = pwd
load Experiment_Information;
load Session_Information;
load Analysis_Information;
load session_list;
clear LFP_Data


lfp_tetrodes = unique([Experiment_Information.theta_refs,Experiment_Information.ripple_refs,Experiment_Information.gamma_refs]);

for i = 1:length(lfp_tetrodes)
    LFP_Data = [];
    for j = 1:length(sessions_to_concatenate)
        session_lfp = load(fullfile(sessions_to_concatenate{j},['LFP' num2str(lfp_tetrodes(i))])).LFP_Data;

        % WARNING: there are repeated LFP timestamps in some sessions. I am
        % currently deleting the repeats, but there might be a better way
        % to deal with this issue.
        [~,rows] = unique(session_lfp(:,1));
        session_lfp = session_lfp(rows,:);

        % Add time offset to the LFP times!
        session_lfp(:,1) = session_lfp(:,1)+Session_Information(j).time_offset;
        LFP_Data = [LFP_Data; session_lfp];

    end
    save(['LFP_Data' num2str(lfp_tetrodes(i))],'LFP_Data')

    % Also create a denoised copy of LFP:
    load(['LFP_Data' num2str(lfp_tetrodes(i))]); %

    artifact_idx = remove_lfp_artifacts_cm(LFP_Data(:,2),LFP_Data(:,1),Experiment_Information.artifact_threshold, artifactTimeBackward, artifactTimeForward, lfpSampRate, 1);
    LFP_Data_denoised = LFP_Data;
    LFP_Data_denoised(artifact_idx,2) = NaN;

    % check if there is a file called "bad_lfp" (made by hand). If so, removed
    % the LFP at the indices indicated.
    if exist('bad_lfp.mat')==2
        load('bad_lfp.mat')
        for n = 1:size(bad_lfp,1)
            LFP_Data_denoised(bad_lfp(n,1):bad_lfp(n,2),2) = nan;
        end
        hold on
        plot(LFP_Data_denoised(:,1),LFP_Data_denoised(:,2))
    end


    save(['LFP_Data' num2str(lfp_tetrodes(i)) '_denoised'],'LFP_Data_denoised','artifact_idx')
end
%%
% 2. Combine the position files into one long file with updated timstamps
% (i.e., with offsets added)
clear Position_Data
Position_Data = [];
for i = 1:length(sessions_to_concatenate)
    if exist(fullfile(sessions_to_concatenate{i},'Position_Data.mat'))
        session_position_data = load(fullfile(sessions_to_concatenate{i},'Position_Data.mat')).Position_Data;
        % Add time offset to position times!
        session_position_data(:,1) = session_position_data(:,1) + Session_Information(i).time_offset;
        Position_Data = [Position_Data; session_position_data];
    end

end
save(fullfile(combined_session_folder,'Position_Data'),'Position_Data')
%%
% 3. Update the Run Times and Sleep Times for the session by adding the time offsets
Experiment_Information.Run_Times = {};
Experiment_Information.Sleep_Times = {};
Experiment_Information.Run_Decoder = {};
Experiment_Information.Sleep_Decoder = {};
run_num = 1;
sleep_num = 1;
for i = 1:length(sessions_to_concatenate)

    if ~isempty(Session_Information(i).Run_Times)
        Experiment_Information.Run_Times{run_num,1} = Session_Information(i).Run_Times + Session_Information(i).time_offset;
        Experiment_Information.Run_Decoder{run_num,1} = run_num;
        run_num = run_num + 1;
    end

    if ~isempty(Session_Information(i).Sleep_Times)
        Experiment_Information.Sleep_Times{sleep_num,1} = Session_Information(i).Sleep_Times + Session_Information(i).time_offset;
        sleep_num = sleep_num + 1;
    end

end

%mannually add in which decoder to use for which sleep session:

disp('Manually add in which decoder to use for which sleep session:')
disp('Manually update the Segments list!')


% add some other relavent things to Experiment_Information:
% clear Experiment_Information.time_clip;
Experiment_Information.spatialDim = 2;
Experiment_Information.posSampRate = 20;
Experiment_Information.spikeSampRate = 30000;
Experiment_Information.Track_Type = 'open_field';
Experiment_Information.trialsPerRun = 1;
save('Experiment_Information','Experiment_Information')

% clear Experiment_Information.time_clip;
% Experiment_Information.spatialDim = 1;
% Experiment_Information.posSampRate = 20;
% Experiment_Information.spikeSampRate = 30000;
% Experiment_Information.Track_Type = 'linear_track';
% Experiment_Information.trialsPerRun = 1;
%
% time_S1 = 26.66; % min
% time_S3 = 34;

% Experiment_Information.Run_Times{4} = Experiment_Information.Run_Times{2};
% Experiment_Information.Run_Times{2} = [Experiment_Information.Run_Times{1}(1) Experiment_Information.Run_Times{1}(1) + time_S1*60*30000];
% Experiment_Information.Run_Times{3} = [Experiment_Information.Run_Times{1}(1) + time_S1*60*30000 Experiment_Information.Run_Times{1}(2)];
% Experiment_Information.Run_Times{5} = [Experiment_Information.Run_Times{4}(1) Experiment_Information.Run_Times{4}(1) + time_S3*60*30000];
% Experiment_Information.Run_Times{6} = [Experiment_Information.Run_Times{4}(1) + time_S3*60*30000 Experiment_Information.Run_Times{4}(2)];
% Experiment_Information.Run_Times{7} = [Experiment_Information.Run_Times{1}(1) Experiment_Information.Run_Times{4}(2)];

% Experiment_Information.Segments(1).Times = Experiment_Information.Run_Times{2};
% Experiment_Information.Segments(2).Times = Experiment_Information.Run_Times{3};
% Experiment_Information.Segments(3).Times  = Experiment_Information.Sleep_Times{1};
% Experiment_Information.Segments(4).Times = Experiment_Information.Run_Times{5};
% Experiment_Information.Segments(5).Times = Experiment_Information.Run_Times{6};
% Experiment_Information.Segments(6).Times  = Experiment_Information.Sleep_Times{2};
%
% Experiment_Information.Segments(1).Decoder = 2;
% Experiment_Information.Segments(2).Decoder = 3;
% % Experiment_Information.Segments(3).Decoder  = 1;
% % Experiment_Information.Segments(4).Decoder = 5;
% % Experiment_Information.Segments(5).Decoder = 6;
% % Experiment_Information.Segments(6).Decoder  = 4;
%
% Experiment_Information.Segments(1).Flags = [11 19];
% Experiment_Information.Segments(2).Flags = [11 19];
% Experiment_Information.Segments(3).Flags = 12;
% Experiment_Information.Segments(4).Flags = 11;
% Experiment_Information.Segments(5).Flags = 11;
% Experiment_Information.Segments(6).Flags = 12;


Experiment_Information.Segments(1).Times = Experiment_Information.Run_Times{1};
Experiment_Information.Segments(1).Decoder = 1;
Experiment_Information.Segments(1).Flags = [14 11 10.5 300];

Experiment_Information.Segments(2).Times = Experiment_Information.Run_Times{2};
Experiment_Information.Segments(2).Decoder = 2;
Experiment_Information.Segments(2).Flags = [14 11 10.5 300];

Experiment_Information.Segments(3).Times = Experiment_Information.Run_Times{3};
Experiment_Information.Segments(3).Decoder = 3;
Experiment_Information.Segments(3).Flags = [14 11 10.5 300];

save('Experiment_Information','Experiment_Information')

%%
% 5. % Load up all the laser timestamp files into one list and add the time offsets
% for each session
% NOTE: right now, 'laser times' is saved in seconds for each session, and needs to be multipled
% by 30,000 to get to timestamps. When saving to the combined session
% folder, I will leave it in timestamps- but check that this doesn't break
% other code (that is expecting times)

laser_timestamps_cat = [];
for i = 1:length(sessions_to_concatenate)
    laser_timestamp_file = dir(fullfile(sessions_to_concatenate{i},'*TimestampData*'));
    if isempty(laser_timestamp_file)
        disp(['No laser timestamps data found for session ' num2str(i)])
        continue
    else

        laser_timestamp_file = laser_timestamp_file.name;
        load(fullfile(fullfile(sessions_to_concatenate{i},laser_timestamp_file)));
        laser_times_sub = laser_times;
%         laser_times_sub(:,1) = (laser_times_sub(:,1)/1000)*params.spike_sampling_rate;
        % remove any laser timetamps that occured before or after the recording
        % segment.
        first_timestamp = Session_Information(i).first_timestamp;
        last_timestamp = Session_Information(i).last_timestamp;
        laser_times_sub(laser_times_sub(:,1)<first_timestamp | laser_times_sub(:,1)>last_timestamp,:) = [];
        % add the session time offset to the laser timestamps
        laser_times_sub(:,1) = laser_times_sub(:,1) + Session_Information(i).time_offset;

        laser_timestamps_cat = [laser_timestamps_cat; laser_times_sub];
    end
end

if ~isempty(laser_timestamps_cat)
    save('LaserTimestampData','laser_timestamps_cat');
end
%%
%%
% 4. Combine the spike data (clustered in mountainsort across multiple sessions)
% Specify the folder containing the mountainsort output: If you are
% creating a clusters file for cells clustered across sessions, use [] as
% the third input. If you are creating a clusters file for cells clustered
% in one of the sessions alone, enter the session number as the third
% input

mountainsort_output_location = '/home/caitlin/Data/clover_20210419_session1/preprocessing';
% Run the following if you used the Franklab pipeline to sort individual
% sessions and then look for overlapping clusters:
load_clusterStruct_mountainsort_tagged_cm(mountainsort_output_location,combined_session_folder,[])

% Run the following if you ran Mountainsort on one huge MDA file with
% sessions concatentated:
load_clusterStruct_Mountainsort_combined_clustering_cm
% you can combine your clusters sorted separately for each session using the following syntax: clusters = [clusters1,clusters2]


%If hand clustered in MClust:
% if you need to combine spikes across sessions (i.e., you didn't export them together in trodes):
time_offsets = [Session_Information.time_offset]';
concatenate_trodes_spike_data_for_mclust(sessions_to_concatenate,time_offsets)


% if the data was clustered separtely for each session:
load_clusterStruct_MClust_cm(combined_session_folder,sessions_to_concatenate,'/media/caitlin/Caitlin_Drive_6/raw/bradTask__20210419_182343/bradTask__20210419_182343_merge.spikes')
%OR



% if the data was clustered on concatenated sessions:
sub_sessions_to_concatenate = sessions_to_concatenate;
load_clusterStruct_MClust_combined_clustering_cm(combined_session_folder,...
    '/media/caitlin/Caitlin_Drive_6/raw/bradTask__20210419_182343/bradTask__20210419_182343_merge.spikes',sessions_to_concatenate)

% ADD ISOLATION METRICS TO MCLUST DATA:
load_cluster_isolation_metrics_mclust(combined_session_folder,'/media/caitlin/My Passport/CM1/20230228_142126/merged_spikes_for_sessions_1-4')

% reformat_kilosort_cluster_struct converts old kilosort structure forma
% to 'clusters' format
% add_contamPct_to_kilosorted_clusters adds the % of contaminated points (a
% metric of cluster isolation in kilosort2)



%% The following requires LFP and Position Data:
cd (combined_session_folder)
zscore_lfp_power % saves scored LFP (beta, theta, gamma, ripple) for the session NEW
zscore_lfp_power_less_smoothing % just saves ripple power z-scored using sigma 1.5 ms
load_laserState
%linear track sessions:
segment_linear_track_passes
% note that there is an issue when using this code on sessions where the
% laser was on during the run
analyze_drinking_period_behavior % adds measures of duration and velocity to each drink period; saves in 'Behavior_Data'
% LFP and session-wide analyses:
load_LFP_power_during_run % compare LFP power in different frequency bands between laser on and off run laps
load_LFP_power_during_reward % compare LFP power in different frequency bands between laser on and off reward zone epochs

%% The following require spikes:
% save_clusters_spike_times_only % load clusters, save a slimmed down
% version containing just spike times [Not working correctly]
load_spikes_behavior_cm % Add position data to the Clusters structure (position of animal at each
% spike time)
load_clusters_rate_maps_cm % Generate rate maps for each session and attach them to the clusters
% struct;
load_rate_maps_by_lap % Generates rate maps for different lap categories (laser on-leftward, laser off-leftward etc) and attaches them to clusters
load_excitatory_designation % add a flag describing if the cell is excitatory or not
load_well_isolated_designation % add a flag describing if the cell is well isolated or not
load_place_field_correlation_across_maps % for each well isolated cell in a session, look at the rate correlation and Pearson's correlation between place maps for each running direction. Save the average to 'session_wide_properties'
% add theta phase to clusters
load_clusters_theta_phase % [TODO make sure this is working correctly! I had introduced an error and fixed it within ripple phase below[
load_clusters_ripple_phase % NEW
% and Modulation to clusters


% analyze_phase_precession
load_firing_rate_during_reward_epochs % computes firing rate in all stationary epochs
load_firing_rate_during_run_epochs

load_spikeDensity_cm % loads and saves spike density across the session. 
load_spikeDensity_pyramidal_only % loads and saves spike density across the session, only using excitatory cells.
load_spikeDensity_inhibitory_only% loads and saves spike density across the session, only using inhibitory cells.

zscore_spikeDensity_pyr(0.0125); 
zscore_spikeDensity_pyr(0.0015);
zscore_spikeDensity_inhibitory(0.0125);
zscore_spikeDensity_inhibitory(0.0015);


load_binDecoding_error_cm
plot_binDecoding_error_cm 
plot_binDecoding_error_segments_cm; % Adds decoding error to Segments
load_binDecoding_cm(0.08,0.005); %?? What did I actually use for the 1D stuff?
load_binDecoding_cm(0.02,0.005); %?? What did I actually use for the 1D stuff?
load_candidateEventTimes % gets ripple times and spike density event times. Also saves the z-scored LFP for the session.
load_theta_phase_locking %For each cell, computes the degree of theta modulation, assigns clusters as 'unimodal' or 'bimodal'

do_filter_candidate_events % Uses a modified version of John's method to find smooth segments within a candidate event
%load_replayEvents_cm % uses John's method to find segments of smooth posterior with minimal criterion

load_candidateEvents_cm_new % computes metrics associated with each candidate event (ripple, spike density, or John's filtering method).
%addLFP_to_decoderEvents % adds 100 ms of raw LFP surrounding the time of peak ripple power to each event in decoder_events

load_session_wide_ripple_sde_rate; %computes the rates of ripples and sdes for a session and adds those metrics to session_wide_properties

% Cell analyses (add's metrics to clusters)
load_ripple_phase_locking %For each cell, computes the degree of ripple modulation
%load_ripple_phase_precession
load_cell_participation_in_ripples_sdes_or_replays %looks at each cell's participation in different ripple types

% Replay analyses:
load_replay_by_stopping_period % counts how many events occured in time bins following reward onset
load_replay_by_stopping_period_aligned_to_departure_v2 % counts how many events occured in time bins aligned to rat's departure from endzone

% Perform on individual sessions:
load_replays_from_individual_session; % function to load replays
analyze_ripple_frequency_across_replay_types
load_LFP_surrounding_ripples % given a list of ripple times, pulls out the LFP around the ripple
plot_place_fields_laser_off_v_laser_on % for sessions where the laser was on during the run, plots place fields in laser off and laser on laps

% group summaries:
combine_replayPropertiesAcrossSessions % compare properties of laser on and off replays
combine_ripple_phase_locking_in_replays % looks at ripple phase locking across cells
combine_theta_phase_modulation

combine_LFP_band_power_while_running % looks at theta, gamma, beta, ripple power during laser on/off run laps
combine_LFP_band_power_while_drinking % looks at theta, gamma, beta, ripple power during reward zone epochs


combine_unimodal_bimodal_properties % looks at basic features of bimodal/unimodal cells
combine_replay_by_stopping_period % visualize replay rate aligned to drink onset (combines across sessions)
combine_replay_by_stopping_period_aligned_to_departure % visualize replay rate aligned to drink onset (combines across sessions)
combine_participation_in_ripples_sdes_replays % look at whether certain cells participate more in certain kinds of replays. assess the impact of laser

combine_firing_rate_during_reward_epochs % compares cells' firing rates during laser off or laser on reward epochs (stationary, within 30 cm of the end-zone);
combine_firing_rate_during_run_epochs % compares cells' firing rates during laser off or laser on run epochs (moving >10 cm/s, not within 50 cm of the end-zone);


% plotting:
plot_candidateEvents_mosaic_1D_directional_cm_only_good_replays
plot_replayEvents_mosaic_cm

% in progress
analyze_ripple_frequency_across_replay_types % currently a mess- but this code runs RhythmSOM_Classifier from Liset M. de la Prida. Also plots spectrograms of ripples of different typese (i.e., reverse, forward, laser on, laser off).
run_pca_on_ripples_and_cluster_according_to_content% perform pca on ripples, plot ripples according to first 2 PC's, color according to replay type
append_candidateEvents_cm % slimmed down version of load_candidateEvents_cm so you can add a new features
in_progress_compare_behavior_laser_on_laser_off % not exactly sure but I think this is looking at behavioral features of on v. off stopping periods
in_progress_4_17_23 % also seems to be comapring behavior between laser on and off epochs

% Analyses:
compare_novelty_measures % combines data from all sessions, looks at the effect of novelty on rate-map correlations for the 2 running directions







% old
% compute_rippleProperties
% compute_cellParticipationInRipples
% compute_firingRateInRunLaserEpochs
% combine_replay_by_stopping_period
% compare_spikeDensityPropertiesAcrossSessions
% compare_ripplePropertiesAcrossSessions
% compare_cellParticipationPropertiesAcrossSessions_v2