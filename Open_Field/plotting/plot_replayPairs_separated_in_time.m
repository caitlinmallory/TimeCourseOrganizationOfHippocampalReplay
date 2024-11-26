load Experiment_Information
load Analysis_Information
% load Behavior_Analysis.mat
load true_drink_periods.mat
load clusters
load Position_Data
load replayEvents
load binDecoding_08
load laser_state

load angles_between_replays.mat
% Filter for nice look examples:
time_difference_thr = [0 3];
dispersionThr = 3;
max_start_distance_from_rat = 7.5 ;%bins = 15 cm.
 numPanels = 32;

fig_path = fullfile(pwd,'figures');
plot_maze_configuration = 0; 
% pathLengthForPastFutureComparison = 100; % cm. % For visualization only
% path_time = 5;
% binSize = 2;
% restrict_to_particular_ring = 1;
ring_range = [7 30]; % set to [1 1] to include just the first ring, etc.
%
% angle_between_past_future_trajectory_Thr = [0 180];
% use_mean_replay_path_angle = 1;
% use_distribution_replay_path_angle = 0;
% restrict_past_future_paths_to_particular_ring = 1;
% past_future_ring_range = [7 30];

% color_future = [ 62 150 81]./255;
% color_past =   [0.6392    0.0118    0.6392];

%load positions
Position_Data_scaled = compute_locsToBins(Position_Data,numSpatialBins,x_edges,y_edges);
%load integrated path length
pathLength = compute_sequenceDistance_cumsum(Position_Data_scaled(:,2:3)); pathLength = [0; pathLength];


%plot spikeDensity/SWR amp
plot_spikeDensityPeak = 0;


Run_Times = Experiment_Information.Run_Times;
Sleep_Times = Experiment_Information.Sleep_Times;

run_segments = [];
for i = 1:length(Experiment_Information.Segments)
    if ismember(11,Experiment_Information.Segments(i).Flags)
        run_segments = [run_segments; i];
    end
end

spikeSampRate = Experiment_Information.spikeSampRate;
spatialDim = Experiment_Information.spatialDim;

unique_Session_Times = vertcat(Experiment_Information.Segments.Times);

Times_day = [min(unique_Session_Times(:)) max(unique_Session_Times(:))];

times = load_timeBins_cm(Times_day,decoder_binDecoding(1).shiftSizeDecoding*spikeSampRate,decoder_binDecoding(1).windowSizeDecoding*spikeSampRate);
if ~isempty(Sleep_Times)
    Position_Data = load_positions_full(Run_Times,Sleep_Times,times,Position_Data,4);
end
Position_Data(isnan(Position_Data(:,5)),5) = 0;

if plot_spikeDensityPeak==1
    %load spike density
    load spikeDensity
    spikeDensity_smoothed = conv(spikeDensity(:,2),setUp_gaussFilt([1 1000],windowSizeDecoding_replay/spikeDensityStepSize),'same');
    spikeDensity = [spikeDensity(:,1),spikeDensity_smoothed];

    %load LFPs

    ripple_tetrode = Experiment_Information.ripple_refs;
    lfp_file_denoised = ['LFP_Data' num2str(ripple_tetrode) '_denoised.mat'];
    lfp_file = ['LFP_Data' num2str(ripple_tetrode) '.mat'];

    load(lfp_file_denoised);
    load(lfp_file);

    LFP = LFP_Data;
    clear LFP_data
    LFPSampRate = 1500;

    [~,LFP_filt_SWRAmp,~] = compute_filteredLFP(SWRFreqRange,LFP(:,2),LFPSampRate);
    LFP = [LFP,LFP_filt_SWRAmp,LFP_filt_SWRAmp];

    LFP_smoothed = conv(LFP(:,4),setUp_gaussFilt([1 1000],decoder_binDecoding(1).windowSizeDecoding/(1/LFPSampRate)),'same');
    LFP(:,4) = LFP_smoothed;
end
%% pull out all future and past trajectories from the drink periods
% For each drink period: pull out the 'past' and 'future' trajectories.
% For each drink period: pull out the 'past' and 'future' trajectories.
% for segment_ID = 1:length(run_segments)
for segment_ID = 3

    run_segment_ind = find(run_segments==segment_ID);
    true_drink_periods = true_drink_periods_summary(run_segment_ind).true_drink_periods;

    allPastPaths = cell(length(true_drink_periods),1);
    allFuturePaths = cell(length(true_drink_periods),1);


    segment_ID_decoder = Experiment_Information.Segments(segment_ID).Decoder;
    figure_page_count = 0;

    if plot_spikeDensityPeak==1
        % position:
        Position_Data_sub = compute_dataTemporalConcatenation(Position_Data,load_timeBounds(unique_Session_Times(segment_ID,:)));

        %load spike density within session and zscore
        spikeDensity_sub = compute_dataTemporalConcatenation(spikeDensity,load_timeBounds(unique_Session_Times(segment_ID,:)));
        Position_Data_sub = compute_dataInterpolation(Position_Data_sub,spikeDensity_sub(:,1),[]);
        moving_ind = find(Position_Data_sub(:,5) > speedThr);
        spikeDensity_sub(moving_ind,2) = nan;
        spikeDensity_sub(:,2) = compute_zscore(spikeDensity_sub(:,2));

        %load LFP ripple amp within session and zscore
        LFP_sub = compute_dataTemporalConcatenation(LFP,load_timeBounds(unique_Session_Times(segment_ID,:)));
        LFP_Data_denoised_sub = compute_dataTemporalConcatenation(LFP_Data_denoised,load_timeBounds(unique_Session_Times(segment_ID,:)));

        % position:
        Position_Data_sub = compute_dataInterpolation(Position_Data_sub, LFP_sub(:,1), []);

        %limit to stationary periods and mask out artifacts prior to z-scoring:
        artifact_ind = find(isnan(LFP_Data_denoised_sub(:,2)));
        moving_ind = find(Position_Data_sub(:,5) > speedThr);
        bad_inds = unique([artifact_ind; moving_ind]);

        LFP_sub(bad_inds,4) = nan;
        LFP_sub(:,4) = compute_zscore(LFP_sub(:,4));

    end

    %load Position_Data
    [~,Position_Data_sub] = compute_dataTemporalConcatenation(Position_Data,unique_Session_Times(segment_ID,:));
    Position_Data_sub = compute_locsToBins(Position_Data_sub,numSpatialBins,x_edges,y_edges);

    %load decoder
    [ind_cluster,rateMap_smoothed,~,rateMap_smoothed_NaN] = load_rateMapsForDecoding_2D_cm(clusters,segment_ID_decoder,Run_Times,numSpatialBins,meanFiringRateThr,spatialInfoThr,noiseOverlapThr,isolationThr,peakSnrThr,numSpatialBins_coarse);
    rateMap_sub = rateMap_smoothed;

    %filter sequences for replays
    %     load_ind_replayEvents_all_session_types
    %     ind_sorted = sortrows([dispersion(ind_replay) ind_replay],'descend');
    %
    %     all_replay_timePoints = mean(load_fieldVec(decoder_replay(segment_ID_decoder).replayEvents,'timePoints',2),2);



    potential_examples = find(angles_between_replays(segment_ID).all_time_differences >= time_difference_thr(1) &  angles_between_replays(segment_ID).all_time_differences<= time_difference_thr(2) & ...
        angles_between_replays(segment_ID).all_dispersion(:,1) > dispersionThr & ...
        angles_between_replays(segment_ID).all_dispersion(:,2) > dispersionThr & ...
        angles_between_replays(segment_ID).all_start_distance_from_rat(:,1) < max_start_distance_from_rat & ...
        angles_between_replays(segment_ID).all_start_distance_from_rat(:,2) < max_start_distance_from_rat);

    % sort the example by time between events:
    [sorted, sort_ind] = sort(angles_between_replays(segment_ID).all_time_differences(potential_examples));
    potential_examples = potential_examples(sort_ind);

    % Inds for the example
    angles_between_replays(segment_ID).all_comparison_inds(potential_examples(1),:);

    replays = angles_between_replays(segment_ID).replays;

    %     for i = 1:min(length(potential_examples),50)
    for i = 1:length(potential_examples)

        %for plotting moaic
        if mod(i,numPanels)==1
           
            subplot_xDim = 8;
            subplot_yDim = (ceil((numPanels)/subplot_xDim));
            figure('Position',[350 474 800 600])
            tiledlayout(subplot_yDim,subplot_xDim,'TileSpacing','compact');
            set(gcf,'PaperUnits','inches');
            set(gcf,'PaperSize', [8 8]);
            set(gcf,'PaperPosition',[0.5 0.5 7 7])
            set(gcf,'PaperPositionMode','Manual');

        end
        inds_example = angles_between_replays(segment_ID).all_comparison_inds(potential_examples(i),:);

        for pair = 1:2
            %load replay properties
            indNaN = replays.indNaN{inds_example(pair)};
            replay = replays.replay{inds_example(pair)};
            replay_NaN = replay; replay_NaN(indNaN,:) = NaN;
            replay_NaNremoved = replay; replay_NaNremoved(indNaN,:) = [];
            timeBins = replays.timeBins{inds_example(pair)};
            timePoints = replays.timePoints(inds_example(pair));

            ratPos = replays.ratPos(inds_example(pair),:);
            ratSpeed = replays.ratSpeed(inds_example(pair));
            ratHD = replays.ratHD(inds_example(pair));
            replay_laser_state = replays.laser_state{inds_example(pair)};

            time_since_stopping = replays.time_since_real_drink_onset(inds_example(pair));

            laser_on_ind = find([0; diff(replay_laser_state)] == 1);
            if replay_laser_state(1) == 1
                laser_on_ind = [1; laser_on_ind];
            end
            laser_off_ind = find([0; diff(replay_laser_state)] == -1);

            %replay decoding
            numSpks = load_numSpks_timeBins(timeBins,clusters,ind_cluster,decoder_binDecoding(segment_ID_decoder).shiftSizeDecoding*spikeSampRate,decoder_binDecoding(segment_ID_decoder).windowSizeDecoding*spikeSampRate);
            [posterior,posteriorCOM,posteriorSpread,posteriorPeak] =  compute_BayesianDecoding_fast(numSpks,rateMap_sub,numSpatialBins,spatialDim,decoder_binDecoding(segment_ID_decoder).windowSizeDecoding,1);
            posterior(indNaN,:) = NaN;
            posteriorCOM(indNaN,:) = NaN;


            [~,posterior_color,posterior_sum,posterior_sum_binarized] = compute_flattenedDecoding(posterior,spatialDim,pair);
            posterior_color = reshape(posterior_color',numSpatialBins(1),numSpatialBins(2),3);


            %plot mosaic
            set(gcf,'color','w')

            if pair == 1
                nexttile()
                hold on
            end

            image(posterior_color,'alphaData',sum(posterior_color,3)~=3), set(gca,'ydir','normal'); 
            axis square
            hold on, plot(replay_NaNremoved(:,1),replay_NaNremoved(:,2),'color',0.5*ones(1,3),'linewidth',1), hold off
            hold on, plot(replay_NaN(:,1),replay_NaN(:,2),'k','linewidth',1), hold on;
box on

if plot_maze_configuration == 1
            well_list = true_drink_periods_summary(run_segment_ind).goal_well_list;
            home_well = mode(well_list);
            wells = (Experiment_Information.rewardLocs{run_segment_ind})./2;
            % mark the location of all wells with open circles
            for j = 1:length(wells)
                well = circle(wells(j,1),wells(j,2),1);
                plot(well(:,1),well(:,2),'k','LineWidth',1)
                hold on
            end

            % mark the location of the home well with a filled circle
            plot(wells(home_well,1),wells(home_well,2),'.k','markerSize',10)

            if isfield(Experiment_Information,'barrierLocs')
                barriers = Experiment_Information.barrierLocs{run_segment_ind};
                % draw in the barriers for this session
                for j = 1:6
                    if ~isempty(barriers{j})
                        line([barriers{j}(1)./2 barriers{j}(2)./2], [barriers{j}(3)./2 barriers{j}(4)./2],'LineWidth',1,'Color','k'); hold on;
                    end
                end
            end
end
            if pair == 2
                plot_ratLoc
                hold on
            end
            xticks([]); yticks([]);

            %time = {floor(((timePoints(1)/spikeSampRate)-(unique_Session_Times(segment_ID,1)/spikeSampRate))/60),round(rem(((timePoints(1)/spikeSampRate)-(unique_Session_Times(segment_ID,1)/spikeSampRate))/60,1)*60)};
            time = {floor(time_since_stopping/60),...
                round(rem(time_since_stopping/60,1)*60)};
            if time{2}>=10
                time_string = strcat(num2str(time{1}),':',num2str(time{2}));
            else
                time_string = strcat(num2str(time{1}),':0',num2str(time{2}));
            end

            all_angles = angles_between_replays(segment_ID).all_angles(potential_examples(i),:);
            mean_angle = nanmean(abs(all_angles(ring_range(1):ring_range(2))));
            time_diff = angles_between_replays(segment_ID).all_time_differences(potential_examples(i));
            LimitsX = xlim; LimitsY = ylim;
            t = title(['t: ' num2str(round(time_diff*10)/10,2) ' angle: ' (num2str(mean_angle,3))], ...
                 'fontsize',6, 'HorizontalAlignment', 'left', 'position', [LimitsX(1), LimitsY(2)], 'Interpreter', 'tex');

            drawnow

            if mod(i,numPanels) == 0

                figure_page_count = figure_page_count + 1;
                % setup_figure_properties
                figure_name = fullfile(fig_path,['Decoder' num2str(segment_ID_decoder) '_Session' num2str(segment_ID) '_page' num2str(figure_page_count)]);

                % saveas(gcf,fullfile(fig_path,figure_name),'tiff');
                set(gcf,'PaperPositionMode','auto')
                saveas(gcf,['S' num2str(segment_ID) '_page' num2str(figure_page_count)],'jpg')
                saveas(gcf,['S' num2str(segment_ID) '_page' num2str(figure_page_count)],'pdf')
            end
        end
    end
end

% plot(angles_between_replays(segment_ID).all_time_differences(potential_examples),nanmean(abs(angles_between_replays(segment_ID).all_angles(potential_examples,ring_range(1):ring_range(2))),2),'ok')
%
% [rho,p] = nancorr(angles_between_replays(segment_ID).all_time_differences(potential_examples),nanmean(abs(angles_between_replays(segment_ID).all_angles(potential_examples,ring_range(1):ring_range(2))),2))
%

