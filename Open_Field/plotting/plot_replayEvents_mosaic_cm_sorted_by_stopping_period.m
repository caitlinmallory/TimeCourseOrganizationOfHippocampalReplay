tic
load Experiment_Information
load Analysis_Information
load true_drink_periods.mat
load clusters
load Position_Data
load replayEvents
load binDecoding_08
load laser_state

if numSpatialBins > 60
    marker_size = 1;
else
    marker_size = 3;
end

% For Figure 2:
% CM1 20230306_01
% specific_stopping_period = [37];
% specific_examples = [1 4 5 6 8 9];

% For Figure 3:
% CM1 20230306_01
% Small angle example:
% specific_stopping_period = [21];
% specific_examples = [2 3 4 6 10];

% Large angle example:
% specific_stopping_period = [41];
% specific_examples = [2 3 4 5 9];

% Supp Fig:
% CM1 20230307_02
% specific_stopping_period = [92]; specific_examples = [3];
% specific_stopping_period = [90]; specific_examples = [2];
% specific_stopping_period = [80]; specific_examples = [1];
% specific_stopping_period = [72]; specific_examples = [2];
% specific_stopping_period = [111]; specific_examples = [];
% specific_stopping_period = [107]; specific_examples = [];
% specific_stopping_period = [79]; specific_examples = [];
% specific_stopping_period = [75]; specific_examples = [];
% specific_stopping_period = [61]; specific_examples = [];
% specific_stopping_period = [56]; specific_examples = [];
% specific_stopping_period = [55]; specific_examples = [];
%
% specific_stopping_period = [51]; specific_examples = [];
% specific_stopping_period = [43]; specific_examples = [];
% specific_stopping_period = [35]; specific_examples = [];
% specific_stopping_period = [33]; specific_examples = [];
% specific_stopping_period = [17]; specific_examples = [];
% specific_stopping_period = [16]; specific_examples = [];
% specific_stopping_period = [13]; specific_examples = [];


% Clover openfield 20210419_session2
% specific_stopping_period = [13]; specific_examples = [2];
% specific_stopping_period = [9]; specific_examples = [5 11];
% specific_stopping_period = [15]; specific_examples = [4];
% specific_stopping_period = [16]; specific_examples = [];


% Clover openfield 20210420
% specific_stopping_period = [3]; specific_examples = [];

% Clover openfield 20210422
% specific_stopping_period = [15]; specific_examples = [];
% specific_stopping_period = [6]; specific_examples = [];

% CM1 20230306_02
% specific_stopping_period = [68]; specific_examples = [];

% Clover openfield 20210419_session1
% specific_stopping_period = [16]; specific_examples = [];

%Bolt
% openfield 20220524
% specific_stopping_period = [4]; specific_examples = [];

specific_stopping_period=[]; specific_examples = [];
plot_titles = 1;

fig_path = fullfile(pwd,'figures');
fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures';
%replay_dispersionThr = 12;
replay_dispersionThr = 0;
startDistThr = inf;
%bins;
pathLengthForPastFutureComparison = 100; % cm. % For visualization only
path_time = 5;
binSize = 2;
restrict_to_particular_ring = 1;
ring_range = [7 60]; % set to [1 1] to include just the first ring, etc.

angle_between_hd_past_trajectory_Thr = [0 60];
angle_between_past_future_trajectory_Thr = [0 inf];
use_mean_replay_path_angle = 1;
use_distribution_replay_path_angle = 0;
restrict_past_future_paths_to_particular_ring = 1;
past_future_ring_range = [7 30];


color_future = [ 62 150 81]./255;
color_past =   [0.6392    0.0118    0.6392];

%load positions
Position_Data_scaled = compute_locsToBins(Position_Data,numSpatialBins,x_edges,y_edges);
%load integrated path length
pathLength = compute_sequenceDistance_cumsum(Position_Data_scaled(:,2:3)); pathLength = [0; pathLength];

%plot spikeDensity/SWR amp
plot_spikeDensityPeak = 0;
replay_durationThr = 0;

% Pull out all unique run sessions in this day:
run_segments = [];
for i = 1:length(Experiment_Information.Segments)
    if ismember(11,Experiment_Information.Segments(i).Flags)
        run_segments = [run_segments; i];
    end
end


Run_Times = Experiment_Information.Run_Times;
Sleep_Times = Experiment_Information.Sleep_Times;


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
for n = 1:length(run_segments)
    run_segment_ID = run_segments(n);
    run_segment_ind = find(run_segments==run_segment_ID);
    true_drink_periods = true_drink_periods_summary(run_segment_ind).true_drink_periods;

    allPastPaths = cell(length(true_drink_periods),1);
    allFuturePaths = cell(length(true_drink_periods),1);

    %   for i = 1:length(true_drink_periods)
    if isempty(specific_stopping_period)
        example_stopping_period = 1:length(true_drink_periods);
    else
        example_stopping_period = specific_stopping_period;
    end

    for i = example_stopping_period
        [~,currentInd] = min(abs(Position_Data_scaled(:,1)-true_drink_periods(i,1)));
        pastInd_dist = find(pathLength(currentInd)-pathLength >= pathLengthForPastFutureComparison/binSize,1,'last');
        [~,pastInd_time] = min(abs(Position_Data_scaled(:,1) - (Position_Data_scaled(currentInd,1) - path_time*30000)));
        pastInd_time_distance_traveled = pathLength(currentInd)-pathLength(pastInd_time); % bins
        if isempty(pastInd_dist)
            pastInd = pastInd_time;
        else
            pastInd = pastInd_dist; % let's compare behavioral trajectories of the same length.
        end

        pastPath = Position_Data_scaled(pastInd:currentInd,:);
        pastPath = pastPath(end:-1:1,:); % reverse the direction so that the first index of past path is the most recent timepoint.)
        allPastPaths{i} = [pastPath(:,2) pastPath(:,3)];

        [~,currentInd] = min(abs(Position_Data_scaled(:,1)-true_drink_periods(i,2)));
        % Find the furture path that will be considered:
        % futurePath will be the path segment beginning 10 seconds after the current timepoint, or pathLengthForPastFutureComparison cm ahead- whichever leads
        % to the larger path.
        futureInd_dist = find(pathLength - pathLength(currentInd) >= pathLengthForPastFutureComparison/binSize,1,'first');
        [val,futureInd_time] = min(abs(Position_Data_scaled(:,1) - (Position_Data_scaled(currentInd,1) + path_time*30000)));
        futureInd_time_distance_traveled = pathLength(futureInd_time) - pathLength(currentInd);
        if isempty(futureInd_dist)
            futureInd = futureInd_time;
        else
            futureInd = futureInd_dist;
        end
        futurePath = Position_Data_scaled(currentInd:futureInd,:);
        allFuturePaths{i} = [futurePath(:,2),futurePath(:,3)];

        %     plot(futurePath_upsampled(:,1),futurePath_upsampled(:,2),'g'); hold on;
        %     plot(pastPath_upsampled(:,1),pastPath_upsampled(:,2),'m'); hold on;
    end

    sessionNum_decoder = Experiment_Information.Segments(run_segment_ID).Decoder;

    figure_page_count = 0;

    if plot_spikeDensityPeak==1
        % position:
        Position_Data_sub = compute_dataTemporalConcatenation(Position_Data,load_timeBounds(unique_Session_Times(run_segment_ID,:)));

        %load spike density within session and zscore
        spikeDensity_sub = compute_dataTemporalConcatenation(spikeDensity,load_timeBounds(unique_Session_Times(run_segment_ID,:)));
        Position_Data_sub = compute_dataInterpolation(Position_Data_sub,spikeDensity_sub(:,1),[]);
        moving_ind = find(Position_Data_sub(:,5) > speedThr);
        spikeDensity_sub(moving_ind,2) = nan;
        spikeDensity_sub(:,2) = compute_zscore(spikeDensity_sub(:,2));

        %load LFP ripple amp within session and zscore
        LFP_sub = compute_dataTemporalConcatenation(LFP,load_timeBounds(unique_Session_Times(run_segment_ID,:)));
        LFP_Data_denoised_sub = compute_dataTemporalConcatenation(LFP_Data_denoised,load_timeBounds(unique_Session_Times(run_segment_ID,:)));

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
    [~,Position_Data_sub] = compute_dataTemporalConcatenation(Position_Data,unique_Session_Times(run_segment_ID,:));
    Position_Data_sub = compute_locsToBins(Position_Data_sub,numSpatialBins,x_edges,y_edges);

    %load decoder
    [ind_cluster,rateMap_smoothed,~,rateMap_smoothed_NaN] = load_rateMapsForDecoding_2D_cm(clusters,sessionNum_decoder,Run_Times,numSpatialBins,meanFiringRateThr,spatialInfoThr,noiseOverlapThr,isolationThr,peakSnrThr,numSpatialBins_coarse);
    rateMap_sub = rateMap_smoothed;

    t_sub = struct2table(decoder_replay(sessionNum_decoder).replayEvents);

    % data cleanup
    if iscell(t_sub.all_angles_between_past_future_trajectory)
        cellfun(@(x) size(x), t_sub.all_angles_between_past_future_trajectory, 'UniformOutput',false)
        mask = ~cellfun(@(x) size(x,2)==90, t_sub.all_angles_between_past_future_trajectory);
        t_sub.all_angles_between_past_future_trajectory(mask) = {nan(1,90)};
        t_sub.all_angles_between_past_future_trajectory = cell2mat(t_sub.all_angles_between_past_future_trajectory);
    end
    if restrict_past_future_paths_to_particular_ring == 1
        t_sub.angle_between_past_future_trajectory = nan(height(t_sub),1);
        t_sub.angle_between_past_future_trajectory = nanmean(abs(t_sub.all_angles_between_past_future_trajectory(:,past_future_ring_range(1):past_future_ring_range(2))),2);
    end

    mask = ~cellfun(@(x) size(x,1)==90, t_sub.angDisplacement_futPath);
    t_sub.angDisplacement_futPath(mask) = {nan(90,1)};
    mask = ~cellfun(@(x) size(x,1)==90, t_sub.angDisplacement_pastPath);
    t_sub.angDisplacement_pastPath(mask) = {nan(90,1)};

    if restrict_to_particular_ring == 1
        t_sub.meanAngDisplacement_futPath = cellfun(@(x) nanmean(abs(x(ring_range(1):ring_range(2)))), t_sub.angDisplacement_futPath);
        t_sub.meanAngDisplacement_pastPath = cellfun(@(x) nanmean(abs(x(ring_range(1):ring_range(2)))), t_sub.angDisplacement_pastPath);
    end

    %for drink_period = 1:length(true_drink_periods)
    if isempty(specific_stopping_period)
        example_stopping_period = 1:length(true_drink_periods);
    else
        example_stopping_period = specific_stopping_period;
    end

    for drink_period = example_stopping_period

        replay_sub = t_sub(t_sub.timePoints(:,1)>= true_drink_periods(drink_period,1) & t_sub.timePoints(:,1) < true_drink_periods(drink_period,2) & ...
            t_sub.startDistFromRat <= startDistThr & ...
            ~isnan(t_sub.meanAngDisplacement_futPath) & ~isnan(t_sub.meanAngDisplacement_pastPath),:);

        pastPath = allPastPaths{drink_period};
        futurePath = allFuturePaths{drink_period};
        if height(replay_sub)==0
            continue
        end
        if ~(replay_sub.angle_between_past_future_trajectory(1) >= angle_between_past_future_trajectory_Thr(1) && replay_sub.angle_between_past_future_trajectory(1) <= angle_between_past_future_trajectory_Thr(2))
            continue
        end

        %for plotting mosaic

        numPanels = 32;
        subplot_xDim = 8;
        subplot_yDim = (ceil((numPanels)/subplot_xDim));
        %figure('Position',[350 474 900 600])
        figure('Position',[350 474 950 600])
        tiledlayout(subplot_yDim,subplot_xDim,'TileSpacing','compact');
        set(gcf,'PaperUnits','inches');
        set(gcf,'PaperSize', [8 8]);
        set(gcf,'PaperPosition',[0.5 0.5 7 7])
        set(gcf,'PaperPositionMode','Manual');


        if isempty(specific_examples)
            replays_to_plot = 1:min(32,height(replay_sub));
        else
            replays_to_plot = specific_examples;
        end
        %for replay_example = 1:min(32,height(replay_sub))
        for replay_example = replays_to_plot
            %load replay properties
            indNaN = replay_sub.indNaN{replay_example};
            replay = replay_sub.replay{replay_example};
            replay_NaN = replay; replay_NaN(indNaN,:) = NaN;
            replay_NaNremoved = replay; replay_NaNremoved(indNaN,:) = [];
            timeBins = replay_sub.timeBins{replay_example};
            timePoints = replay_sub.timePoints(replay_example,:);
            ratPos = replay_sub.ratPos(replay_example,:);
            ratHD = replay_sub.ratHD(replay_example,:);
            replay_laser_state = replay_sub.laser_state{replay_example,:};
            time_since_stopping = replay_sub.time_since_real_drink_onset(replay_example);

            laser_on_ind = find([0; diff(replay_laser_state)] == 1);
            if replay_laser_state(1) == 1
                laser_on_ind = [1; laser_on_ind];
            end
            laser_off_ind = find([0; diff(replay_laser_state)] == -1);

            if plot_spikeDensityPeak==1
                spikeDensity_replay = compute_dataTemporalConcatenation(spikeDensity_sub,timePoints);
                ripplePower_replay = compute_dataTemporalConcatenation(LFP_sub,timePoints);

                spikeDensityPeak = max(spikeDensity_replay(:,end));
                ripplePowerPeak = max(ripplePower_replay(:,end));
            end

            %replay decoding
            numSpks = load_numSpks_timeBins(timeBins,clusters,ind_cluster,decoder_binDecoding(sessionNum_decoder).shiftSizeDecoding*spikeSampRate,decoder_binDecoding(sessionNum_decoder).windowSizeDecoding*spikeSampRate);
            [posterior,posteriorCOM,posteriorSpread,posteriorPeak] =  compute_BayesianDecoding_fast(numSpks,rateMap_sub,numSpatialBins,spatialDim,decoder_binDecoding(sessionNum_decoder).windowSizeDecoding,1);
            %                     jumps = compute_sequenceJumps(posteriorCOM); jumps = [jumps; jumps(end)];
            %                     indNaN = find(posteriorSpread>sequence_posteriorSpreadThr);% & jumps>sequence_jumpThr);
            posterior(indNaN,:) = NaN;
            posteriorCOM(indNaN,:) = NaN;
            %                     posteriorCOM_NaNremoved = posteriorCOM; posteriorCOM_NaNremoved(indNaN,:) = [];
            [~,posterior_color,posterior_sum,posterior_sum_binarized] = compute_flattenedDecoding(posterior,spatialDim,1);
            posterior_color = reshape(posterior_color',numSpatialBins(1),numSpatialBins(2),3);

            %plot mosaic
            set(gcf,'color','w')
            nexttile()

            image(posterior_color,'alphaData',sum(posterior_color,3)~=3), set(gca,'ydir','normal')
            axis square
            hold on, plot(replay_NaNremoved(:,1),replay_NaNremoved(:,2),'color',0.5*ones(1,3),'linewidth',1), hold off
            hold on, plot(replay_NaN(:,1),replay_NaN(:,2),'k','linewidth',1), hold on;
            %                 hold on
            %                 if ~isempty(laser_on_ind)
            %                     hold on, plot(replay(laser_on_ind,1),replay(laser_on_ind,2),'or','markerFaceColor','r','markersize', 6,'linewidth',1), hold off
            %                 end
            %
            %                 if ~isempty(laser_off_ind)
            %                     hold on, plot(replay(laser_off_ind,1),replay(laser_off_ind,2),'ok','markerFaceColor','k','markersize',6,'linewidth',1), hold off
            %                 end

            % plot_behavior_barriers_wells_johns_task

            well_list = true_drink_periods_summary(run_segment_ind).goal_well_list;
            home_well = mode(well_list);
            wells = (Experiment_Information.rewardLocs{run_segment_ind})./2;
            % mark the location of all wells with open circles
            %             for j = 1:length(wells)
            %                 well = circle(wells(j,1),wells(j,2),1);
            %                 plot(well(:,1),well(:,2),'k','LineWidth',1)
            %                 hold on
            %             end


            plot(wells(:,1),wells(:,2),'ok','LineWidth',1,'MarkerSize',marker_size)
            hold on


            % mark the location of the home well with a filled circle
            % plot(wells(home_well,1),wells(home_well,2),'ok','markerSize',10)
            plot(wells(home_well,1),wells(home_well,2),'ok','LineWidth',1,'MarkerSize',marker_size,'MarkerFaceColor','k')
            hold on

            if isfield(Experiment_Information,'barrierLocs')
                barriers = Experiment_Information.barrierLocs{run_segment_ind};
                % draw in the barriers for this session
                for j = 1:6
                    if ~isempty(barriers{j})
                        line([barriers{j}(1)./2 barriers{j}(2)./2], [barriers{j}(3)./2 barriers{j}(4)./2],'LineWidth',1,'Color','k'); hold on;
                    end
                end
            end

            xticks([]); yticks([]);

            if ~isempty(futurePath), hold on, p1 = plot(futurePath(:,1),futurePath(:,2),'color',color_future,'linewidth',1.5); p1.Color(4) = 1; hold off, end
            if ~isempty(pastPath), hold on, p2 = plot(pastPath(:,1),pastPath(:,2),'color',color_past,'linewidth',1.5); p2.Color(4) = 1; hold off, end
            plot_ratLoc
            %time = {floor(((timePoints(1)/spikeSampRate)-(unique_Session_Times(sessionNum,1)/spikeSampRate))/60),round(rem(((timePoints(1)/spikeSampRate)-(unique_Session_Times(sessionNum,1)/spikeSampRate))/60,1)*60)};
            time = {floor(time_since_stopping/60),...
                round(rem(time_since_stopping/60,1)*60)};
            if time{2}>=10
                time_string = strcat(num2str(time{1}),':',num2str(time{2}));
            else
                time_string = strcat(num2str(time{1}),':0',num2str(time{2}));
            end

            if plot_titles == 1
                p.col=[0 0 0 ; color_future; color_past];
                LimitsX = xlim; LimitsY = ylim;
                t = title([ ...
                    '\color[rgb]{' sprintf('%1.2f,%1.2f,%1.2f', p.col(1,:)) '} ' time_string, ...
                    ' \color[rgb]{' sprintf('%f,%f,%f', p.col(2,:)) '} F:' num2str(round(replay_sub.meanAngDisplacement_futPath(replay_example)),3), ...
                    ' \color[rgb]{' sprintf('%f,%f,%f', p.col(3,:)) '} P:' num2str(round(replay_sub.meanAngDisplacement_pastPath(replay_example)),3)], ...
                    'fontsize',6, 'HorizontalAlignment', 'left', 'position', [LimitsX(1), LimitsY(2)], 'Interpreter', 'tex');
            end

            if plot_spikeDensityPeak==1
                text(numSpatialBins(2),5,num2str(compute_round(spikeDensityPeak,10)),'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','Color','K','fontsize',6)
                text(numSpatialBins(2),1,num2str(compute_round(ripplePowerPeak,10)),'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom','Color','K','fontsize',6)
            end
            drawnow
        end
        % setup_figure_properties
        figure_name = ['Decoder' num2str(sessionNum_decoder) '_Session' num2str(run_segment_ID) 'stopping_period_' num2str(drink_period), 'ang_' num2str(round(replay_sub.angle_between_past_future_trajectory(1)),3)];

        set(gcf, 'renderer','painters','Color', 'white','PaperPositionMode','auto');
        saveas(gcf,fullfile(fig_path,figure_name),'pdf')

        %         figure_name = [figure_name '.png'];
        %         exportgraphics(gcf,figure_name,'Resolution',700)
        %                  export_fig(fullfile(fig_path,figure_name),'-pdf','-nocrop')
        % hgexport(gcf,[figure_name '.pdf'],figure_property); %Set desired file name
        disp([num2str(drink_period) '/' num2str(length(true_drink_periods))])
    end
end


toc

