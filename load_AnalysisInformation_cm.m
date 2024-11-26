binSize = 2;
binSize_coarse = 4;

[x_edges, y_edges, x_centers, y_centers, x_edges_coarse, y_edges_coarse, x_centers_coarse, y_centers_coarse]...
    = compute_bin_edges(0,Experiment_Information.maze_size,0,Experiment_Information.maze_size,binSize,binSize_coarse);

BetaFreqRange = [12,30];
DeltaFreqRange = [2,4];
HighGammaFreqRange = [60,100];
ind_NaN_unvistedBins = {};
ind_NaN_unvistedBins_coarse = {};
isolationThr = 0.95;
LowGammaFreqRange = [30,60];
meanFiringRateThr = 0.01;
noiseOverlapThr = 0.03;
numSpatialBins = [length(y_centers),length(x_centers)];
numSpatialBins_coarse = [length(y_centers_coarse),length(x_centers_coarse)];
numSpatialBins_coarse_smoothing = 2;
numSpatialBins_smoothing = 8;
peakSnrThr = 1.5;
replay_dispersionThr = 12;
replay_durationThr = 0.1; %sec
sequence_deltThr = 0.05; %sec
sequence_deltxThr = 20; %cm
sequence_durationThr = 0.1; %sec
sequence_jumpThr = 20;
sequence_posteriorSpreadThr = 10;
sequence_spikeDensityThr = -Inf;
shiftSizeDecodingReplay = 0.0050;
spatialInfoThr = 0.5;
speedThr = 5;
spikeDensityStepSize = 1e-3;
SWRFreqRange = [150,250]; % John used 120-170
ThetaFreqRange = [4,12];
Times_day = [];
windowSizeDecoding_replay = 0.08;
windowSizeDecoding_timeLag = 0.02;
windowSizeDecoding_behavior = 0.4;

spikeSampRate = 30000;
lfpSampRate = 1500;

artifactTimeBackward = 0.2; %sec
artifactTimeForward = 0.2;
rippleThr = 3;
spikeDensityThr = 3;

% 
% if ~exist('Analysis_Information.mat')
    save('Analysis_Information.mat','binSize','binSize_coarse','x_edges','y_edges','x_centers','y_centers','x_edges_coarse','y_edges_coarse','x_centers_coarse','y_centers_coarse'...
    ,'BetaFreqRange','DeltaFreqRange','HighGammaFreqRange','ind_NaN_unvistedBins','ind_NaN_unvistedBins_coarse','isolationThr','LowGammaFreqRange','meanFiringRateThr','noiseOverlapThr'...
    ,'numSpatialBins','numSpatialBins_coarse','numSpatialBins_coarse_smoothing','numSpatialBins_smoothing','peakSnrThr','replay_dispersionThr','replay_durationThr','sequence_deltThr'...
    ,'sequence_deltxThr','sequence_durationThr','sequence_jumpThr','sequence_posteriorSpreadThr','sequence_spikeDensityThr','shiftSizeDecodingReplay','spatialInfoThr','speedThr'...
    ,'spikeDensityStepSize','SWRFreqRange','ThetaFreqRange','Times_day','windowSizeDecoding_replay','windowSizeDecoding_replay','windowSizeDecoding_timeLag','windowSizeDecoding_behavior'...
    ,'spikeSampRate','lfpSampRate','artifactTimeBackward','artifactTimeForward','rippleThr','spikeDensityThr');

% else
% save('Analysis_Information.mat','binSize','binSize_coarse','x_edges','y_edges','x_centers','y_centers','x_edges_coarse','y_edges_coarse','x_centers_coarse','y_centers_coarse'...
%     ,'BetaFreqRange','DeltaFreqRange','HighGammaFreqRange','ind_NaN_unvistedBins','ind_NaN_unvistedBins_coarse','isolationThr','LowGammaFreqRange','meanFiringRateThr','noiseOverlapThr'...
%     ,'numSpatialBins','numSpatialBins_coarse','numSpatialBins_coarse_smoothing','numSpatialBins_smoothing','peakSnrThr','replay_dispersionThr','replay_durationThr','sequence_deltThr'...
%     ,'sequence_deltxThr','sequence_durationThr','sequence_jumpThr','sequence_posteriorSpreadThr','sequence_spikeDensityThr','shiftSizeDecodingReplay','spatialInfoThr','speedThr'...
%     ,'spikeDensityStepSize','SWRFreqRange','ThetaFreqRange','Times_day','windowSizeDecoding_replay','windowSizeDecoding_replay','windowSizeDecoding_timeLag','windowSizeDecoding_behavior'...
%     ,'spikeSampRate','lfpSampRate','artifactTimeBackward','artifactTimeForward','rippleThr','spikeDensityThr','-append');
% end

