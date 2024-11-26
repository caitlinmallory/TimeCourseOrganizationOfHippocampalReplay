function [] = load_binDecoding_cm(windowSizeDecoding,shiftSizeDecoding)


clearvars -except dayFiles day directory rat windows hand_clustered_only windowSizeDecoding shiftSizeDecoding

plot_figure = 0;

load Experiment_Information
load Analysis_Information
load clusters
load spikeDensity
load Position_Data

%sliding window size
name_ext = split(string(windowSizeDecoding),'.'); name_ext = name_ext(2);

%load time bounds
% a = spikeDensity; a(a(:,2)==0,:) = [];
% Times_day = [min(a(:,1)) max(a(:,1))];

Run_Times = Experiment_Information.Run_Times;
Sleep_Times = Experiment_Information.Sleep_Times;
times_list = cat(2, (cat(2, Run_Times{:})), Sleep_Times{:});
Times_day = [min(times_list) max(times_list)];
spikeSampRate = Experiment_Information.spikeSampRate;
spatialDim = Experiment_Information.spatialDim;


%time bins
%timeBins = load_timeBins(Times_day,shiftSizeDecoding*spikeSampRate,windowSizeDecoding*spikeSampRate);
timeBins = load_timeBins_cm(Times_day,shiftSizeDecoding*spikeSampRate,windowSizeDecoding*spikeSampRate);

if spatialDim == 1
    directionalDecoding = 1;
    numSpatialBins = [1 numSpatialBins(2)];
end

for sessionNum_decoder = 1:size(Run_Times,1)

    %select clusters to use for decoding
    if spatialDim == 2
        [ind_cluster,rateMap_smoothed]  = load_rateMapsForDecoding_2D_cm(clusters,sessionNum_decoder,Run_Times,numSpatialBins,meanFiringRateThr,spatialInfoThr,noiseOverlapThr,isolationThr,peakSnrThr,numSpatialBins_coarse);
    else
        [ind_cluster,rateMap_smoothed,resort_index] = load_rateMapsForDecoding_1D_cm(clusters,sessionNum_decoder,Run_Times,meanFiringRateThr,spatialInfoThr,noiseOverlapThr,isolationThr,peakSnrThr,numSpatialBins,directionalDecoding);
    end

    rateMap_sub = rateMap_smoothed;

    %collect spikes within timeBin windows
    if spatialDim == 1 && directionalDecoding == 1
        numSpks = load_numSpks_timeBins(timeBins,clusters,(1:length(clusters)),shiftSizeDecoding*spikeSampRate,windowSizeDecoding*spikeSampRate);
    else
        numSpks = load_numSpks_timeBins(timeBins,clusters,ind_cluster,shiftSizeDecoding*spikeSampRate,windowSizeDecoding*spikeSampRate);
    end

    %decoding
    if spatialDim == 2
        return_fullPosterior = 0;
    else
        return_fullPosterior = 1;
    end

    if spatialDim == 2
        [~,posteriorCOM,posteriorSpread,posteriorPeak] =  compute_BayesianDecoding_fast(numSpks,rateMap_sub,numSpatialBins,spatialDim,windowSizeDecoding,return_fullPosterior);
    else
        [posterior,posteriorCOM,posteriorSpread,posteriorPeak] =  compute_BayesianDecoding_fast_1D(numSpks,rateMap_sub,numSpatialBins,windowSizeDecoding,return_fullPosterior,directionalDecoding);
    end
    %load into struct and save
    decoder_binDecoding(sessionNum_decoder).sessionNum_decoder = sessionNum_decoder;
    decoder_binDecoding(sessionNum_decoder).shiftSizeDecoding = shiftSizeDecoding;
    decoder_binDecoding(sessionNum_decoder).windowSizeDecoding = windowSizeDecoding;
    decoder_binDecoding(sessionNum_decoder).posteriorCOM = posteriorCOM;
    decoder_binDecoding(sessionNum_decoder).posteriorSpread = posteriorSpread;
    decoder_binDecoding(sessionNum_decoder).posteriorPeak = posteriorPeak;
    decoder_binDecoding(sessionNum_decoder).timeBins = timeBins;
    decoder_binDecoding(sessionNum_decoder).binSize = binSize;
    if return_fullPosterior == 1
        decoder_binDecoding(sessionNum_decoder).posterior = posterior;

        if plot_figure==1
            figure()
            ax1 = subplot(2,1,1);
            imagesc(posterior(:,1:size(posterior,2)/2)');
            set(gca,'YDir','normal')
            colormap hot
            caxis([0 0.1])
            ax2 = subplot(2,1,2);
            imagesc(posterior(:,size(posterior,2)/2+1:end)');
            set(gca,'YDir','normal')
            linkaxes([ax1,ax2],'x')
            colormap hot
            caxis([0 0.1])
        end
    end

    dt=whos('decoder_binDecoding');
    if dt.bytes < 2e+09
        save(strcat('binDecoding_',name_ext,'.mat'),'decoder_binDecoding')
    else
        save(strcat('binDecoding_',name_ext,'.mat'),'decoder_binDecoding','-v7.3')
    end

    clear numSpks posteriorCOM posteriorSpread posteriorPeak
end

close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BINNING BY PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%variable del_t, fixed del_phi

%     %prepare LFP data for extracting phase bins
%         [~,thetaLocs] = findpeaks(LFP_filtered(:,3),'minpeakheight',2*pi-0.5);
%         LFP_theta = [LFP_filtered(thetaLocs,1),2*pi*(0:length(thetaLocs)-1)'];
%
%     %phase bins
%         startPhases = (startPhase:phiShiftSize:(max(LFP_theta(:,2))-phiStepSize))';
%         endPhases = startPhases + phiStepSize;
%
%     %convert phase bins to time bins
%         starttimes = compute_dataInterpolation([LFP_theta(:,2) LFP_theta(:,1)],startPhases,[]); starttimes = starttimes(:,2);
%         endtimes = compute_dataInterpolation([LFP_theta(:,2) LFP_theta(:,1)],endPhases,[]); endtimes = endtimes(:,2);
%
%     %eliminate bins outside times
%         ind = [];
%         for i = 1:size(times,1)
%             ind = [ind; find(starttimes>times(i,1) & starttimes<times(i,2))];
%         end
%         starttimes = starttimes(ind);
%         endtimes = endtimes(ind);
%
%     %convert to cell (ith cell contains time windows that contribute to ith step)
%         timeBins = cell(length(starttimes)-2*pi/phiShiftSize*(numCycToAvg-1),1);
%         for i = 1:length(starttimes)-2*pi/phiShiftSize*(numCycToAvg-1)
%             ind = i:2*pi/phiShiftSize:i+2*pi/phiShiftSize*(numCycToAvg-1);
%             timeBins{i} = [starttimes(ind) endtimes(ind)];
%
% %             plot(LFP_filtered(1:300,1),LFP_filtered(1:300,3),'.-'), hold on, vline(timeBins{i}(:,1),'k'), vline(timeBins{i}(:,2),'r--'), hold off
% %             keyboard
%         end
%
%     %decoding
%         [timeBins,~,posteriorCoverage,decodedPosPeak,decodedPosCOM,numSpks,clusterIDs] = compute_BayesianDecoding(timeBins,clusters,numSpatialBins,spatialFiringRateThr,spatialCoverageThr,posteriorCoverageThr,0);
%
%     %saving
%         save('binDecoding_phase','timeBins','decodedPosPeak','decodedPosCOM','numSpks','clusterIDs','posteriorCoverage')

