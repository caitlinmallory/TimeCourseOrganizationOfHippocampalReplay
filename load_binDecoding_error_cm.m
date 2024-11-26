clearvars -except dayFiles day directory rat windows hand_clustered_only
%decodes activity (COM and spread of posterior) from sliding window in time

load Experiment_Information
load Analysis_Information
load clusters
load Position_Data

return_fullPosterior = 1;

times_list = [];
if isfield(Experiment_Information,'Run_Times')
    Run_Times = Experiment_Information.Run_Times;
    run_times_list = cat(2, Experiment_Information.Run_Times{:})';
    times_list = [times_list; run_times_list];
end
if isfield(Experiment_Information,'Sleep_Times')
    sleep_times_list = cat(2, Experiment_Information.Sleep_Times{:})';
    times_list = [times_list; sleep_times_list];
end

Times_day = [min(times_list) max(times_list)];

spikeSampRate = Experiment_Information.spikeSampRate;
spatialDim = Experiment_Information.spatialDim;


%sliding window size
windowSizeDecoding = 0.4;
shiftSizeDecoding = windowSizeDecoding;%shiftSizeDecoding*10;

name_ext = split(string(windowSizeDecoding),'.'); name_ext = name_ext(2);

%time bins
timeBins = load_timeBins(Times_day,shiftSizeDecoding*spikeSampRate,windowSizeDecoding*spikeSampRate);

if spatialDim == 1
    directionalDecoding = 1;
    numSpatialBins = [1 numSpatialBins(2)];
end

decoder_binDecoding = struct;
for sessionNum_decoder = 1:length(Run_Times)
    disp(strcat('session decoder ',num2str(sessionNum_decoder)))

    %select clusters to use for decoding (for 2D and 1D non-directional
    %decoding, ind_cluster pulls out the clusters to use. For directional
    %decoding, ind_cluster is a vector 1:n_clusters- but the rate_maps are
    %returned as all 1's if the cell didn't meet the criterion for a
    %particular map.


    if spatialDim == 2
        [ind_cluster,rateMap_smoothed]  = load_rateMapsForDecoding_2D_cm(clusters,sessionNum_decoder,Run_Times,numSpatialBins,meanFiringRateThr,spatialInfoThr,noiseOverlapThr,isolationThr,peakSnrThr,numSpatialBins_coarse);
   
   
    
    else
        [ind_cluster,rateMap_smoothed,~] = load_rateMapsForDecoding_1D_cm(clusters,sessionNum_decoder,Run_Times,meanFiringRateThr,spatialInfoThr,noiseOverlapThr,isolationThr,peakSnrThr,numSpatialBins,directionalDecoding);

        number_place_cells = length(unique([ind_cluster{1};ind_cluster{2}]));
        if exist('session_wide_properties.mat')==2
            save('session_wide_properties.mat','number_place_cells','-append')
        else
            save('session_wide_properties.mat','number_place_cells')
        end
    end

    %collect spikes within timeBin windows
    if spatialDim == 1 && directionalDecoding == 1
        numSpks = load_numSpks_timeBins(timeBins,clusters,(1:length(clusters)),shiftSizeDecoding*spikeSampRate,windowSizeDecoding*spikeSampRate);
    else
        numSpks = load_numSpks_timeBins(timeBins,clusters,ind_cluster,shiftSizeDecoding*spikeSampRate,windowSizeDecoding*spikeSampRate);
    end



    if spatialDim == 2
        [~,posteriorCOM,posteriorSpread,posteriorPeak] =  compute_BayesianDecoding_fast(numSpks,rateMap_smoothed,numSpatialBins,spatialDim,windowSizeDecoding,return_fullPosterior);
    else
        [posterior,posteriorCOM,posteriorSpread,posteriorPeak] =  compute_BayesianDecoding_fast_1D(numSpks,rateMap_smoothed,numSpatialBins,windowSizeDecoding,return_fullPosterior,directionalDecoding);
    end

    if return_fullPosterior == 1 && spatialDim == 1
        figure()
        ax1 = subplot(2,1,1);
        imagesc(posterior(:,1:size(posterior,2)/2)');
        set(gca,'YDir','normal')
        colormap(ax1,hot)
        caxis([0 0.1])
        ax2 = subplot(2,1,2);
        imagesc(posterior(:,size(posterior,2)/2+1:end)');
        set(gca,'YDir','normal')
        linkaxes([ax1,ax2],'x')
        colormap(ax2,hot)
        caxis([0 0.1])
    end


%load into struct and save
decoder_binDecoding(sessionNum_decoder).sessionNum_decoder = sessionNum_decoder;
decoder_binDecoding(sessionNum_decoder).shiftSizeDecoding = shiftSizeDecoding;
decoder_binDecoding(sessionNum_decoder).windowSizeDecoding = windowSizeDecoding;
decoder_binDecoding(sessionNum_decoder).posteriorCOM = posteriorCOM;
decoder_binDecoding(sessionNum_decoder).posteriorSpread = posteriorSpread;
decoder_binDecoding(sessionNum_decoder).posteriorPeak = posteriorPeak;
decoder_binDecoding(sessionNum_decoder).timeBins = timeBins;
if return_fullPosterior == 1 && spatialDim==1
    decoder_binDecoding(sessionNum_decoder).posterior = posterior;
end

end
%     %shuffles
%     numShuffles = 10;
%     for i = 1:numShuffles
%
%         %TODO: I'm not sure that this is quite correct for the 1D
%         %directional case:
%
%         disp(strcat('shuffle ',num2str(i)))
%
%         %shuffle cell ids
%         ind_shuffle = randperm(size(rateMap_smoothed,1));
%         rateMap_smoothed_shuffle = rateMap_smoothed(ind_shuffle,:);
%
%         %decoding
%         return_fullPosterior = 0;
%
%         if spatialDim == 2
%             [~,posteriorCOM,posteriorSpread,posteriorPeak] =  compute_BayesianDecoding_fast(numSpks,rateMap_smoothed_shuffle,numSpatialBins,spatialDim,windowSizeDecoding,return_fullPosterior);
%         else
%             [~,posteriorCOM,posteriorSpread,posteriorPeak] =  compute_BayesianDecoding_fast_1D(numSpks,rateMap_smoothed_shuffle,numSpatialBins,windowSizeDecoding,return_fullPosterior,directionalDecoding);
%         end
%
%         %load into struct and save
%         decoder_binDecoding(sessionNum_decoder).shuffle(i).posteriorCOM = posteriorCOM;
%         decoder_binDecoding(sessionNum_decoder).shuffle(i).posteriorSpread = posteriorSpread;
%         decoder_binDecoding(sessionNum_decoder).shuffle(i).posteriorPeak = posteriorPeak;
%     end

save('binDecoding_error.mat','decoder_binDecoding','-v7.3')
clear numSpks posteriorCOM posteriorSpread posteriorPeak



