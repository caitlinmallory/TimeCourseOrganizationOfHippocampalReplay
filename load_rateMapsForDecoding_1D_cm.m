function [ind_cluster,rateMap,resort_index] = load_rateMapsForDecoding_1D_cm(clusters,sessionNum_decoder,Run_Times,meanFiringRateThr,spatialInfoThr,noiseOverlapThr,isolationThr,peakSnrThr,numSpatialBins,directionalDecoding)
% NOTE: unlike for load_rateMapsForDecoding_2D_cm, here I'm going to remove
% cells that don't meet our criterion for inclusion in decoding by setting
% their rate maps to 1.
plot_figure=0;

isolationThr = 0;
noiseOverlapThr = inf;
peakSnrThr = 0;

MeanFiringRate_MaxThr = 20;
%Max had been set to 20
PeakFiringRate_MinThr = 0.1;
%Min had had been set to 0.1

%select which clusters to use for decoding
clusters_runs = [clusters.runs]; clusters_runs = clusters_runs(sessionNum_decoder:length(Run_Times):end);
spatialInfo = load_fieldVec(clusters_runs,'spatialInfo',1);
meanFiringRate = load_fieldVec(clusters_runs,'meanFiringRate',1);
maxSpatialFiringBinLoc = load_fieldVec(clusters_runs,'maxSpatialFiringBinLoc',1);
maxSpatialFiringRate = load_fieldVec(clusters_runs, 'maxSpatialFiringRate',1);

if directionalDecoding == 0
    rateMap = load_fieldVec(clusters_runs,'rateMap_smoothed_vectorized',numSpatialBins(2));
    peak_firing_rate = max(rateMap,[],2);
    mean_firing_rate = mean(rateMap,2);
else

    clusters_directions = [clusters_runs.directions];
    clusters_direction_1 = clusters_directions(1:2:end);
    clusters_direction_2 = clusters_directions(2:2:end);

    rateMap_direction_1 = load_fieldVec(clusters_direction_1,'rateMap_smoothed_vectorized',numSpatialBins(2));
    rateMap_direction_2 = load_fieldVec(clusters_direction_2,'rateMap_smoothed_vectorized',numSpatialBins(2));

    peak_firing_rate_direction_1 = max(rateMap_direction_1,[],2);
    peak_firing_rate_direction_2 = max(rateMap_direction_2,[],2);
    mean_firing_rate_direction_1 = mean(rateMap_direction_1,2);
    mean_firing_rate_direction_2 = mean(rateMap_direction_2,2);

end


if ~isfield(clusters,'tag') %MClust
    if directionalDecoding == 0
        ind_cluster = find(peak_firing_rate > PeakFiringRate_MinThr & mean_firing_rate < MeanFiringRate_MaxThr);
    else
        ind_cluster_direction_1 = find(peak_firing_rate_direction_1 > PeakFiringRate_MinThr & mean_firing_rate_direction_1 < MeanFiringRate_MaxThr);
        ind_cluster_direction_2 = find(peak_firing_rate_direction_2 > PeakFiringRate_MinThr & mean_firing_rate_direction_2 < MeanFiringRate_MaxThr);
    end

    %ignoring spatialInfoThr because directional fields might have low
    %spatialInformation

else %mountainsort
    tag = load_fieldVec(clusters,'tag_abbrev',1);
    noiseOverlap = load_fieldVec(clusters,'noise_overlap',1);
    isolation = load_fieldVec(clusters,'isolation',1);
    peakSnr = load_fieldVec(clusters,'peak_snr',1);

    if isfield(clusters,'tag')
        %ind_cluster = find(noiseOverlap<noiseOverlapThr & isolation>isolationThr & peakSnr>peakSnrThr & meanFiringRate>meanFiringRateThr & spatialInfo>spatialInfoThr & (strcmp(string(tag),'A')|strcmp(string(tag),'B')));
        if directionalDecoding == 0
            ind_cluster = find(noiseOverlap<noiseOverlapThr & isolation>isolationThr & peakSnr>peakSnrThr & (strcmp(string(tag),'A')|strcmp(string(tag),'B'))...
                & peak_firing_rate > PeakFiringRate_MinThr & mean_firing_rate < MeanFiringRate_MaxThr);
        else
            ind_cluster_direction_1 = find(noiseOverlap<noiseOverlapThr & isolation>isolationThr & peakSnr>peakSnrThr & (strcmp(string(tag),'A')|strcmp(string(tag),'B'))...
                & peak_firing_rate_direction_1 > PeakFiringRate_MinThr & mean_firing_rate_direction_1 < MeanFiringRate_MaxThr);
            ind_cluster_direction_2 = find(noiseOverlap<noiseOverlapThr & isolation>isolationThr & peakSnr>peakSnrThr & (strcmp(string(tag),'A')|strcmp(string(tag),'B'))...
                & peak_firing_rate_direction_2 > PeakFiringRate_MinThr & mean_firing_rate_direction_2 < MeanFiringRate_MaxThr);
        end

    else
        %ind_cluster = find(noiseOverlap<noiseOverlapThr & isolation>isolationThr & peakSnr>peakSnrThr & meanFiringRate>meanFiringRateThr & spatialInfo>spatialInfoThr);
        if directionalDecoding == 0
            ind_cluster = find(noiseOverlap<noiseOverlapThr & isolation>isolationThr & peakSnr>peakSnrThr & peak_firing_rate > PeakFiringRate_MinThr & mean_firing_rate < MeanFiringRate_MaxThr);
        else
            ind_cluster_direction_1 = find(noiseOverlap<noiseOverlapThr & isolation>isolationThr & peakSnr>peakSnrThr ...
                & peak_firing_rate_direction_1 > PeakFiringRate_MinThr & mean_firing_rate_direction_1 < MeanFiringRate_MaxThr);
            ind_cluster_direction_2 = find(noiseOverlap<noiseOverlapThr & isolation>isolationThr & peakSnr>peakSnrThr...
                & peak_firing_rate_direction_2 > PeakFiringRate_MinThr & mean_firing_rate_direction_2 < MeanFiringRate_MaxThr);
        end
    end
end



if directionalDecoding==0
    resort_index = sortrows([maxSpatialFiringBinLoc(ind_cluster) (1:length(ind_cluster))']);
    resort_index = resort_index(isnan(resort_index(:,1))==0,:);
    rateMap_for_plot = rateMap(ind_cluster,:);
    %             %Spatial tuning curves plot (used for decoding)
    if plot_figure == 1
        figure
        rateMap_norm =  rateMap_for_plot./max( rateMap_for_plot,[],2);
        imagesc(rateMap_norm(resort_index(:,2),:)); set(gca,'ydir','normal')
    end
    %hold on, vline(numSpatialBins(1)*cumsum(trackLengths)/sum(trackLengths),'w--'), hold off

elseif directionalDecoding==1

    rateMap_direction_1(setdiff(1:length(clusters),ind_cluster_direction_1),:) = 1;
    rateMap_direction_2(setdiff(1:length(clusters),ind_cluster_direction_2),:) = 1;

    [~, max_bin_map_1] = max(rateMap_direction_1,[],2);
    [~, max_bin_map_2] = max(rateMap_direction_2,[],2);
    resort_index_1 = sortrows([max_bin_map_1, (1:length(clusters))']);
    resort_index_1 = resort_index_1(isnan(resort_index_1(:,1))==0,:);
    resort_index_2 = sortrows([max_bin_map_2, (1:length(clusters))']);
    resort_index_2 = resort_index_2(isnan(resort_index_2(:,1))==0,:);
    rateMap_direction_1_norm = rateMap_direction_1./max(rateMap_direction_1,[],2);
    rateMap_direction_2_norm = rateMap_direction_2./max(rateMap_direction_2,[],2);

    if plot_figure==1
        clusters_to_remove = setdiff(1:length(clusters),[ind_cluster_direction_1; ind_cluster_direction_2])';
        clusters_for_plotting = setdiff(1:length(clusters),clusters_to_remove);

        figure()
        % for plotting, take out the clusters that didn't fire at all
        rateMap_direction_1_norm_plotting = rateMap_direction_1_norm;
        rateMap_direction_2_norm_plotting = rateMap_direction_2_norm;
        rateMap_direction_1_norm_plotting(prod(rateMap_direction_1_norm_plotting,2)==1,:) = 0;
        rateMap_direction_2_norm_plotting(prod(rateMap_direction_2_norm_plotting,2)==1,:) = 0;
        resort_index_1_plotting = resort_index_1(ismember(resort_index_1(:,2),clusters_for_plotting),:);
        resort_index_2_plotting = resort_index_2(ismember(resort_index_2(:,2),clusters_for_plotting),:);

        subplot(2,2,1), imagesc(rateMap_direction_1_norm_plotting(resort_index_1_plotting(:,2),:)), set(gca,'ydir','normal')
        ylabel('Cell num')
        xticks([0 25 50 75 100 125])
        xticklabels({'0','50','100','150','200','250'})
        subplot(2,2,2), imagesc(rateMap_direction_2_norm_plotting(resort_index_1_plotting(:,2),:)), set(gca,'ydir','normal')
        xticks([0 25 50 75 100 125])
        xticklabels({'0','50','100','150','200','250'})
        subplot(2,2,3), imagesc(rateMap_direction_1_norm_plotting(resort_index_2_plotting(:,2),:)), set(gca,'ydir','normal')
        xticks([0 25 50 75 100 125])
        xticklabels({'0','50','100','150','200','250'})
        ylabel('Cell num')
        xlabel('Track position (cm)')
        subplot(2,2,4), imagesc(rateMap_direction_2_norm_plotting(resort_index_2_plotting(:,2),:)), set(gca,'ydir','normal')
        xlabel('Track position (cm)')
        xticks([0 25 50 75 100 125])
        xticklabels({'0','50','100','150','200','250'})
        colormap parula

    end

end

if directionalDecoding == 1
    ind_cluster = {ind_cluster_direction_1 ind_cluster_direction_2};
    rateMap = [rateMap_direction_1 rateMap_direction_2];
    resort_index = {resort_index_1, resort_index_2};
end
