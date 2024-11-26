function [numSpks] = load_numSpks_timeBins(timeBins,clusters,cluster_ind,shiftSizeDecoding,windowSizeDecoding)

%Collects spiking from all clusters within windows prescribed by timeBins
%TimeBins is unfiltered (continuous tiling of bins)


times = min(timeBins(:)):shiftSizeDecoding:max(timeBins(:));
numSpks = zeros(length(cluster_ind),size(timeBins,1));

%break up into subsets of clusters, since numSpks sometimes gets too big
cluster_ind_vec = [(1:10:length(cluster_ind))'; length(cluster_ind)];



for i = 1:length(cluster_ind_vec)-1
    %subset of clusters to deal with
        ind = cluster_ind_vec(i):cluster_ind_vec(i+1);
        cluster_ind_sub = cluster_ind(ind);
    
    %collects spikes in windows of size shiftSizeDecoding
        numSpks_micro_sub = zeros(length(cluster_ind_sub),length(times));
        for j = 1:length(cluster_ind_sub)
            numSpks_micro_sub(j,1:end-1) = histcounts(clusters(cluster_ind_sub(j)).spkTime,times);
        end

    %combines windows to collect spikes in window of size windowSizeDecoding
        numBins = round(windowSizeDecoding/shiftSizeDecoding);
        numSpks_sub = zeros(length(cluster_ind_sub),size(timeBins,1));
        for j = 1:size(timeBins,1)
            numSpks_sub(:,j) = nansum(numSpks_micro_sub(:,j:j+numBins-1),2);
        end


    %Caitlin changed the above code slightly, so that the spike count is
    %centered evenly on the windowSizeInQuestion. Pad the edges with extra
    %nans.
    
%         numBins = round(windowSizeDecoding/shiftSizeDecoding);
%         numSpks_sub = zeros(length(cluster_ind_sub),size(timeBins,1));
%         
%         num_nan_pad_bins = 5000;
%         numSpks_micro_sub_nan_padded = [nan(size(numSpks_micro_sub,1), num_nan_pad_bins) numSpks_micro_sub nan(size(numSpks_micro_sub,1), num_nan_pad_bins)];
%         
%         ind_start = num_nan_pad_bins+1;
%         ind_stop = num_nan_pad_bins + size(timeBins,1);
%        
%         for j = ind_start:ind_stop            
%             numSpks_sub(:,j-num_nan_pad_bins) = nansum(numSpks_micro_sub_nan_padded(:,j-(numBins-1)/2:j+(numBins-1)/2),2);
%         end   
    
        
        
    %add to larger numspk matrix
        numSpks(ind,:) = numSpks_sub;
    
    clear numSpks_micro_sub numSpks_sub    
end
    
    