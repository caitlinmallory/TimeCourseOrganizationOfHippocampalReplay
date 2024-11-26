clearvars -except dayFiles day directory rat windows hand_clustered_only

plot_fig = 0;
%Load files
load Experiment_Information
%load Analysis_Information
load clusters
load Position_Data


spatialDim = Experiment_Information.spatialDim;

%for each spike of each cluster, interpolate behavior at spike time
Run_Times_concat = cat(1,Experiment_Information.Run_Times{:});
Run_Times_concat = [min(min(Run_Times_concat)) max(max(Run_Times_concat))];close 

for i = 1:length(clusters)
    [i,length(clusters)]
    
    spkTime = clusters(i).spkTime;
    if isempty(spkTime)
        continue
    end
    
    spkTimes_Run = compute_dataTemporalConcatenation(spkTime,Run_Times_concat); %eliminate spikes outside of run times
    ind = find(ismember(spkTime,spkTimes_Run));
    spkPos = NaN(length(spkTime),spatialDim);
    spkSpeed = NaN(length(spkTime),1);
    spkHD = NaN(length(spkTime),1);
    spkAccMag = NaN(length(spkTime),1);
    spkAccHD = NaN(length(spkTime),1);
 
    if size(Position_Data,2)<9
        Position_Data_sub = compute_dataInterpolation(Position_Data,spkTimes_Run,[4]);
    else
        Position_Data_sub = compute_dataInterpolation(Position_Data,spkTimes_Run,[4 9]);
    end
    
    if spatialDim==2
        spkPos(ind,:) = Position_Data_sub(:,2:3);
    elseif spatialDim==1
        spkPos(ind) = Position_Data_sub(:,2);
    end
    spkSpeed(ind) = Position_Data_sub(:,5);
    spkHD(ind) = Position_Data_sub(:,4);

    
    clusters(i).spkPos = spkPos;
    clusters(i).spkSpeed = spkSpeed;
    clusters(i).spkHD = spkHD;
    %clusters(i).spkAccMag = spkAccMag;
    %clusters(i).spkAccHD = spkAccHD;
    
    if plot_fig ==1
    if spatialDim==1
        %linear
        figure()
        plot(Position_Data(:,1),Position_Data(:,2))
        hold on, plot(spkTime,spkPos(:,1),'.'), hold off

    else
        figure()
        subplot(211)
        %open field
        plot(Position_Data(:,2),Position_Data(:,3))
        hold on, plot(spkPos(:,1),spkPos(:,2),'.'), hold off
        
        subplot(212)
        plot(Position_Data(:,1),Position_Data(:,2))
        hold on, plot(spkTime,spkPos(:,1),'.'), hold off
        drawnow
    end
    end
    
end
save('clusters.mat','clusters','-append')



