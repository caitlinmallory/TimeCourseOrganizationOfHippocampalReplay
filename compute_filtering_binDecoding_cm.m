function [x_NaN,x,timeBins,posteriorSpread,spikeDensity,positions,laser_state,timeBins_NaNremoved,posteriorPeak,jumps,timeBins_NaNremoved_merge] = compute_filtering_binDecoding_cm(binDecoding,sequence_posteriorSpreadThr,speedRangeThr,sequence_spikeDensityThr,positions,spikeDensity,spikeDensityStepSize,sequence_jumpThr,Times,binSize,spikeSampRate, laser_state)

    %decoded position
        x = binDecoding.posteriorCOM;
        timeBins = binDecoding.timeBins; 
        windowSizeDecoding = binDecoding.windowSizeDecoding;
        shiftSizeDecoding = binDecoding.shiftSizeDecoding;
        x = [mean(timeBins,2), x];
                
    %decoded position jumps

        jumps = [0; abs(diff(x(:,2)))];

    %posterior spread
        posteriorSpread = binDecoding.posteriorSpread;
        
    %posterior peak
        posteriorPeak = binDecoding.posteriorPeak;
                
    %spike density
        squareFilt = ones(1,round(windowSizeDecoding/spikeDensityStepSize))/round(windowSizeDecoding/spikeDensityStepSize); 
        gaussFilt = setUp_gaussFilt([1 500],windowSizeDecoding/spikeDensityStepSize);
        spikeDensity_squareFilt = zscore(conv(spikeDensity(:,2),squareFilt,'same'));
        spikeDensity_gaussFilt = zscore(conv(spikeDensity(:,2),gaussFilt,'same'));
        spikeDensity = [spikeDensity(:,1:2),spikeDensity_gaussFilt,spikeDensity_squareFilt];
        spikeDensity = compute_dataInterpolation(spikeDensity,x(:,1),[]);
        

    %laser state
        if ~isempty(laser_state)
            laser_state = compute_dataInterpolation(laser_state, x(:,1), []);
        else
            laser_state = [x(:,1) nan(size(x(:,1)))];
        end
                        
    %rat trajectory
        
        positions = compute_dataInterpolation(positions,x(:,1),[]);
       
    %cut to Time window
        data = compute_dataTemporalConcatenation([x,timeBins,jumps,posteriorSpread,posteriorPeak,spikeDensity,positions, laser_state],Times);
        
        x = data(:,1:3);
        timeBins = data(:,4:5);
        jumps = data(:,6);
        posteriorSpread = data(:,7);
        posteriorPeak = data(:,8);
        spikeDensity = data(:,9:12);
        positions = data(:,13:end-3);
        laser_state = data(:,end-1:end);
        clear data
        
    %filter for bins based spike density, posterior spread, and rat speed
        ind = find(positions(:,5)>=speedRangeThr(1) & positions(:,5)<speedRangeThr(2) & spikeDensity(:,3)>sequence_spikeDensityThr & posteriorSpread<(sequence_posteriorSpreadThr/binSize) & jumps<(sequence_jumpThr/binSize));
        
        if isempty(ind) 
            x_NaN = nan; x = nan; timeBins = nan; posteriorSpread = nan; spikeDensity = nan; positions = nan; laser_state = nan;...
                timeBins_NaNremoved = nan; posteriorPeak = nan; jump = nan; timeBins_NaNremoved_merge = nan;
        else
            
        x_NaN = x; 
        x_NaN(:,2:end) = NaN;
        x_NaN(ind,:) = x(ind,:);
        
        timeBins_NaNremoved = timeBins(ind,:);
        
        
        %merge overlapping time bins
        ind = [0; find(diff(timeBins_NaNremoved(:,1))>windowSizeDecoding*spikeSampRate);length(timeBins_NaNremoved(:,1))];            
        timeBins_NaNremoved_merge = [timeBins_NaNremoved(ind(1:end-1)+1,1) timeBins_NaNremoved(ind(2:end),2)];
        end
    %plots
%         subplot(121), histogram(spikeDensity(:,2),linspace(-2,6,50)), hold on, vline(spikeDensityThr,'k'), hold off
%         subplot(122), histogram(posteriorSpread,linspace(0,50,50)), hold on, vline(posteriorSpreadThr,'k'), hold off
%         keyboard
%         close
%         