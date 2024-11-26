function [timeBins,centers] = load_timeBins(Times,windowShift,windowSize)

 %time bins
        startTimes = (min(Times(:)):windowShift:(max(Times(:))-windowSize))';
        endTimes = startTimes+windowSize;
        centers = mean([startTimes endTimes],2)';

    %eliminate bins outside times
%         ind = [];
%         for i = 1:size(Times,1)
%             ind = [ind; find(middleTimes>Times(i,1) & middleTimes<Times(i,2))];
%         end
%         startTimes = startTimes(ind);
%         endTimes = endTimes(ind);

    %convert to cell (ith cell contains time windows that contribute to ith step)      
%         timeBins = cell(length(startTimes),1);
%         for i = 1:length(startTimes)
%             timeBins{i} = [startTimes(i) endTimes(i)];
%         end

timeBins = [startTimes endTimes];
