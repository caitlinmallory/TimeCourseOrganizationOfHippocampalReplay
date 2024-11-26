function [Times_bounds] = load_timeBounds(Times)

if ~iscell(Times)
    Times_bounds = [min(Times(:)) max(Times(:))];
else
    Times_bounds = zeros(length(Times),2);
    for i = 1:length(Times)
        Times_bounds(i,:) = [min(Times{i}(:)) max(Times{i}(:))]; 
    end
end