function [positions_full] = load_positions_full(Run_Times,Sleep_Times,times,positions,circular_columns)


positions = compute_dataInterpolation(positions,mean(times,2),circular_columns);
positions_full = positions; positions_full(:,2:end) = NaN;

for i = 1:length(Run_Times)

    [positions_sub,positions_NaN_sub] = compute_dataTemporalConcatenation(positions,Run_Times{i});
 
    [~,ind] = ismember(positions_sub(:,1),positions_NaN_sub);
    
    positions_full(ind,2:end) = positions_sub(:,2:end);
end

if ~isempty(Sleep_Times)
    for i = 1:length(Sleep_Times)
        [positions_sub,positions_NaN_sub] = compute_dataTemporalConcatenation(positions,Sleep_Times{i});

        [~,ind] = ismember(positions_sub(:,1),positions_NaN_sub);
        positions_full(ind,2:end) = 0;
    end
end


% figure()
% subplot(211), plot(positions(:,1),positions(:,4),'.-')
% subplot(212), plot(positions_full(:,1),positions_full(:,4),'.-')
