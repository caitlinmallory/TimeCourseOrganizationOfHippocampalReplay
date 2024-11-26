function [Data_concat,Data_NaN] = compute_dataTemporalConcatenation(Data,Times)


% If Times is just 2 values, the code below is much faster:
if size(Times,1)==1
    Data_concat=Data(Data(:,1)>=Times(1) & Data(:,1)<Times(2),:);

    % Data_NaN returns a matrix the same size as the input data. Data
    % outside the requested times will be NaNed out. This is slow, so only
    % do it if requested.
    if nargout > 1
        Data_NaN = Data; Data_NaN(:,2:end) = NaN;
        ind = find(Data(:,1)>=Times(1) & Data(:,1)<Times(2));
        Data_NaN(ind,:) = Data(ind,:);
    end
else

    Data_concat = [];
    if iscell(Times)==1
        for i = 1:length(Times)
            for j = 1:size(Times{i})
                Data_concat = [Data_concat; Data((Data(:,1)>=Times{i}(j,1) & Data(:,1)<Times{i}(j,2)),:)];
            end
        end
    else
        for i = 1:size(Times,1)
            Data_concat = [Data_concat; Data((Data(:,1)>=Times(i,1) & Data(:,1)<Times(i,2)),:)];
        end
    end

    % Data_NaN returns a matrix the same size as the input data. Data
    % outside the requested times will be NaNed out. This is slow, so only
    % do it if requested.
    if nargout > 1
        Data_NaN = Data; Data_NaN(:,2:end) = NaN;
        if iscell(Times)==1
            for i = 1:length(Times)
                for j = 1:size(Times{i})
                    ind = find(Data(:,1)>=Times{i}(j,1) & Data(:,1)<Times{i}(j,2));
                    Data_NaN(ind,:) = Data(ind,:);
                end
            end
        else
            for i = 1:size(Times,1)
                ind = find(Data(:,1)>=Times(i,1) & Data(:,1)<Times(i,2));
                Data_NaN(ind,:) = Data(ind,:);
            end
        end
    end
end
% Data_concat = unique(Data_concat,'rows');
% Data_NaN = unique(Data_NaN,'rows');