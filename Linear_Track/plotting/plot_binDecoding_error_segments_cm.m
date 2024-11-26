clearvars -except dayFiles day directory rat windows hand_clustered_only

plot_figure = 0;

load Experiment_Information
load Analysis_Information
load binDecoding_error
load Position_Data

speedThr = 10;
spatialDim = Experiment_Information.spatialDim;

figure_folder = 'figures';
if ~isfolder(figure_folder)
    mkdir(figure_folder);
end
for sessionNum = 1:length(Experiment_Information.Segments)
    if ismember(12,Experiment_Information.Segments(sessionNum).Flags)
         Experiment_Information.Segments(sessionNum).decodingError = nan;
    continue
    end
    sessionNum_decoder = Experiment_Information.Segments(sessionNum).Decoder;

    timeBins = decoder_binDecoding(sessionNum_decoder).timeBins;
    posteriorCOM_sub = decoder_binDecoding(sessionNum_decoder).posteriorCOM;
    timePoints = mean(timeBins,2);

    ind_session = ismember(timePoints,compute_dataTemporalConcatenation(timePoints,Experiment_Information.Segments(sessionNum).Times));
    timePoints = timePoints(ind_session);
    posteriorCOM_sub = posteriorCOM_sub(ind_session,:);

    Position_Data_sub = compute_dataTemporalConcatenation(Position_Data,Experiment_Information.Segments(sessionNum).Times);
    Position_Data_sub = compute_locsToBins(Position_Data_sub,numSpatialBins,x_edges,y_edges);
    Position_Data_sub = compute_dataInterpolation(Position_Data_sub,timePoints,[]);

    ind_highSpeed = (Position_Data_sub(:,5)>speedThr);

    if spatialDim == 2
        error = sqrt((Position_Data_sub(:,2)-posteriorCOM_sub(:,1)).^2+(Position_Data_sub(:,3)-posteriorCOM_sub(:,2)).^2);
        map = nan(length(Position_Data_sub),1);

    elseif spatialDim == 1

        posteriorPeak = decoder_binDecoding(sessionNum_decoder).posteriorPeak;
        posteriorCOM_directional = nan(length(posteriorCOM_sub),1);
        [~,map] = max(posteriorPeak,[],2);
        for i = 1:length(map)
            posteriorCOM_directional(i) = posteriorCOM_sub(i,map(i));
        end
        posteriorCOM_sub = posteriorCOM_directional;

        error = abs(Position_Data_sub(:,2)-posteriorCOM_sub);
    end

    true_direction = nan(length(Position_Data_sub),1);
    true_direction(Position_Data_sub(:,6)<0) = 1;
    true_direction(Position_Data_sub(:,6)>0) = 2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ind = find(ind_highSpeed == 1);
    bad_ind = setdiff([1:length(Position_Data_sub)],ind);

    segment_error = error(ind);
    disp(['Segment' num2str(sessionNum) ' using decoder ' num2str(sessionNum_decoder)])
    disp(['Median error = ' num2str(nanmedian(segment_error),2)]);
    median_decoding_error = nanmedian(segment_error);

    map_copy(bad_ind) = NaN;
    true_direction_copy(bad_ind) = NaN;
    percent_correct_directional_assigment = sum(map_copy == true_direction_copy)/sum(~isnan(true_direction_copy));
%     disp(['Session ' num2str(sessionNum) ' using decoder ' num2str(sessionNum_decoder)])
%     disp(['% correct direction assigment = ' num2str(percent_correct_directional_assigment)]);

    Position_Data_copy = Position_Data_sub;
    posteriorCOM_copy = posteriorCOM_sub;
    Position_Data_copy(bad_ind,:) = NaN;
    posteriorCOM_copy(bad_ind,:) = NaN;

    if plot_figure ==1
        ha_fig();
        axes(ha(figure_panel));
        plot(Position_Data_copy(:,1),Position_Data_copy(:,2))
        hold on
        plot(Position_Data_copy(:,1),posteriorCOM_copy(:,1))
    end

    Experiment_Information.Segments(sessionNum).decodingError = median_decoding_error;
end

save('Experiment_Information','Experiment_Information')