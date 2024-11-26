clearvars -except dayFiles day directory rat windows hand_clustered_only

plot_figure = 0;

load Experiment_Information
load Analysis_Information
load binDecoding_error
load Position_Data

speedThr = 10;
spatialDim = Experiment_Information.spatialDim;
directionalDecoding = 1;
Run_Times = Experiment_Information.Run_Times;
numSessions = length(Run_Times);
edges = linspace(0,30,100);

if plot_figure == 1
    ha_fig = figure('Position', [50, 50, 1500, 860]);
    ha = tight_subplot(numSessions,numSessions,[.02 .02],[.02 .02],[.02 .02]);
    ha_fig2 = figure('Position', [50, 50, 1500, 860]);
    ha2 = tight_subplot(numSessions,numSessions,[.02 .02],[.02 .02],[.02 .02]);
end

% if isfield(Experiment_Information,'homeLoc')
%     if sum([Experiment_Information.homeLoc{1}]) > 0
%         homeLocs = nan(length(Experiment_Information.Segments),2);
%         for run = 1:length(Experiment_Information.Segments)
%             homeLocs(run,1) = Experiment_Information.homeLoc{run}(1,1);
%             homeLocs(run,2) = Experiment_Information.homeLoc{run}(1,2);
%         end
% 
%         homeLocs = compute_locsToBins(homeLocs,numSpatialBins,x_edges,y_edges);
%         homeLocs = unique(homeLocs,'rows','stable');
% 
%         th = 0:pi/50:2*pi;
%         xunit = nan(length(homeLocs),length(th));
%         yunit = nan(length(homeLocs),length(th));
%         for home = 1:length(homeLocs)
%             xunit(home,:) = 15/binSize*cos(th) + homeLocs(home,1);
%             yunit(home,:) = 15/binSize*sin(th) + homeLocs(home,2);
%         end
% 
%         homeLoc_colors = cool(length(homeLocs));
% 
%     end
% end

figure_folder = 'figures';
if ~isfolder(figure_folder)
    mkdir(figure_folder);
end
for sessionNum_decoder = 1:length(Run_Times)

    timeBins = decoder_binDecoding(sessionNum_decoder).timeBins;
    times = mean(timeBins,2);

    Position_Data_scaled = compute_locsToBins(Position_Data,numSpatialBins,x_edges,y_edges);
    Position_Data_sub = compute_dataInterpolation(Position_Data_scaled,times,[]);

    ind_highSpeed = (Position_Data_sub(:,5)>speedThr);

    posteriorCOM = decoder_binDecoding(sessionNum_decoder).posteriorCOM;

    if spatialDim == 2
        error = sqrt((Position_Data_sub(:,2)-posteriorCOM(:,1)).^2+(Position_Data_sub(:,3)-posteriorCOM(:,2)).^2);
        map = nan(length(Position_Data_sub),1);

    elseif spatialDim == 1

        posteriorPeak = decoder_binDecoding(sessionNum_decoder).posteriorPeak;
        posteriorCOM_directional = nan(length(posteriorCOM),1);
        [~,map] = max(posteriorPeak,[],2);
        for i = 1:length(map)
            posteriorCOM_directional(i) = posteriorCOM(i,map(i));
        end
        posteriorCOM = posteriorCOM_directional;

        error = abs(Position_Data_sub(:,2)-posteriorCOM);
    end

    true_direction = nan(length(Position_Data_sub),1);
    true_direction(Position_Data_sub(:,6)<0) = 1;
    true_direction(Position_Data_sub(:,6)>0) = 2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a = lines(length(Run_Times));
    for sessionNum = 1:length(Run_Times)

        figure_panel = sub2ind([length(Run_Times) length(Run_Times)],sessionNum_decoder,sessionNum);


        Position_Data_copy = Position_Data_sub;
        true_direction_copy = true_direction;
        posteriorCOM_copy = posteriorCOM;
        map_copy = map;


        ind_session = ismember(times,compute_dataTemporalConcatenation(times,load_timeBounds(Run_Times{sessionNum})),'rows');
        ind = find(ind_highSpeed == 1 & ind_session == 1);
        bad_ind = setdiff([1:length(Position_Data_sub)],ind);

        session_error = error(ind);
        disp(['Session ' num2str(sessionNum) ' using decoder ' num2str(sessionNum_decoder)])
        disp(['Median error = ' num2str(nanmedian(session_error),2)]);
        median_decoding_error = nanmedian(session_error);

        map_copy(bad_ind) = NaN;
        true_direction_copy(bad_ind) = NaN;
        percent_correct_directional_assigment = sum(map_copy == true_direction_copy)/sum(~isnan(true_direction_copy));
        disp(['Session ' num2str(sessionNum) ' using decoder ' num2str(sessionNum_decoder)])
        disp(['% correct direction assigment = ' num2str(percent_correct_directional_assigment)]);

        Position_Data_copy(bad_ind,:) = NaN;
        posteriorCOM_copy(bad_ind,:) = NaN;

        if plot_figure ==1
        ha_fig();
        axes(ha(figure_panel));
        plot(Position_Data_copy(:,1),Position_Data_copy(:,2))
        hold on
        plot(Position_Data_copy(:,1),posteriorCOM_copy(:,1))
        end
        %mean decoding error across decoders
        %         figure(), set(gcf,'color','w')
        %         subplot(1,length(Run_Times),sessionNum)
        %         hold on, bar(sessionNum_decoder,nanmean(error(ind)),'facecolor',a(sessionNum_decoder,:),'linewidth',1)
        %         hold on, errorbar(sessionNum_decoder,nanmean(error(ind)),nanstd(error(ind))/sqrt(length(ind)),'linestyle','--')
        %         if sessionNum_decoder==1
        %             axis square
        %             set(gcf,'color','w')
        %             set(gca,'xtick',1:length(Run_Times),'xticklabel',1:length(Run_Times))
        %             title(strcat('session ',num2str(sessionNum)))
        %             xlabel('decoder')
        %             ylabel('error (bins)')
        %             ylim([0 60])
        %         end

        %         %shuffles
        %         numShuffles = length(decoder_binDecoding(sessionNum_decoder).shuffle);
        %         shuffle_mean = zeros(numShuffles,1);
        %         for i = 1:numShuffles
        %             posteriorCOM_shuffle = decoder_binDecoding(sessionNum_decoder).shuffle(i).posteriorCOM;
        %             error_shuffle = sqrt((Position_Data_sub(:,2)-posteriorCOM_shuffle(:,1)).^2+(Position_Data_sub(:,3)-posteriorCOM_shuffle(:,2)).^2);
        %             shuffle_mean(i) = nanmean(error_shuffle(ind));
        %         end
        %         bar(sessionNum_decoder,mean(shuffle_mean),'facecolor','none','linestyle','--')
    end



    
end

if spatialDim == 2
    %spatial histograms of decoding error to see how error distributes in space
    %         figure(),
    %         subplot(length(Run_Times),length(Run_Times),(sessionNum-1)*length(Run_Times)+sessionNum_decoder),

    [M,N,SEM] = compute_positionHist(Position_Data_sub(ind,2:3),error(ind),numSpatialBins,0);
    if plot_figure==1
        figure(ha_fig2)
        axes(ha2(figure_panel));


        imagesc(M,'alphadata',~isnan(M)),
        colormap parula
        %hold on, plot_barrierLocs(barrierLocs{sessionNum_decoder},1,1,numSpatialBins,4), hold off
        set(gcf,'color','w')
        set(gca,'ydir','normal','ytick',[],'xtick',[]), axis square, caxis([0 max(edges)])
        hold on
        %
        if isfield(Experiment_Information,'homeLoc')
            if sum([Experiment_Information.homeLoc{1}]) > 0
                for home = 1:length(homeLocs)
                    plot(xunit(home,:),yunit(home,:),'color',homeLoc_colors(home,:),'Linewidth',3)
                end
            end
        end
        title(strcat('S',num2str(sessionNum),'_D',num2str(sessionNum_decoder),'_',num2str(nanmedian(session_error),2)),'Interpreter', 'none')
    end
end


%cumsum histograms of decoding error
%         figure(3),
%         subplot(1,length(Run_Times),sessionNum)
%         h = histc(error(ind),edges);
%         hold on, plot(edges,cumsum(h)/sum(h),'linewidth',2,'color',a(sessionNum,:)), hold off




if plot_figure==1
    saveas(ha_fig2,[pwd '\figures\decodingError'],'tiff')
end
if exist('session_wide_properties.mat')==2
    save('session_wide_properties','median_decoding_error','percent_correct_directional_assigment','-append')
else
    save('session_wide_properties','median_decoding_error','percent_correct_directional_assigment')
end