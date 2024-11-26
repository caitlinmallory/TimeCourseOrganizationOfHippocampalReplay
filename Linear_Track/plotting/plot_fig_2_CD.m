%% Create a nice figure for publication-- OFF Data
sz = combined_hists.reverse_congruent_replays_spike;
sz = size(sz,2);
if align_to_stopping_period_start==1
    time_bins = 0.25:0.5:(0.5*sz-0.25);
    bins_to_plot = find(time_bins>=0 & time_bins<=10);
else
    time_bins = (-1*sz*0.5+0.25):0.5:-0.25;
    bins_to_plot = find(time_bins>= -10 & time_bins<=0);
end

hist_bin_centers_plot = time_bins;

% To plot ALL events (laser off and on):
num_laps_to_plot = 30;
events_to_plot = [{'reverse_congruent_replays_spike'},{'forward_congruent_replays_spike'}; {'reverse_congruent_replays_ripple'},{'forward_congruent_replays_ripple'}];...
    laser_state_to_plot = {'off'; 'off'};
replay_detection_type = {'spike';'ripple'};
figure_titles = {'spike, laser off'; 'ripple, laser off'};
for i = 1:size(events_to_plot,1)
    reverse_time_in_stopping_period_by_lap = nan(num_laps_to_plot,sz);
    forward_time_in_stopping_period_by_lap = nan(num_laps_to_plot,sz);
    for lap = 1:num_laps_to_plot
        inds_0 = find(combined_table.laser_state == 0 & combined_table.laser_state_pass_count == lap);
        inds_1 = find(combined_table.laser_state == 1 & combined_table.laser_state_pass_count == lap);
        inds_all = find(combined_table.laser_state_pass_count == lap);
        if strcmp(laser_state_to_plot{i}, 'off')
            inds = inds_0;
        elseif strcmp(laser_state_to_plot{i},'on')
            inds = inds_1;
        elseif strcmp(laser_state_to_plot{i},'both')
            inds = inds_all;
        elseif strcmp(laser_state_to_plot{i},'both_combined')
            inds = inds_all;
        end
        reverse_time_in_stopping_period_by_lap(lap,:) = nansum(combined_hists.(events_to_plot{i,1})(inds,:),1);
        forward_time_in_stopping_period_by_lap(lap,:) = nansum(combined_hists.(events_to_plot{i,2})(inds,:),1);
    end
    green_colormap = customcolormap(linspace(0,1,5), {'#014419','#1c7735','#5aae60','#a6db9d','#d7f1d6'});
    purple_colormap = customcolormap(linspace(0,1,5), {'#410149','#762a84','#9b6fac','#c1a5cd','#e7d4e8'});
    colorlimits = [0 2];

    fig1 = figure();
    fig1.Position = [600 600 300 250];
    tiledlayout(2,2);
    ax1=nexttile;
    imagesc(forward_time_in_stopping_period_by_lap);
    box off
    colormap(ax1,green_colormap)
    caxis(colorlimits)
    xticks (linspace(1,length(hist_bin_centers_plot),6))
    xticklabels({})
    xtickangle(0);
    yticks([1,10,20,num_laps_to_plot]);
    ylabel('Trial')
    ax1.FontSize = 8;
    xticks(0.5:4:sz+0.5)
    xlim([bins_to_plot(1)-0.5 bins_to_plot(end)+0.5])

    ax3=nexttile(3);
    imagesc(reverse_time_in_stopping_period_by_lap);
    box off
    colormap(ax3,purple_colormap)
    caxis(colorlimits)
    xticks(0.5:4:sz+0.5)
    xticklabels({})
    xtickangle(0);
    yticks([1,10,20,num_laps_to_plot]);
    ylabel('Trial')
    ax3.FontSize = 8;
    if align_to_reward_onset == 1 || align_to_stopping_period_start == 1
        xticks(0.5:4:sz+0.5)
        %       xticklabels(num2cell(-0.5*sz:2:0))
        xticklabels(num2cell(0:2:sz))
        xtickangle(0);
        xlabel('Time since arrival (s)')
        xlim([bins_to_plot(1)-0.5 bins_to_plot(end)+0.5])
    elseif align_to_stopping_period_departure==1||align_to_reward_offset==1
        xticks(0.5:4:sz+0.5)
        xticklabels(num2cell(-0.5*sz:2:0))
        xtickangle(0);
        xlabel('Time till departure (s)')
        xlim([bins_to_plot(1)-0.5 bins_to_plot(end)+0.5])
    end

    ax2 = nexttile(2);
    % Also plot forward rates and reverse rates on top of each other
    event_types = events_to_plot(i,1:2);

    mean_off_rev = mean_combined_hists_rates_off.(event_types{1});
    mean_on_rev = mean_combined_hists_rates_on.(event_types{1});
    mean_rev = mean_combined_hists_rates.(event_types{1});
    sem_off_rev = sem_combined_hists_rates_off.(event_types{1});
    sem_on_rev = sem_combined_hists_rates_on.(event_types{1});
    sem_rev = sem_combined_hists_rates.(event_types{1});

    mean_off_for = mean_combined_hists_rates_off.(event_types{2});
    mean_on_for = mean_combined_hists_rates_on.(event_types{2});
    mean_for = mean_combined_hists_rates.(event_types{2});
    sem_off_for = sem_combined_hists_rates_off.(event_types{2});
    sem_on_for = sem_combined_hists_rates_on.(event_types{2});
    sem_for = sem_combined_hists_rates.(event_types{2});

    colors = [.4660 0.6740 0.1880;0.4940 0.1840 0.5560];
    transparency_pcnt = 0.99;
    colors_transparent = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
        [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];

    if strcmp(laser_state_to_plot{i},'off')
        % Plot forward, off
        if align_to_stopping_period_departure==1 || align_to_reward_offset==1
            if smooth_rate_plots==1
                h = shadedErrorBar(hist_bin_centers_plot(1:end-(filter_length-1)/2),mean_off_for(1:end-(filter_length-1)/2),...
                    sem_off_for(1:end-(filter_length-1)/2),'lineprops','g'); hold on;
            else
                h = shadedErrorBar(hist_bin_centers_plot,mean_off_for,...
                    sem_off_for,'lineprops','g'); hold on;
            end
            xticks(-10:2:0)
            xlim([-10 0]);
        else
            h =  shadedErrorBar(hist_bin_centers_plot,mean_off_for,sem_off_for,'lineprops','g'); hold on
            xticks(0:2:10)
            xlim([0 10]);
        end
        h.patch.FaceColor = colors(1,:);
        h.mainLine.Color = colors(1,:);
        h.patch.FaceColor = colors_transparent(1,:);
        h.mainLine.LineWidth = 2;
        h.edge(1).Color = 'none';
        h.edge(2).Color = 'none';
    elseif strcmp(laser_state_to_plot{i},'on')
        % Plot forward, on
        if align_to_stopping_period_departure==1||align_to_reward_offset==1
            if smooth_rate_plots==1
                h = shadedErrorBar(hist_bin_centers_plot(1:end-(filter_length-1)/2),mean_on_for(1:end-(filter_length-1)/2),...
                    sem_on_for(1:end-(filter_length-1)/2),'lineprops','--g'); hold on;
            else
                h = shadedErrorBar(hist_bin_centers_plot,mean_on_for,...
                    sem_on_for,'lineprops','--g'); hold on;
            end
            xticks(-10:2:0)
            xlim([-10 0]);
        else
            h =  shadedErrorBar(hist_bin_centers_plot,mean_on_for,sem_on_for,'lineprops','--g'); hold on
            xticks(0:2:10)
            xlim([0 10]);
        end
        h.patch.FaceColor = colors(1,:);
        h.mainLine.Color = colors(1,:);
        h.patch.FaceColor = colors_transparent(1,:);
        h.mainLine.LineWidth = 2;
        h.edge(1).Color = 'none';
        h.edge(2).Color = 'none';
    elseif strcmp(laser_state_to_plot{i},'both')
        % Plot forward, both
        % Plot forward, on
        if align_to_stopping_period_departure==1||align_to_reward_offset==1
            if smooth_rate_plots==1
                h = shadedErrorBar(hist_bin_centers_plot(1:end-(filter_length-1)/2),mean_off_for(1:end-(filter_length-1)/2),...
                    sem_off_for(1:end-(filter_length-1)/2),'lineprops','-g'); hold on;
            else
                h = shadedErrorBar(hist_bin_centers_plot,mean_off_for,...
                    sem_off_for,'lineprops','-g'); hold on;
            end
            xticks(-10:2:0)
            xlim([-10 0]);
        else
            h =  shadedErrorBar(hist_bin_centers_plot,mean_off_for,sem_off_for,'lineprops','-g'); hold on
            xticks(0:2:10)
            xlim([0 10]);
        end
        h.mainLine.Color = colors(1,:);
        h.mainLine.LineWidth = 2;
        h.patch.FaceColor = colors_transparent(1,:);
        h.edge(1).Color = 'none';
        h.edge(2).Color = 'none';
        hold on
        % Plot forward, on
        if align_to_stopping_period_departure==1||align_to_reward_offset==1
            if smooth_rate_plots
                h = shadedErrorBar(hist_bin_centers_plot(1:end-(filter_length-1)/2),mean_on_for(1:end-(filter_length-1)/2),...
                    sem_on_for(1:end-(filter_length-1)/2),'lineprops','--g'); hold on;
            else
                h = shadedErrorBar(hist_bin_centers_plot,mean_on_for,...
                    sem_on_for,'lineprops','--g'); hold on;
            end
            xticks(-10:2:0)
            xlim([-10 0]);
        else
            h =  shadedErrorBar(hist_bin_centers_plot,mean_on_for,sem_on_for,'lineprops','--g'); hold on
            xticks(0:2:10)
            xlim([0 10]);
        end
        h.patch.FaceColor = colors(1,:);
        h.mainLine.Color = colors(1,:);
        h.mainLine.LineWidth = 2;
        h.patch.FaceColor = colors_transparent(1,:);
        h.edge(1).Color = 'none';
        h.edge(2).Color = 'none';
    elseif strcmp(laser_state_to_plot{i},'both_combined')
        % Plot forward, both
        if align_to_stopping_period_departure==1||align_to_reward_offset==1
            if smooth_rate_plots
                h = shadedErrorBar(hist_bin_centers_plot(1:end-(filter_length-1)/2),mean_for(1:end-(filter_length-1)/2),...
                    sem_for(1:end-(filter_length-1)/2),'lineprops','-g'); hold on;
            else
                h = shadedErrorBar(hist_bin_centers_plot,mean_for,...
                    sem_for,'lineprops','-g'); hold on;
            end
            xticks(-10:2:0)
            xlim([-10 0]);
        else
            h =  shadedErrorBar(hist_bin_centers_plot,mean_for,sem_for,'lineprops','-g'); hold on
            xticks(0:2:10)
            xlim([0 10]);
        end
        h.mainLine.Color = colors(1,:);
        h.mainLine.LineWidth = 2;
        h.patch.FaceColor = colors_transparent(1,:);
        h.edge(1).Color = 'none';
        h.edge(2).Color = 'none';
        hold on
    end

    if strcmp(laser_state_to_plot{i},'off')
        % Plot reverse, off
        if align_to_stopping_period_departure==1||align_to_reward_offset==1
            if smooth_rate_plots
                h = shadedErrorBar(hist_bin_centers_plot(1:end-(filter_length-1)/2),mean_off_rev(1:end-(filter_length-1)/2),...
                    sem_off_rev(1:end-(filter_length-1)/2),'lineprops','-r'); hold on;
            else
                h = shadedErrorBar(hist_bin_centers_plot,mean_off_rev,...
                    sem_off_rev,'lineprops','-r'); hold on;
            end
            xticks(-10:2:0)
            xlim([-10 0]);
        else
            h =  shadedErrorBar(hist_bin_centers_plot,mean_off_rev,sem_off_rev,'lineprops','-r'); hold on
            xticks(0:2:10)
            xlim([0 10]);
        end
        h.patch.FaceColor = colors(2,:);
        h.mainLine.Color = colors(2,:);
        h.mainLine.LineWidth = 2;
        h.patch.FaceColor = colors_transparent(2,:);
        h.edge(1).Color = 'none';
        h.edge(2).Color = 'none';
    elseif strcmp(laser_state_to_plot{i},'on')
        % Plot reverse on
        if align_to_stopping_period_departure==1||align_to_reward_offset==1
            if smooth_rate_plots==1
                h = shadedErrorBar(hist_bin_centers_plot(1:end-(filter_length-1)/2),mean_on_rev(1:end-(filter_length-1)/2),...
                    sem_on_rev(1:end-(filter_length-1)/2),'lineprops','--r'); hold on;
            else
                h = shadedErrorBar(hist_bin_centers_plot,mean_on_rev,...
                    sem_on_rev,'lineprops','--r'); hold on;
            end
            xticks(-10:2:0)
            xlim([-10 0]);
        else
            h =  shadedErrorBar(hist_bin_centers_plot,mean_on_rev,sem_on_rev,'lineprops','--r'); hold on
            xticks(0:2:10)
            xlim([0 10]);
        end
        h.patch.FaceColor = colors(2,:);
        h.mainLine.Color = colors(2,:);
        h.mainLine.LineWidth = 2;
        h.patch.FaceColor = colors_transparent(2,:);
        h.edge(1).Color = 'none';
        h.edge(2).Color = 'none';
    elseif strcmp(laser_state_to_plot{i},'both')
        % Plot reverse, off and reverse, on
        if align_to_stopping_period_departure==1||align_to_reward_offset==1
            if smooth_rate_plots==1
                h = shadedErrorBar(hist_bin_centers_plot(1:end-(filter_length-1)/2),mean_off_rev(1:end-(filter_length-1)/2),...
                    sem_off_rev(1:end-(filter_length-1)/2),'lineprops','-r'); hold on;
            else
                h = shadedErrorBar(hist_bin_centers_plot,mean_off_rev,...
                    sem_off_rev,'lineprops','-r'); hold on;
            end
            xticks(-10:2:0)
            xlim([-10 0]);
        else
            h =  shadedErrorBar(hist_bin_centers_plot,mean_off_rev,sem_off_rev,'lineprops','-r'); hold on
            xticks(0:2:10)
            xlim([0 10]);
        end
        h.patch.FaceColor = colors(2,:);
        h.mainLine.Color = colors(2,:);
        h.mainLine.LineWidth = 2;
        h.patch.FaceColor = colors_transparent(2,:);
        h.edge(1).Color = 'none';
        h.edge(2).Color = 'none';
        hold on
        if align_to_stopping_period_departure==1||align_to_reward_offset==1
            if smooth_rate_plots==1
                h = shadedErrorBar(hist_bin_centers_plot(1:end-(filter_length-1)/2),mean_on_rev(1:end-(filter_length-1)/2),...
                    sem_on_rev(1:end-(filter_length-1)/2),'lineprops','--r'); hold on;
            else
                h = shadedErrorBar(hist_bin_centers_plot,mean_on_rev,...
                    sem_on_rev,'lineprops','--r'); hold on;
            end
            xticks(-10:2:0)
            xlim([-10 0]);
        else
            h =  shadedErrorBar(hist_bin_centers_plot,mean_on_rev,sem_on_rev,'lineprops','--r'); hold on
            xticks(0:2:10)
            xlim([0 10]);
        end
        h.patch.FaceColor = colors(2,:);
        h.mainLine.Color = colors(2,:);
        h.mainLine.LineWidth = 2;
        h.patch.FaceColor = colors_transparent(2,:);
        h.edge(1).Color = 'none';
        h.edge(2).Color = 'none';
    elseif strcmp(laser_state_to_plot{i},'both_combined')
        % Plot reverse, off and reverse, on
        if align_to_stopping_period_departure==1||align_to_reward_offset==1
            if smooth_rate_plots==1
                h = shadedErrorBar(hist_bin_centers_plot(1:end-(filter_length-1)/2),mean_rev(1:end-(filter_length-1)/2),...
                    sem_rev(1:end-(filter_length-1)/2),'lineprops','-r'); hold on;
            else
                h = shadedErrorBar(hist_bin_centers_plot,mean_rev,...
                    sem_rev,'lineprops','-r'); hold on;
            end
            xticks(-10:2:0)
            xlim([-10 0]);
        else
            h =  shadedErrorBar(hist_bin_centers_plot,mean_rev,sem_rev,'lineprops','-r'); hold on
            xticks(0:2:10)
            xlim([0 10]);
        end
        h.patch.FaceColor = colors(2,:);
        h.mainLine.Color = colors(2,:);
        h.mainLine.LineWidth = 2;
        h.patch.FaceColor = colors_transparent(2,:);
        h.edge(1).Color = 'none';
        h.edge(2).Color = 'none';
        hold on
    end
    xticklabels({})
    ylabel('Events/s')
    hold on
    ylim(fig_ylim)
    xticklabels({})
    xtickangle(0);

    ax4 = nexttile(4);
    if strcmp(replay_detection_type{i},'spike')
        event_type = 'congruent_forward_minus_reverse_spike';
    elseif strcmp(replay_detection_type{i},'ripple')
        event_type = 'congruent_forward_minus_reverse_ripple';
    end
    colors = [0 0 0; 0 0 0];
    transparency_pcnt = 0.99;
    colors_transparent = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
        [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];


    mean_off = mean_combined_hists_rates_off.(event_type)(:);
    mean_on = mean_combined_hists_rates_on.(event_type)(:);
    mean_combined = mean_combined_hists_rates.(event_type)(:);
    sem_off = sem_combined_hists_rates_off.(event_type)(:);
    sem_on = sem_combined_hists_rates_on.(event_type)(:);
    sem_combined = sem_combined_hists_rates.(event_type)(:);

    if strcmp(laser_state_to_plot{i},'off')==1
        if align_to_stopping_period_departure==1||align_to_reward_offset==1
            if smooth_rate_plots==1
                h = shadedErrorBar(hist_bin_centers_plot(1:end-(filter_length-1)/2),mean_off(1:end-(filter_length-1)/2),...
                    sem_off(1:end-(filter_length-1)/2),'lineprops','-k'); hold on;
            else
                h = shadedErrorBar(hist_bin_centers_plot,mean_off,...
                    sem_off,'lineprops','-k'); hold on;
            end
            xticks(-10:2:0)
            xlim([-10 0]);
        else
            h =  shadedErrorBar(hist_bin_centers_plot,mean_off,sem_off,'lineprops','-k'); hold on
            xticks(0:2:10)
            xlim([0 10]);
        end
        h.mainLine.LineWidth = 2;
        h.edge(1).Color = 'none';
        h.edge(2).Color = 'none';
        hold on
        yline(0)
    elseif strcmp(laser_state_to_plot{i},'on')==1

        if align_to_stopping_period_departure==1||align_to_reward_offset==1
            if smooth_rate_plots==1
                h = shadedErrorBar(hist_bin_centers_plot(1:end-(filter_length-1)/2),mean_on(1:end-(filter_length-1)/2),...
                    sem_on(1:end-(filter_length-1)/2),'lineprops','--k'); hold on;
            else
                h = shadedErrorBar(hist_bin_centers_plot,mean_on,...
                    sem_on,'lineprops','--k'); hold on;
            end
            xticks(-10:2:0)
            xlim([-10 0]);
        else
            h =  shadedErrorBar(hist_bin_centers_plot,mean_on,sem_on,'lineprops','--k'); hold on
            xticks(0:2:10)
            xlim([0 10]);
        end
        h.mainLine.LineWidth = 2;
        h.edge(1).Color = 'none';
        h.edge(2).Color = 'none';
        hold on
        yline(0)
    elseif strcmp(laser_state_to_plot{i},'both')==1
        if align_to_stopping_period_departure==1||align_to_reward_offset==1
            if smooth_rate_plots==1
                h = shadedErrorBar(hist_bin_centers_plot(1:end-(filter_length-1)/2),mean_off(1:end-(filter_length-1)/2),...
                    sem_off(1:end-(filter_length-1)/2),'lineprops','-k'); hold on;
            else
                h = shadedErrorBar(hist_bin_centers_plot,mean_off,...
                    sem_off,'lineprops','-k'); hold on;
            end
            xticks(-10:2:0)
            xlim([-10 0]);
        else
            h =  shadedErrorBar(hist_bin_centers_plot,mean_off,sem_off,'lineprops','-k'); hold on
            xticks(0:2:10)
            xlim([0 10]);
        end
        h.mainLine.LineWidth = 2;
        h.edge(1).Color = 'none';
        h.edge(2).Color = 'none';
        hold on
        if align_to_stopping_period_departure==1||align_to_reward_offset==1
            if smooth_rate_plots == 1
                h = shadedErrorBar(hist_bin_centers_plot(1:end-(filter_length-1)/2),mean_on(1:end-(filter_length-1)/2),...
                    sem_on(1:end-(filter_length-1)/2),'lineprops','--k'); hold on;
            else
                h = shadedErrorBar(hist_bin_centers_plot,mean_on,...
                    sem_on,'lineprops','--k'); hold on;
            end
            xticks(-10:2:0)
            xlim([-10 0]);
        else
            h =  shadedErrorBar(hist_bin_centers_plot,mean_on,sem_on,'lineprops','--k'); hold on
            xticks(0:2:10)
            xlim([0 10]);
        end
        h.mainLine.LineWidth = 2;
        h.edge(1).Color = 'none';
        h.edge(2).Color = 'none';
        hold on
        yline(0)
    elseif strcmp(laser_state_to_plot{i},'both_combined')==1
        if align_to_stopping_period_departure==1||align_to_reward_offset==1
            if smooth_rate_plots==1
                h = shadedErrorBar(hist_bin_centers_plot(1:end-(filter_length-1)/2),mean_combined(1:end-(filter_length-1)/2),...
                    sem_combined(1:end-(filter_length-1)/2),'lineprops','-k'); hold on;
            else
                h = shadedErrorBar(hist_bin_centers_plot,mean_combined,...
                    sem_combined,'lineprops','-k'); hold on;
            end
            xticks(-10:2:0)
            xlim([-10 0]);
        else
            h =  shadedErrorBar(hist_bin_centers_plot,mean_combined,sem_combined,'lineprops','-k'); hold on
            xticks(0:2:10)
            xlim([0 10]);
        end
        h.mainLine.LineWidth = 2;
        h.edge(1).Color = 'none';
        h.edge(2).Color = 'none';
        hold on
        yline(0)
    end

    ylabel('Events/s')
    ylim(fig_ylim_diff)
    xtickangle(0);
    y_limit = gca().YLim;
    ylim([y_limit(1) (y_limit(2) + 0.01)]);
    y_limit = gca().YLim;
    y_pos_for_astrisks = y_limit(2);
    hold on
    if strcmp(laser_state_to_plot{i},'off')==1
        if strcmp(replay_detection_type{i},'spike')==1
            plot(hist_bin_centers_plot(congruent_for_minus_rev_spike_off_sig ==1),repmat(y_pos_for_astrisks,[sum(congruent_for_minus_rev_spike_off_sig ) 1]), '.k')
        elseif strcmp(replay_detection_type{i},'ripple')==1
            plot(hist_bin_centers_plot(congruent_for_minus_rev_ripple_off_sig ==1),repmat(y_pos_for_astrisks,[sum(congruent_for_minus_rev_ripple_off_sig ) 1]), '.k')
        end
    elseif strcmp(laser_state_to_plot{i},'on')==1
        if strcmp(replay_detection_type{i},'spike')==1
            plot(hist_bin_centers_plot(congruent_for_minus_rev_spike_on_sig ==1),repmat(y_pos_for_astrisks,[sum(congruent_for_minus_rev_spike_on_sig) 1]), '.k')
        elseif strcmp(replay_detection_type{i},'ripple')==1
            plot(hist_bin_centers_plot(congruent_for_minus_rev_ripple_on_sig ==1),repmat(y_pos_for_astrisks,[sum(congruent_for_minus_rev_ripple_on_sig) 1]), '.k')
        end
    elseif (strcmp(laser_state_to_plot{i},'both')==1 || strcmp(laser_state_to_plot{i},'both_combined'))
        if strcmp(replay_detection_type{i},'spike')==1
            plot(hist_bin_centers_plot(congruent_for_minus_rev_spike_sig ==1),repmat(y_pos_for_astrisks,[sum(congruent_for_minus_rev_spike_sig) 1]), '.k')
        elseif strcmp(replay_detection_type{i},'ripple')==1
            plot(hist_bin_centers_plot(congruent_for_minus_rev_ripple_sig ==1),repmat(y_pos_for_astrisks,[sum(congruent_for_minus_rev_ripple_sig) 1]), '.k')
        end
    end
    ax4.FontSize = 8;
    set(gcf, 'Color', 'white','Renderer','painters');
    set(gcf, 'PaperPositionMode', 'auto');
    if align_to_reward_onset==1 || align_to_stopping_period_start==1
        set(gcf, 'PaperPositionMode', 'auto');
        fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures';
        saveas(gcf,fullfile(fig_path,['Fig1_aligned_to_reward_' replay_detection_type{i} '_laser' laser_state_to_plot{i}]),'pdf')
        saveas(gcf,fullfile(fig_path,['Fig1_aligned_to_reward_' replay_detection_type{i} '_laser' laser_state_to_plot{i}]),'jpg')
    else
        set(gcf, 'PaperPositionMode', 'auto');
        fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures';
        saveas(gcf,fullfile(fig_path,['Fig1_aligned_to_departure_' replay_detection_type{i} '_laser' laser_state_to_plot{i}]),'pdf')
        saveas(gcf,fullfile(fig_path,['Fig1_aligned_to_departure_' replay_detection_type{i} '_laser' laser_state_to_plot{i}]),'jpg')
    end
end
