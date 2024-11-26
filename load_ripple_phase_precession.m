
clearvars -except dayFiles day directory rat windows hand_clustered_only

fig_folder = 'ripple_phase_precession';
if exist(fig_folder) == 0
    mkdir(fig_folder)
end



plot_fig = 0;
analyze_time_versus_phase = 0;
analyze_distance_versus_phase = 1;
% This script adds a table containing 9 entries (ripple-locking in all
% phase locking within all candidate events (ripples)
% phase locking within all laser off candidate events
% phase locking within all laser on candidate events
% phase locking within all laser off forward replays events
% phase locking within all laser off reverse replays events
% phase locking within all laser on forward replays events
% phase locking within all laser on reverse replays events
% phase locking within all forward replays
% phase locking within all reverse replays

%Each entry has 5 columns:
% 1. The pval from a Rayleigh test- tells whether the cell is modulated by
% ripple phase
% 2. The mean resultant vector
% 3. The number of spikes
% 4. The mean phase at which the cell fires
% 5. pval for a shuffle (not currently done because it takes a very long
% time).
min_time_into_stopping_period = 0;
max_time_into_stopping_period = inf;
coverage_thr = 0.0; % for replays
weighted_r_thr = 0.6; % for replays
posterior_diff_thr = 0.33;
do_shuffle = 0;
phase_interpolation_method = 1; % 1 = fast, but less accurate. 2 = slow, much more accurate.

if do_shuffle==1
    % load the ripple phase for all timepoints sampled in the session
    load zscored_lfp_power.mat
    session_ripple_phase = [zscored_lfp_power.Ripple(:,1) zscored_lfp_power.Ripple(:,4)];
end

load zscored_ripple_power.mat
load clusters
load Experiment_Information
spikeSampRate = Experiment_Information.spikeSampRate;

% get the times associated with all forward, reverse replays, laser on and
% laser off.

% for each cluster, calculate the degree of ripple phase locking


% phase locking within all candidate events (ripples)
% phase locking within all laser off candidate events
% phase locking within all laser on candidate events

% phase locking within all laser off forward replays events
% phase locking within all laser off reverse replays events
% phase locking within all laser on forward replays events
% phase locking within all laser on reverse replays events
times_to_look_at_phase_locking = cell(10,1);

% Find your ripples
replay_selection = 'ripple_filtered';
sde_thr = -inf;
ripple_thr = 3;
use_duration_og = 1; % set to 1 if looking at ripples or sde's. set to 0 if looking at replays.
t_ripples = load_replays_from_individual_session(replay_selection,use_duration_og,-inf,-inf,posterior_diff_thr,sde_thr,ripple_thr,min_time_into_stopping_period,max_time_into_stopping_period);

% Find your replays
replay_selection = 'spike_filtered';
sde_thr = 3;
ripple_thr = -inf;
use_duration_og = 0; % set to 1 if looking at ripples or sde's. set to 0 if looking at replays.
t = load_replays_from_individual_session(replay_selection,use_duration_og,coverage_thr,weighted_r_thr,posterior_diff_thr,sde_thr,ripple_thr,min_time_into_stopping_period,max_time_into_stopping_period);
t_replay = t(t.replay==1,:);

%% Phase precession analyses:

if analyze_time_versus_phase == 1
% First, look at phase precession as a function of normalized time in event
% versus spike phase (this is done on a cell by cell basis, regardless of
% whether the cell has a field).

a_all_fractional_time_in_event = [];
a_all_ripple_phase_in_event = [];
a_all_unwrapped_ripple_phase_in_event = [];
a_all_event_number = [];
a_all_event_laser_state = [];
a_all_event_direction = [];
a_all_event_best_map = [];
a_all_event_is_replay = [];
a_all_event_cluster_id = [];

% Plot this to make unwrapping clear
% a = pi + (0:20).*(0.9*2*pi)
% figure()
% plot(wrapTo2Pi(a))
% plot(unwrap(a)) 

for i = 1:length(clusters)


    % In progress: look at phase precession to ripples for each cell.
    spkTimes = clusters(i).spkTime;
    spkLaserStatus = clusters(i).spkLaserState;
    spkRipplePhase = clusters(i).spkRipplePhase;
    spkPos = clusters(i).spkPos(:,1);
    spkData = [spkTimes spkLaserStatus spkRipplePhase spkPos];

    % Todo: loop over candidate events (SD, ripple, forward, reverse, etc)
%     sub_replays = t;
    sub_replays = t_replay;

    for sub_event_num = 1:height(sub_replays)

        % for each sub_replay, you need the time and phase of each spike
        t_replay_og = sub_replays.timePoints_og(sub_event_num ,:);

        % Limit spikes to this CE. Only consider events that had at least 3
        % spikes
        spkData_in_event = compute_dataTemporalConcatenation(spkData,t_replay_og);
        if size(spkData_in_event,1)<3
            continue
        end

            % Look at the ripple during this time
     
         time_start = spkData_in_event(1,1);
         time_end = spkData_in_event(end,1);
         fractional_time_in_event = (spkData_in_event(:,1)-time_start)./(time_end-time_start);
         lfp_during_replay = compute_dataTemporalConcatenation(zscored_ripple_power,t_replay_og);
         ripple_at_spikes = compute_dataInterpolation(lfp_during_replay,spkData_in_event(:,1),[4]);
%          ax1 = subplot(3,1,1);
%          plot(fractional_time_in_event,spkData_in_event(:,3),'ok','MarkerFaceColor','k')
%          ax2 = subplot(3,1,2);
%          plot(fractional_time_in_event,unwrap(spkData_in_event(:,3)),'ok','MarkerFaceColor','k')
%          ax3 = subplot(3,1,3);
%          plot((lfp_during_replay(:,1)-time_start)./(time_end-time_start),lfp_during_replay(:,3))
%          hold on
%          plot((ripple_at_spikes(:,1)-time_start)./(time_end-time_start),ripple_at_spikes(:,3),'.k','MarkerSize',10);
          
         numTimeBins = size(t_replay.posterior_full_left{sub_event_num},1);
         ripple_pseudo_time = linspace(1,numTimeBins,size(lfp_during_replay,1));
         spikes_pseudo_time = compute_dataInterpolation([lfp_during_replay(:,1) ripple_pseudo_time'],spkData_in_event(:,1),[]);
         spikes_pseudo_time = spikes_pseudo_time(:,2);


        figure('Position',[473 1 1025 943])

        ax1 = subplot(3,2,1);
        plot(ripple_pseudo_time,lfp_during_replay(:,3))
        hold on
        plot(spikes_pseudo_time,ripple_at_spikes(:,3),'.k','MarkerSize',10)
        xlim([1 max(ripple_pseudo_time)])
        ax2 = subplot(3,2,2);
        plot(spikes_pseudo_time,spkData_in_event(:,3),'.k','MarkerSize',10)
        ylim([0 2*pi]);
        ax3 = subplot(3,2,3);
        imagesc(t_replay.posterior_full_left{sub_event_num}'); set(gca,'YDir','normal')
        ax5 = subplot(3,2,5);
        imagesc(t_replay.posterior_full_right{sub_event_num}'); set(gca,'YDir','normal')
        ax4 = subplot(3,2,4);
        plot(clusters(i).runs(1).directions(1).rateMap_smoothed,[1:length(clusters(i).runs(1).directions(1).rateMap_smoothed)])
        ylim([0 length(clusters(i).runs(1).directions(1).rateMap_smoothed)])
        ax6 = subplot(3,2,6);
        plot(clusters(i).runs(1).directions(2).rateMap_smoothed,[1:length(clusters(i).runs(1).directions(2).rateMap_smoothed)])


        linkaxes([ax3,ax4,ax5,ax6],'y')
        linkaxes([ax1,ax3,ax5],'x')



         linkaxes([ax1,ax2,ax3],'x')
        drawnow
        pause
       

        a_all_fractional_time_in_event = [a_all_fractional_time_in_event; fractional_time_in_event];
        a_all_ripple_phase_in_event = [a_all_ripple_phase_in_event; spkData_in_event(:,3)];
        a_all_unwrapped_ripple_phase_in_event = [a_all_unwrapped_ripple_phase_in_event; unwrap(spkData_in_event(:,3))];
        a_all_event_number = [a_all_event_number; sub_event_num*ones(size(spkData_in_event,1),1)];
        a_all_event_direction = [a_all_event_direction; sub_replays.direction(sub_event_num)*ones(size(spkData_in_event,1),1)];
        a_all_event_laser_state = [a_all_event_laser_state; sub_replays.laser_state(sub_event_num)*ones(size(spkData_in_event,1),1)];
        a_all_event_best_map = [a_all_event_best_map; sub_replays.best_map(sub_event_num)*ones(size(spkData_in_event,1),1)];
        a_all_event_is_replay = [a_all_event_is_replay; sub_replays.replay(sub_event_num)*ones(size(spkData_in_event,1),1)];
        a_all_event_cluster_id = [a_all_event_cluster_id; i*ones(size(spkData_in_event,1),1)];
    end
end
cluster_phase_precession = table();
cluster_phase_precession.cluster_id = a_all_event_cluster_id;
cluster_phase_precession.time_in_event = a_all_fractional_time_in_event;
cluster_phase_precession.phase = a_all_ripple_phase_in_event;
cluster_phase_precession.unwrapped_phase = a_all_unwrapped_ripple_phase_in_event;
cluster_phase_precession.event_number = a_all_event_number;
cluster_phase_precession.event_direction = a_all_event_direction;
cluster_phase_precession.laser_state = a_all_event_laser_state;
cluster_phase_precession.is_replay = a_all_event_is_replay;
cluster_phase_precession.best_map = a_all_event_best_map;
end
%% In progress: look at phase precession to ripples for each field

if analyze_distance_versus_phase==1
b_all_distance_through_field = [];
b_all_ripple_phase_in_event = [];
b_all_unwrapped_ripple_phase_in_event = [];
b_all_event_number = [];
b_all_subfield_num = [];
b_all_replay_direction  = [];
b_all_event_cluster_id = [];
b_all_event_laser_state = [];


for i = 1:length(clusters)

    spkTimes = clusters(i).spkTime;
    spkLaserStatus = clusters(i).spkLaserState;
    spkRipplePhase = clusters(i).spkRipplePhase;
    spkPos = clusters(i).spkPos(:,1);
    spkData = [spkTimes spkLaserStatus spkRipplePhase spkPos];

    % B) Now plot as distance through field (0 to 1) versus spike
    % ripple phase
    % For each cell, make a table containing its place fields and the
    % direction of movement for that field.
    fields = table();
    sub_field_number = 0;
    for direction=1:2
        num_fields = max(clusters(i).runs(1).directions(direction).Place_Field);
        if num_fields>0
            for sub_field = 1:num_fields
                sub_field_number = sub_field_number + 1;
                edge1=find(clusters(i).runs(1).directions(direction).Place_Field == sub_field,1,'first');
                edge2=find(clusters(i).runs(1).directions(direction).Place_Field == sub_field,1,'last');
                fields.edges(sub_field_number,:) = [edge1 edge2];
                fields.direction(sub_field_number) = direction;
            end
        end
    end

    for subfield = 1:height(fields)
        map = fields.direction(subfield);
        % For this field, only consider replays that utilized this map.
        sub_replays = t_replay(t_replay.best_map == map,:);

     
        for sub_event_num = 1:height(sub_replays)
            % Loop through all events of this type:
            if map == 1
                posterior = sub_replays.posterior_left{sub_event_num};
            else
                posterior = sub_replays.posterior_right{sub_event_num};
            end
            % Decoded position: using com
            x = sub_replays.com{sub_event_num};
            % Decoded position: using peak
%             [~,x] = max(posterior,[],2);
            decoding_timeBin_centers = mean(sub_replays.timeBins{sub_event_num,:},2);
            sub_replay_timePoints = sub_replays.timePoints(sub_event_num,:);
            %     figure()
            %     plot(decoding_timeBin_centers,x);

            % Find the time that the replay entered the cell's field
            if  ((sub_replays.direction(sub_event_num) == 1 && sub_replays.best_map(sub_event_num) == 1) ...
                    || (sub_replays.direction(sub_event_num) == 2 && sub_replays.best_map(sub_event_num) == 2))
                % replay is going from right to left (or top to bottom)
                field_entry = fields.edges(subfield,2);
                field_exit = fields.edges(subfield,1);
                entry_ind = find(x<=field_entry,1,'first');
                exit_ind = find(x<=field_exit,1,'first');
            else
                % replay is going from left to right (forward using right map, or
                % reverse using left map)
                field_entry = fields.edges(subfield,1);
                field_exit = fields.edges(subfield,2);
                entry_ind = find(x>=field_entry,1,'first');
                exit_ind = find(x>=field_exit,1,'first');
            end

            spkData_in_replay =  compute_dataTemporalConcatenation(spkData,sub_replay_timePoints);
            decoded_position_data = [decoding_timeBin_centers,x];
            decoded_position_at_time_of_spks= compute_dataInterpolation(decoded_position_data,spkData_in_replay(:,1),[]);
            spkData_in_subfield = [spkData_in_replay decoded_position_at_time_of_spks];
            % Throw away spikes outside this subfield
            spkData_in_subfield(spkData_in_subfield(:,6)<fields.edges(subfield,1) | ...
                spkData_in_subfield(:,6)>fields.edges(subfield,2),:) = [];

            distance_through_field = abs( spkData_in_subfield(:,6)-field_entry)./abs(field_entry-field_exit);

            b_all_distance_through_field = [b_all_distance_through_field; distance_through_field];
            b_all_ripple_phase_in_event = [b_all_ripple_phase_in_event; spkData_in_subfield(:,3)];
            b_all_unwrapped_ripple_phase_in_event = [b_all_unwrapped_ripple_phase_in_event; unwrap(spkData_in_subfield(:,3))];
            b_all_event_number = [b_all_event_number; sub_event_num*(ones(size(spkData_in_subfield,1),1))];
            b_all_subfield_num = [b_all_subfield_num; subfield*(ones(size(spkData_in_subfield,1),1))];
            b_all_replay_direction = [b_all_replay_direction; sub_replays.direction(sub_event_num)*ones(size(spkData_in_subfield,1),1)];
            b_all_event_cluster_id = [b_all_event_cluster_id; i*ones(size(spkData_in_subfield,1),1)];
            b_all_event_laser_state = [b_all_event_laser_state; sub_replays.laser_state(sub_event_num)*ones(size(spkData_in_subfield,1),1)];
        end
    end
end

cluster_field_phase_precession = table();
cluster_field_phase_precession.cluster_id = b_all_event_cluster_id;
cluster_field_phase_precession.distance_through_field =b_all_distance_through_field;
cluster_field_phase_precession.phase = b_all_ripple_phase_in_event;
cluster_field_phase_precession.unwrapped_phase = b_all_unwrapped_ripple_phase_in_event;
cluster_field_phase_precession.event_number = b_all_event_number;
cluster_field_phase_precession.event_direction = b_all_replay_direction;
cluster_field_phase_precession.laser_state = b_all_event_laser_state;
cluster_field_phase_precession.subfield = b_all_subfield_num;

cluster_field_phase_precession.unique_field_id(1) = 1;
for i = 2:height(cluster_field_phase_precession)
    if cluster_field_phase_precession.cluster_id(i) == cluster_field_phase_precession.cluster_id(i-1) && ...
            cluster_field_phase_precession.subfield(i) == cluster_field_phase_precession.subfield(i-1)
        cluster_field_phase_precession.unique_field_id(i) = cluster_field_phase_precession.unique_field_id(i-1);
    else
        cluster_field_phase_precession.unique_field_id(i) = cluster_field_phase_precession.unique_field_id(i-1)+1;
    end
end
cluster_field_phase_precession(isnan(cluster_field_phase_precession.distance_through_field),:) = [];
end
%%
numShuff = 0;
x = linspace(0,1,100);
field_precession_summary = table();

all_replay2_table = table(); all_replay_row = 0;
all_replays_laser_off_table = table(); all_replays_laser_off_row = 0;
all_replays_laser_on_table = table(); all_replays_laser_on_row = 0;
f_table = table(); f_row = 0;
f1_table = table(); f1_row = 0;
f0_table = table(); f0_row = 0;
r_table = table(); r_row = 0;
r1_table = table(); r1_row = 0;
r0_table = table(); r0_row = 0;

for subfield = 1:max(cluster_field_phase_precession.unique_field_id)
    plot_fig = 0;
    % All replays::
    sub_table = cluster_field_phase_precession(cluster_field_phase_precession.unique_field_id == subfield,:);
    if isempty(sub_table)
        continue
    end


    cluster = unique(sub_table.cluster_id);
    cluster_subfield = unique(sub_table.subfield);
        % B) Now plot as distance through field (0 to 1) versus spike
    % ripple phase
    % For each cell, make a table containing its place fields and the
    % direction of movement for that field.
    fields = table();
    sub_field_number = 0;
    for direction=1:2
        num_fields = max(clusters(cluster).runs(1).directions(direction).Place_Field);
        if num_fields>0
            for sub_field = 1:num_fields
                sub_field_number = sub_field_number + 1;
                edge1=find(clusters(cluster).runs(1).directions(direction).Place_Field == sub_field,1,'first');
                edge2=find(clusters(cluster).runs(1).directions(direction).Place_Field == sub_field,1,'last');
                fields.edges(sub_field_number,:) = [edge1 edge2];
                fields.direction(sub_field_number) = direction;
            end
        end
    end

    
    num_spks = height(sub_table);
    [slope,intercept,rho, pval, pval_shuffle,slope_CI] = circular_linear_analysis(sub_table.phase, sub_table.distance_through_field, numShuff);
    regLine = mod(2*pi*slope*x + intercept,2*pi);

    if num_spks > 4 && max(sub_table.distance_through_field) - min(sub_table.distance_through_field) >= 0.5
        plot_fig=1;
        all_replay_row = all_replay_row  + 1;
        all_replay2_table.slope(all_replay_row ) = slope;
        all_replay2_table.intercept(all_replay_row ) = intercept;
        all_replay2_table.rho(all_replay_row )= rho;
        all_replay2_table.pval(all_replay_row ) = pval;
        all_replay2_table.field_id(all_replay_row ) = subfield;
        all_replay2_table.cluster_id(all_replay_row ) = unique(sub_table.cluster_id);
        all_replay2_table.excitatory(all_replay_row ) = clusters(unique(sub_table.cluster_id)).Excitatory;
    else
        plot_fig =0;
    end
    if plot_fig==1
        figure('Position',[234 151 1585 755])
        subplot(2,5,2)
        plot([sub_table.distance_through_field; sub_table.distance_through_field],[sub_table.phase; sub_table.phase + 2*pi],'ok','MarkerFaceColor','k')
        xlim([0 1])
        ylim([0 4*pi])
        hold on
        plot([x x],[regLine regLine + 2*pi],'.r')
        hold off
        title(['All replays p=' num2str(pval,2)])
    end


    if plot_fig==1
        subplot(2,5,1);
        plot(clusters(cluster).runs(1).directions(1).rateMap_smoothed,'-b'); hold on
        plot(clusters(cluster).runs(1).directions(2).rateMap_smoothed,'-r'); hold on
        
        cluster_subfield_edges = fields.edges(cluster_subfield,:);
        hold on
        if fields.direction(cluster_subfield) == 1
        xline(cluster_subfield_edges(1),'b'); hold on
        xline(cluster_subfield_edges(2),'b');
        else
        xline(cluster_subfield_edges(1),'r'); hold on
        xline(cluster_subfield_edges(2),'r'); 
        end
    end
    
    if plot_fig==1
        subplot(2,5,6)
        rose(clusters(cluster).spkRipplePhase)
    end







        % All laser off replays:
    sub_table = cluster_field_phase_precession(cluster_field_phase_precession.unique_field_id == subfield & ...
        cluster_field_phase_precession.laser_state == 0,:);
    num_spks = height(sub_table);
    [slope,intercept,rho, pval, pval_shuffle,slope_CI] = circular_linear_analysis(sub_table.phase, sub_table.distance_through_field, numShuff);
    regLine = mod(2*pi*slope*x + intercept,2*pi);

    if plot_fig==1
        subplot(2,5,7)
        plot([sub_table.distance_through_field; sub_table.distance_through_field],[sub_table.phase; sub_table.phase + 2*pi],'ok','MarkerFaceColor','k')
        xlim([0 1])
        ylim([0 4*pi])
        hold on
        plot([x x],[regLine regLine + 2*pi],'.r')
        title(['All laser off p=' num2str(pval,2)])
    end
    if num_spks > 4 && max(sub_table.distance_through_field) - min(sub_table.distance_through_field) >= 0.5
        all_replays_laser_off_row = all_replays_laser_off_row + 1;
        all_replays_laser_off_table.slope(all_replays_laser_off_row) = slope;
        all_replays_laser_off_table.intercept(all_replays_laser_off_row) = intercept;
        all_replays_laser_off_table.rho(all_replays_laser_off_row)= rho;
        all_replays_laser_off_table.pval(all_replays_laser_off_row) = pval;
        all_replays_laser_off_table.field_id(all_replays_laser_off_row) = subfield;
        all_replays_laser_off_table.cluster_id(all_replays_laser_off_row) = unique(sub_table.cluster_id);
        all_replays_laser_off_table.excitatory(all_replays_laser_off_row) = clusters(unique(sub_table.cluster_id)).Excitatory;
    end

        % All laser on replays:
    sub_table = cluster_field_phase_precession(cluster_field_phase_precession.unique_field_id == subfield & ...
        cluster_field_phase_precession.laser_state == 1,:);
    num_spks = height(sub_table);
    [slope,intercept,rho, pval, pval_shuffle,slope_CI] = circular_linear_analysis(sub_table.phase, sub_table.distance_through_field, numShuff);
    regLine = mod(2*pi*slope*x + intercept,2*pi);

%     if plot_fig==1
%         subplot(2,4,5)
%         plot([sub_table.distance_through_field; sub_table.distance_through_field],[sub_table.phase; sub_table.phase + 2*pi],'ok','MarkerFaceColor','k')
%         xlim([0 1])
%         ylim([0 4*pi])
%         hold on
%         plot([x x],[regLine regLine + 2*pi],'.r')
%         title(['All laser on p=' num2str(pval,2)])
%     end
    if num_spks > 4 && max(sub_table.distance_through_field) - min(sub_table.distance_through_field) >= 0.5
        all_replays_laser_on_row = all_replays_laser_on_row + 1;
        all_replays_laser_on_table.slope(all_replays_laser_on_row) = slope;
        all_replays_laser_on_table.intercept(all_replays_laser_on_row) = intercept;
        all_replays_laser_on_table.rho(all_replays_laser_on_row)= rho;
        all_replays_laser_on_table.pval(all_replays_laser_on_row) = pval;
        all_replays_laser_on_table.field_id(all_replays_laser_on_row) = subfield;
        all_replays_laser_on_table.cluster_id(all_replays_laser_on_row) = unique(sub_table.cluster_id);
        all_replays_laser_on_table.excitatory(all_replays_laser_on_row) = clusters(unique(sub_table.cluster_id)).Excitatory;
    end

    % All forward replays:
    sub_table = cluster_field_phase_precession(cluster_field_phase_precession.unique_field_id == subfield & ...
        cluster_field_phase_precession.event_direction==1,:);
    num_spks = height(sub_table);
    [slope,intercept,rho, pval, pval_shuffle,slope_CI] = circular_linear_analysis(sub_table.phase, sub_table.distance_through_field, numShuff);
    regLine = mod(2*pi*slope*x + intercept,2*pi);
    if plot_fig==1
        subplot(2,5,3)
        plot([sub_table.distance_through_field; sub_table.distance_through_field],[sub_table.phase; sub_table.phase + 2*pi],'ok','MarkerFaceColor','k')
        xlim([0 1])
        ylim([0 4*pi])
        hold on
        plot([x x],[regLine regLine + 2*pi],'.r')
        hold off
        title(['Forward p=' num2str(pval,2)])
    end

    if num_spks > 4 && max(sub_table.distance_through_field) - min(sub_table.distance_through_field) >= 0.5
        f_row = f_row + 1;
        f_table.slope(f_row) = slope;
        f_table.intercept(f_row) = intercept;
        f_table.rho(f_row)= rho;
        f_table.pval(f_row) = pval;
        f_table.field_id(f_row) = subfield;
        f_table.cluster_id(f_row) = unique(sub_table.cluster_id);
        f_table.excitatory(f_row) = clusters(unique(sub_table.cluster_id)).Excitatory;
    end

    % All forward off replays:
    sub_table = cluster_field_phase_precession(cluster_field_phase_precession.unique_field_id == subfield & ...
        cluster_field_phase_precession.event_direction==1 & cluster_field_phase_precession.laser_state==0,:);
    num_spks = height(sub_table);
    [slope,intercept,rho, pval, pval_shuffle,slope_CI] = circular_linear_analysis(sub_table.phase, sub_table.distance_through_field, numShuff);
    regLine = mod(2*pi*slope*x + intercept,2*pi);
    if plot_fig==1
        subplot(2,5,4)
        plot([sub_table.distance_through_field; sub_table.distance_through_field],[sub_table.phase; sub_table.phase + 2*pi],'ok','MarkerFaceColor','k')
        xlim([0 1])
        ylim([0 4*pi])
        hold on
        plot([x x],[regLine regLine + 2*pi],'.r')
        hold off
        title(['Forward off p=' num2str(pval,2)])
    end
    if num_spks > 4 && max(sub_table.distance_through_field) - min(sub_table.distance_through_field) >= 0.5
        f0_row = f0_row + 1;
        f0_table.slope(f0_row) = slope;
        f0_table.intercept(f0_row) = intercept;
        f0_table.rho(f0_row)= rho;
        f0_table.pval(f0_row) = pval;
        f0_table.field_id(f0_row) = subfield;
        f0_table.cluster_id(f0_row) = unique(sub_table.cluster_id);
        f0_table.excitatory(f0_row) = clusters(unique(sub_table.cluster_id)).Excitatory;
    end

    % Forward On Replays:
    sub_table = cluster_field_phase_precession(cluster_field_phase_precession.unique_field_id == subfield & ...
        cluster_field_phase_precession.event_direction==1 & cluster_field_phase_precession.laser_state==1,:);
    num_spks = height(sub_table);
    [slope,intercept,rho, pval, pval_shuffle,slope_CI] = circular_linear_analysis(sub_table.phase, sub_table.distance_through_field, numShuff);
    regLine = mod(2*pi*slope*x + intercept,2*pi);

    if plot_fig==1
        subplot(2,5,5)
        plot([sub_table.distance_through_field; sub_table.distance_through_field],[sub_table.phase; sub_table.phase + 2*pi],'ok','MarkerFaceColor','k')
        xlim([0 1])
        ylim([0 4*pi])
        hold on
        plot([x x],[regLine regLine + 2*pi],'.r')
        hold off
        title(['Forward on p=' num2str(pval,2)])
    end
    if num_spks > 4 && max(sub_table.distance_through_field) - min(sub_table.distance_through_field) >= 0.5
        f1_row = f1_row+1;
        f1_table.slope(f1_row) = slope;
        f1_table.intercept(f1_row) = intercept;
        f1_table.rho(f1_row)= rho;
        f1_table.pval(f1_row) = pval;
        f1_table.field_id(f1_row) = subfield;
        f1_table.cluster_id(f1_row) = unique(sub_table.cluster_id);
        f1_table.excitatory(f1_row) = clusters(unique(sub_table.cluster_id)).Excitatory;
    end


    % All reverse replays:
    sub_table = cluster_field_phase_precession(cluster_field_phase_precession.unique_field_id == subfield & ...
        cluster_field_phase_precession.event_direction==2,:);
    num_spks = height(sub_table);
    [slope,intercept,rho, pval, pval_shuffle,slope_CI] = circular_linear_analysis(sub_table.phase, sub_table.distance_through_field, numShuff);
    regLine = mod(2*pi*slope*x + intercept,2*pi);

    if plot_fig==1
        subplot(2,5,8)
        plot([sub_table.distance_through_field; sub_table.distance_through_field],[sub_table.phase; sub_table.phase + 2*pi],'ok','MarkerFaceColor','k')
        xlim([0 1])
        ylim([0 4*pi])
        hold on
        plot([x x],[regLine regLine + 2*pi],'.r')
        title(['Reverse p=' num2str(pval,2)])
    end
    if num_spks > 4 && max(sub_table.distance_through_field) - min(sub_table.distance_through_field) >= 0.5
        r_row = r_row + 1;
        r_table.slope(r_row) = slope;
        r_table.intercept(r_row) = intercept;
        r_table.rho(r_row)= rho;
        r_table.pval(r_row) = pval;
        r_table.field_id(r_row) = subfield;
        r_table.cluster_id(r_row) = unique(sub_table.cluster_id);
        r_table.excitatory(r_row) = clusters(unique(sub_table.cluster_id)).Excitatory;
    end

    % Reverse off replays
    sub_table = cluster_field_phase_precession(cluster_field_phase_precession.unique_field_id == subfield & ...
        cluster_field_phase_precession.event_direction==2  & cluster_field_phase_precession.laser_state==0,:);
    num_spks = height(sub_table);
    [slope,intercept,rho, pval, pval_shuffle,slope_CI] = circular_linear_analysis(sub_table.phase, sub_table.distance_through_field, numShuff);
    regLine = mod(2*pi*slope*x + intercept,2*pi);
    if plot_fig==1
        subplot(2,5,9)
        plot([sub_table.distance_through_field; sub_table.distance_through_field],[sub_table.phase; sub_table.phase + 2*pi],'ok','MarkerFaceColor','k')
        xlim([0 1])
        ylim([0 4*pi])
        hold on
        plot([x  x],[regLine regLine + 2*pi],'.r')
        hold off
        title(['Reverse off p=' num2str(pval,2)])
    end
    if num_spks > 4 && max(sub_table.distance_through_field) - min(sub_table.distance_through_field) >= 0.5
        r0_row = r0_row + 1;
        r0_table.slope(r0_row) = slope;
        r0_table.intercept(r0_row) = intercept;
        r0_table.rho(r0_row)= rho;
        r0_table.pval(r0_row) = pval;
        r0_table.field_id(r0_row) = subfield;
        r0_table.cluster_id(r0_row) = unique(sub_table.cluster_id);
        r0_table.excitatory(r0_row) = clusters(unique(sub_table.cluster_id)).Excitatory;
    end

    % Reverse on repalys:
    sub_table = cluster_field_phase_precession(cluster_field_phase_precession.unique_field_id == subfield & ...
        cluster_field_phase_precession.event_direction==2  & cluster_field_phase_precession.laser_state==1,:);
    num_spks = height(sub_table);
    [slope,intercept,rho, pval, pval_shuffle,slope_CI] = circular_linear_analysis(sub_table.phase, sub_table.distance_through_field, numShuff);
    regLine = mod(2*pi*slope*x + intercept,2*pi);
    if plot_fig==1
        subplot(2,5,10)
        plot([sub_table.distance_through_field; sub_table.distance_through_field],[sub_table.phase; sub_table.phase + 2*pi],'ok','MarkerFaceColor','k')
        xlim([0 1])
        ylim([0 4*pi])
        hold on
        plot([x x],[regLine regLine + 2*pi],'.r')
        hold off
        title(['Reverse on p=' num2str(pval,2)])
        sgtitle(['Cluster ' num2str(unique(sub_table.cluster_id)) ' field ' num2str(unique(sub_table.subfield))])
    
    export_fig(gcf,fullfile(fig_folder,['Cluster ' num2str(unique(sub_table.cluster_id)) ' field ' num2str(unique(sub_table.subfield))]),'-jpeg')


    end
    if num_spks > 4 && max(sub_table.distance_through_field) - min(sub_table.distance_through_field) >= 0.5
        r1_row = r1_row + 1;
        r1_table.slope(r1_row) = slope;
        r1_table.intercept(r1_row) = intercept;
        r1_table.rho(r1_row)= rho;
        r1_table.pval(r1_row) = pval;
        r1_table.field_id(r1_row) = subfield;
        r1_table.cluster_id(r1_row) = unique(sub_table.cluster_id);
        r1_table.excitatory(r1_row) = clusters(unique(sub_table.cluster_id)).Excitatory;
    end
    drawnow
end

% figure()
% hist_bin_edges = linspace(-2,2,10);
% subplot(2,4,1)
% hist(all_replay_table.slope(all_replay_table.excitatory==1 & all_replay_table.pval<0.05),hist_bin_edges)
% subplot(2,4,2)
% hist(f_table.slope(f_table.excitatory==1 & f_table.pval<0.05),hist_bin_edges)
% subplot(2,4,3)
% hist(f0_table.slope(f0_table.excitatory==1 & f0_table.pval<0.05),hist_bin_edges)
% subplot(2,4,4)
% hist(f1_table.slope(f1_table.excitatory==1 & f1_table.pval<0.05),hist_bin_edges)
% subplot(2,4,5)
% hist(all_replays_laser_off_table.slope(all_replays_laser_off_table.excitatory==1 & all_replays_laser_off_table.pval<0.05),hist_bin_edges)
% subplot(2,4,6)
% hist(r_table.slope(r_table.excitatory==1 & r_table.pval<0.05),hist_bin_edges)
% subplot(2,4,7)
% hist(r0_table.slope(r0_table.excitatory==1 & r0_table.pval<0.05),hist_bin_edges)
% subplot(2,4,8)
% hist(r1_table.slope(r1_table.excitatory==1 & r1_table.pval<0.05),hist_bin_edges)

ripple_precession = struct();
ripple_precession.all_replays = all_replay2_table;
ripple_precession.all_replays_laser_off = all_replays_laser_off_table;
ripple_precession.all_replays_laser_on = all_replays_laser_on_table;
ripple_precession.forward_replays = f_table;
ripple_precession.forward_off_replays = f0_table;
ripple_precession.forward_on_replays = f1_table;
ripple_precession.reverse_replays = r_table;
ripple_precession.reverse_off_replays = r0_table;
ripple_precession.reverse_on_replays = r1_table;

save('ripple_phase_precession','cluster_field_phase_precession','ripple_precession')
%%
% plot(a_all_fractional_time_in_event(a_all_event_laser_state==0 & a_all_event_is_replay==1),a_all_unwrapped_ripple_phase_in_event(a_all_event_laser_state==0 & a_all_event_is_replay==1),'ok')
% 
% plot(a_all_fractional_time_in_event(a_all_event_laser_state==0 & a_all_event_is_replay==1 & a_all_event_direction == 1), ...
%     a_all_unwrapped_ripple_phase_in_event(a_all_event_laser_state==0 & a_all_event_is_replay==1 & a_all_event_direction == 1),'ok')
% plot(a_all_fractional_time_in_event(a_all_event_laser_state==1 & a_all_event_is_replay==1 & a_all_event_direction == 1), ...
%     a_all_unwrapped_ripple_phase_in_event(a_all_event_laser_state==1 & a_all_event_is_replay==1 & a_all_event_direction == 1),'ok')
% plot(a_all_fractional_time_in_event(a_all_event_laser_state==0 & a_all_event_is_replay==1 & a_all_event_direction == 2), ...
%     a_all_unwrapped_ripple_phase_in_event(a_all_event_laser_state==0 & a_all_event_is_replay==1 & a_all_event_direction == 2),'ok')
% plot(a_all_fractional_time_in_event(a_all_event_laser_state==1 & a_all_event_is_replay==1 & a_all_event_direction == 2), ...
%     a_all_unwrapped_ripple_phase_in_event(a_all_event_laser_state==1 & a_all_event_is_replay==1 & a_all_event_direction == 2),'ok')
% 
% 
% plot(b_all_distance_through_field,b_all_unwrapped_ripple_phase_in_event,'ok')
% 
% plot(b_all_distance_through_field(b_all_event_laser_state==0 & b_all_replay_direction == 1), ...
%     b_all_unwrapped_ripple_phase_in_event(b_all_event_laser_state==0 & b_all_replay_direction == 1),'ok')
% plot(b_all_distance_through_field(b_all_event_laser_state==1 & b_all_replay_direction == 1), ...
%     b_all_unwrapped_ripple_phase_in_event(b_all_event_laser_state==1  & b_all_replay_direction == 1),'ok')
% plot(b_all_distance_through_field(b_all_event_laser_state==0 & b_all_replay_direction == 2), ...
%     b_all_unwrapped_ripple_phase_in_event(b_all_event_laser_state==0  & b_all_replay_direction == 2),'ok')
% plot(b_all_distance_through_field(b_all_event_laser_state==1 &  b_all_replay_direction == 2), ...
%     b_all_unwrapped_ripple_phase_in_event(b_all_event_laser_state==1 & b_all_replay_direction == 2),'ok')
% 
% plot(b_all_distance_through_field(b_all_replay_direction == 1), ...
%     b_all_unwrapped_ripple_phase_in_event(b_all_replay_direction == 1),'ok')
% plot(b_all_distance_through_field(b_all_replay_direction == 2), ...
%     b_all_unwrapped_ripple_phase_in_event(b_all_replay_direction == 2),'ok')
