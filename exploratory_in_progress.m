% Loops through all sessions, loading clusters. Looks at the participation
% of different cell types (pyramidal, inhibitory, unimodal, or bimodal) in
% either ripples, sdes or replays of different types. Also examines the
% impact of laser.

restrict_to_well_isolated_cells = 0;
num_fields_thr = -inf;
min_peak_fr_thr = 0;

windows = 0;
rats = [1 2 3 6];

if windows == 1
    fig_path = 'D:/Dropbox/Foster Lab/Data/Replay summary/cell_type_participation';
else
    fig_path = '/home/caitlin/Insync/caitlinmallory@berkeley.edu/Google Drive/Processed_Data/cell_type_participation';
end

must_include_all_flags = 1; % 1: requires that all flags in Flags_to_include are present. 0: requires that at least one flag in Flags_to_include are present.
novel_range = [10.1 10.2];

flags_to_include = 11;
flags_to_exclude = [0];
flags.flags_to_include = flags_to_include;
flags.flags_to_exclude = flags_to_exclude;

Rat_Names = {'Clover','Bo','MEC1','Bolt','Dash','CM1'};
table_properties = {'spkTime','Tetrode','spkPos','spkSpeed','spkHD','Modality','Excitatory','well_isolated','mvl_laser_off','mvl_laser_on','Rayleigh_P_laser_off',...
    'cell_participation_metrics','ripple_phase_locking'};

all_clusters = table;

session_id = 0;

for rat = rats
    if rat == 1 %Clover

        if windows == 1
            directory = 'G:/My Drive/Processed_Data/Clover/linear_track';
        else
            directory = '/home/caitlin/Insync/caitlinmallory@berkeley.edu/Google Drive/Processed_Data/Clover/linear_track';
        end
        % use the sessions with the best decoding:
        dayFiles = {'20210427/hand_clustered','20210428/hand_clustered','20210429/kilosort_clustered','20210430/kilosort_clustered','20210501/hand_clustered','20210507/kilosort_clustered','20210508/hand_clustered','20210509/kilosort_clustered','20210512/kilosort_clustered','20210517/hand_clustered','20210519/kilosort_clustered'};

    end
    if rat == 2 %Bo
        if windows == 1
            directory = 'G:/My Drive/Processed_Data/Bo/linear_track';
        else
            directory = '/home/caitlin/Insync/caitlinmallory@berkeley.edu/Google Drive/Processed_Data/Bo/linear_track';
        end
        %best decoding
        dayFiles = {'20210511/kilosort_clustered','20210512/kilosort_clustered','20210513/kilosort_clustered','20210517/kilosort_clustered'};

    end
    if rat == 3 % MEC1
        if windows == 1
            directory = 'G:/My Drive/Processed_Data/MEC1';
        else
            directory = '/home/caitlin/Insync/caitlinmallory@berkeley.edu/Google Drive/Processed_Data/MEC1';
        end
        dayFiles = {'20200926/kilosort_clustered','20200927/hand_clustered2'};
    end
    if rat == 4
        if windows == 1
            % directory = 'L:/';
            directory = 'G:/My Drive/Processed_Data/Bolt/linear_track';
        else
            directory = '/home/caitlin/Insync/caitlinmallory@berkeley.edu/Google Drive/Processed_Data/Bolt/linear_track';
        end
        dayFiles = {'20220714/msort_clustered','20220715/msort_clustered','20220718_1/msort_clustered','20220718_2/msort_clustered','20220719_1_opto/msort_clustered_full','20220720/msort_clustered','20220722/msort_clustered_full','20220729/msort_clustered'};
    elseif rat == 5 % Dash
        if windows == 1
            directory = 'G:/My Drive/Processed_Data/Dash';
        else
            directory = '/home/caitlin/Insync/caitlinmallory@berkeley.edu/Google Drive/Processed_Data/Dash';
        end
        dayFiles = {'20221025_113913'; '20221025_170200'};
    elseif rat == 6
        if windows == 1
            directory = 'G:/My Drive/Processed_Data/CM1';
        else
            directory = '/home/caitlin/Insync/caitlinmallory@berkeley.edu/Google Drive/Processed_Data/CM1';
        end
        dayFiles = ...
            {'20230217','20230219','20230219_01','20230220','20230221','20230221_2_bad','20230221_3','20230222_1','20230222_2','20230224_2','20230228','20230301_01','20230302_01','20230302_02','20230303_01','20230303_02','20230305_01','20230305_02'}; 
    end

    for day = 1:length(dayFiles)
        display([Rat_Names{rat},' Day ', dayFiles{day}])
        cd(fullfile(directory,dayFiles{day}))

        load Experiment_Information

        sessions_that_meet_criterion_day = check_for_relevant_sessions(flags,1,Experiment_Information,must_include_all_flags,1,novel_range);

        if ~isempty(sessions_that_meet_criterion_day)

            for session_count = 1:length(sessions_that_meet_criterion_day)
                session_id = session_id + 1;
                load clusters
                clusters_sub = struct2table(clusters);

                if any("tetrode" == string(clusters_sub.Properties.VariableNames))
                    clusters_sub = renamevars(clusters_sub,'tetrode','Tetrode');
                end
                clusters_sub = clusters_sub(:,table_properties);
                % Add rat and session info to clusters
                clusters_sub.Rat = repmat(rat,[height(clusters_sub),1]);
                clusters_sub.session_path = repmat(dayFiles(day),[height(clusters_sub),1]);
                clusters_sub.unique_session_id = repmat(session_id,[height(clusters_sub),1]);

                for n = 1:length(clusters)
                    clusters_sub.peak_fr(n) = clusters(n).runs(1).maxSpatialFiringRate;
                    clusters_sub.num_subFields(n) = clusters(n).runs(1).num_subFields;
                end
                all_clusters = [all_clusters; clusters_sub];

            end

        end
    end
    cd ..
end

%%
if restrict_to_well_isolated_cells == 0
    all_clusters.meets_well_isolated_requirement = ones(height(all_clusters),1);
else
    all_clusters.meets_well_isolated_requirement = all_clusters.well_isolated;
end

% Require unimodal and bimodal cells to have at least one place field for
% consideration
excitatory = all_clusters.Excitatory == 1 & all_clusters.meets_well_isolated_requirement == 1  & all_clusters.peak_fr > min_peak_fr_thr;
inhibitory = all_clusters.Excitatory == 0 & all_clusters.meets_well_isolated_requirement == 1  & all_clusters.peak_fr > min_peak_fr_thr;
unimodal = all_clusters.Modality == 1 & excitatory & all_clusters.num_subFields > num_fields_thr;
bimodal = all_clusters.Modality == 2 & excitatory > num_fields_thr;

inds = find(all_clusters.Excitatory==1);


inds = find(all_clusters.cell_participation_metrics.laser_off_reverse_replays.pcnt_participation > 0.2 & ...
    all_clusters.cell_participation_metrics.laser_on_reverse_replays.pcnt_participation > 0.2 & ... 
    all_clusters.Excitatory == 1);


inds = find(all_clusters.cell_participation_metrics.all_reverse_replays.firing_rate_in_event_type./...
    all_clusters.cell_participation_metrics.all_forward_replays.firing_rate_in_event_type > 2 ...
    & all_clusters.Excitatory == 1);






a = all_clusters.ripple_phase_locking.laser_off_reverse.mvl(inds);
b = all_clusters.ripple_phase_locking.laser_on_reverse.mvl(inds);

a = all_clusters.ripple_phase_locking.laser_off_reverse.peak_phase(inds);
b = all_clusters.ripple_phase_locking.laser_on_reverse.peak_phase(inds);



a = all_clusters.cell_participation_metrics.laser_off_reverse_replays.spikes_per_event(inds)
b = all_clusters.cell_participation_metrics.laser_on_reverse_replays.spikes_per_event(inds)

a = all_clusters.cell_participation_metrics.laser_off_reverse_replays.firing_rate_in_event_type(inds)
b = all_clusters.cell_participation_metrics.laser_on_reverse_replays.firing_rate_in_event_type(inds)

a = all_clusters.cell_participation_metrics.laser_off_reverse_replays.pcnt_participation(inds)
b = all_clusters.cell_participation_metrics.laser_on_reverse_replays.pcnt_participation(inds)

signrank(a,b)

figure();
plot(a,b,'.k'); hold on
axis equal
axis square
ax = gca
plot([ax.XLim(1) ax.XLim(2)],[ax.YLim(1) ax.YLim(2)],'r')

figure()
centers = linspace(-1,1,20)
hist(a-b,centers)
xline(0)

figure();
ecdf(a); hold on;
ecdf(b);


rose(a); hold on; rose(b);


























































