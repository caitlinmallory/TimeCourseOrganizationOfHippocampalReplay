% Set default axis colors for plots to black
set(groot, {'DefaultAxesXColor', 'DefaultAxesYColor', 'DefaultAxesZColor'}, {'k', 'k', 'k'})

% Define the rat IDs for the experiment
rats = [1 4 6 12 13 14];  % List of rats in the experiment
flags_to_include = [11];   % Flags to include in analysis (e.g., specific conditions)
% Define the flags to exclude, excluding data with MEC (Medial Entorhinal Cortex) activity
flags_to_exclude = [100 20 18 14];  % Flags indicating invalid or excluded data

% Structure to store the flags for inclusion and exclusion
flags.flags_to_include = flags_to_include;
flags.flags_to_exclude = flags_to_exclude;

% Define additional parameters for filtering data
novel_range = [];  % Not used in this version of the code, could be for a novel stimulus range
session_decoding_accuracy_thr = 5;  % Threshold for session decoding accuracy
must_include_all_flags = 1;  % 1 means all flags must be present for inclusion

% Rat names for easier identification in the results
Rat_Names = {'Clover','Bo','MEC1','Bolt','Dash','CM1','Janni','Harpy','Imp','W18','W19','Billy3','Curly2','Goethe2'};

% Categorize rats into experimental groups
experimental_rats = [1 2 3 6];  % Optogenetics with Jaws
control_rats1 = [4 5];  % GFP only (control group)
control_rats2 = [7 8 9 10 11 12 13 14];  % No laser or injections (control group)

% Define parameters for binning the time and angle data
binSize = 2;  % Bin size for time-based analysis (e.g., in seconds)
data_tbl = table();  % Initialize an empty table to store the event data

% Initialize session_id to keep track of the total number of valid sessions
session_id = 0;

% Loop over the rats and their respective data files
for rat_num = 1:length(rats)
    rat = rats(rat_num);  % Select the current rat
    load_open_field_session_list  % Load the session list for this rat

        % Loop over each day's data files
    for day = 1:length(dayFiles)
        display([Rat_Names{rat},' Day ', dayFiles{day}])  % Display the current rat and day
        cd(fullfile(directory,dayFiles{day}))  % Change to the directory of the day's data

        load Experiment_Information  % Load experimental information for the day

        
        % Identify the segments from the day's data that are relevant
        sessions_that_meet_criterion_day = check_for_relevant_sessions(flags,1,Experiment_Information,must_include_all_flags,1,novel_range);
        % Pull out all unique run sessions in this day:
        run_segments = [];
        for i = 1:length(Experiment_Information.Segments)
            if ismember(11,Experiment_Information.Segments(i).Flags) % If the segment meets the 'run' criteria
                run_segments = [run_segments; i];  % Add the segment index
            end
        end

      % If any valid sessions were found, proceed to analyze the data
        if ~isempty(sessions_that_meet_criterion_day)
            for session_count = 1:length(sessions_that_meet_criterion_day)
                session_id = session_id + 1;
                sub_session_num = sessions_that_meet_criterion_day(session_count); % Get the session number

                run_segment_num = find(run_segments == sub_session_num); % Find the relevant segment number
                sessionNum = sub_session_num; % Assign session number
                sessionNum_decoder = Experiment_Information.Segments(sessionNum).Decoder; % Get the decoder used
                sessionDecoding_error = Experiment_Information.Segments(sessionNum).decodingError;% Get the decoding error for the session
     
                % If the decoding error is below the threshold, proceed with analysis
                if sessionDecoding_error < session_decoding_accuracy_thr

                   load angles_between_replays.mat % Load data on the angles between replays
                    % Create a table to store the session data
                   t_sub = table();
                   t_sub.angles = angles_between_replays(run_segment_num).all_angles;
                   t_sub.comparison_inds = angles_between_replays(run_segment_num).all_comparison_inds;
                   t_sub.time_difference = angles_between_replays(run_segment_num).all_time_differences;
                   t_sub.time_into_stopping_period = angles_between_replays(run_segment_num).all_times_into_stopping_period;
                   t_sub.dispersion = angles_between_replays(run_segment_num).all_dispersion;
                   t_sub.start_distance_from_rat = angles_between_replays(run_segment_num).all_start_distance_from_rat;
                     % Append the session data to the main table
                    data_tbl = [data_tbl; t_sub];

                 
                end
            end
        end
    end
end
%%
% Filter the events based on certain conditions
events = data_tbl;

% Define thresholds for filtering events based on spatial and temporal criteria
start_distance_from_rat_thr = 7.5;  % Threshold for the start distance from the rat
dispersionThr = 0;  % Threshold for dispersion
max_time_difference = 4;  % Maximum allowed time difference between replays
restricted_stop_times = [1 inf];  % Restrict the time of the stopping period
ring_range = [7 90];  % Define the angular range for analysis

% Take the mean absolute angles within the range of interest
events.angles = nanmean(abs(events.angles(:,ring_range(1):ring_range(end))),2);

% Apply filtering conditions to select relevant events
sub_inds = find(events.time_into_stopping_period(:,1)>= restricted_stop_times(1) & events.time_into_stopping_period(:,1) <= restricted_stop_times(2) & ...
    events.time_into_stopping_period(:,2)>= restricted_stop_times(1) & events.time_into_stopping_period(:,2) <= restricted_stop_times(2) & ...
    events.dispersion(:,1)>dispersionThr & events.dispersion(:,2)>dispersionThr & ...
    events.time_difference < max_time_difference & events.start_distance_from_rat(:,1) < start_distance_from_rat_thr & events.start_distance_from_rat(:,2) < start_distance_from_rat_thr);

% Filter the events based on the selected indices
events = events(sub_inds,:);


%%
windowSize = 0.4;
windowShift = 0.1;
start_time = 0;
end_time = 4;
bin_start = (start_time:windowShift:(end_time-windowSize))';
bin_end = bin_start + windowSize;
bin_edges = [bin_start bin_end];
time_bin_centers = mean(bin_edges,2);
time_angle_hist_mean = nan(length(bin_edges),1);
time_angle_hist_sem = nan(length(bin_edges),1);

% Perform bootstrapping to estimate the mean and SEM of angles within each time window
num_boots=1000;
bootstrapped_time_angle_hist = nan(length(bin_edges),num_boots);
data_n = nan(length(bin_edges),1);
angles_binned = cell(length(bin_edges),1);

for i = 1:length(bin_edges)
    angles = events.angles(events.time_difference >= bin_edges(i,1) & events.time_difference < bin_edges(i,2));
    angles(isnan(angles))=[];
    data_n(i) = height(angles);
    angles_binned{i} = angles;
    time_angle_hist_mean(i) = nanmean(angles);
    time_angle_hist_sem(i) = nanstd(angles)./sqrt(data_n(i));

    % Perform bootstrapping by resampling angles with replacement
    for j = 1:num_boots
        rand_inds = randi(length(angles),[length(angles),1]);
        bootstrapped_time_angle_hist(i,j) = nanmean(angles(rand_inds));
    end
end

% Perform shuffling to generate a distribution of time angle histograms under the null hypothesis
num_shuffles = 100;
shuffled_time_angle_hist_mean = nan(num_shuffles,length(bin_edges));
for shuffle = 1:num_shuffles
    shuffled_events = events;
    shuffled_events.time_difference =  events.time_difference(randperm(height(events),height(events)));
    for i = 1:length(bin_edges)
        shuffled_angles = shuffled_events.angles(shuffled_events.time_difference >= bin_edges(i,1) & shuffled_events.time_difference < bin_edges(i,2));
        shuffled_angles(isnan(shuffled_angles))=[];
        shuffled_time_angle_hist_mean(shuffle,i) = nanmean(shuffled_angles);
    end
end

% Calculate 95% confidence intervals for the shuffled data
low = quantile(shuffled_time_angle_hist_mean, 0.025, 1);  % 2.5% quantile
high = quantile(shuffled_time_angle_hist_mean, 0.975, 1);  % 97.5% quantile

% Compute p-values by comparing the observed data to the shuffled data
pvals = nan(length(time_bin_centers),1);
for i = 1:length(time_bin_centers)
  right_side = sum(shuffled_time_angle_hist_mean(:,i) > time_angle_hist_mean(i));
  left_side = sum(shuffled_time_angle_hist_mean(:,i) < time_angle_hist_mean(i));
  count = min(right_side,left_side);
  pvals(i) = (count+1)/(num_shuffles+1);
end


sig_bins = find(time_angle_hist_mean' > high | time_angle_hist_mean' < low);

min_data_n = min(data_n);
max_data_n = max(data_n);
total_data_n_unbinned = sum(~isnan(events.angles));
%%

% Code for figure creation follows...

figure('Position',[2313 876 175 125])
plot(time_bin_centers,time_angle_hist_mean,'k','LineWidth',2)
xlabel('Time between replays (s)');
ylabel('|Displacement|')
box off
xticks(0:1:end_time)
hold on
xtickangle(0)
ylim([75 95])
yticks([75 85 95]);
yticklabels({})
h1 = plot(time_bin_centers,low,'--k'); hold on; h1.MarkerSize = 2;
h2 = plot(time_bin_centers,high,'--k'); hold on; h2.MarkerSize = 2; 
h = shadedErrorBar(time_bin_centers,time_angle_hist_mean,time_angle_hist_sem); hold on
h.edge(1).Color = 'none'; h.edge(2).Color = 'none';
plot(time_bin_centers(sig_bins),95*ones(length(sig_bins),1),'.k')
%  shadedErrorBar(time_bin_centers,shuffled_time_angle_hist_mean,shuffled_time_angle_hist_sem);
axis square
set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures';
saveas(gcf,fullfile(fig_path,'Replays avoid repeating themselves open field'),'jpg')
saveas(gcf,fullfile(fig_path,'Replays avoid repeating themselves open field'),'pdf')


% figure()
% plot(events.time_difference,events.angles,'.k')
% [rho,p] = nancorr(events.time_difference,events.angles)
% 
% x = events.time_difference; y = events.angles;
% x(isnan(y)) = []; y(isnan(y))= [];
% 
% g = fittype('a*exp(-b*x)+c');
% f0 = fit(x,y,g,'StartPoint',[[ones(size(x)),exp(x)]\y; 1]);
% 
% tbl = table(x(:),y(:));
% modelfun = @(b,x) b(1)*exp(-b(2)*x(:,1))+b(3);
% aGuessed = 30; bGuessed = 2.5; cGuessed = 90;
% beta0 = [aGuessed, bGuessed, cGuessed];
% md1 = fitnlm(tbl, modelfun, beta0);
% 
% coefficients = md1.Coefficients{:,'Estimate'};
% yFitted = coefficients(1)*exp(-coefficients(2)*linspace(0,max_time_difference,50)) + coefficients(3);
% hold on
% plot(linspace(0,max_time_difference,50),yFitted,'r-','LineWidth',4);
% 
% % xx = linspace(0,max_time_difference,50);
% % hold on;
% % plot(xx,f0(xx),'r-','linewidth',4)
% xlabel('Time between events, s')
% ylabel('Angle between events, deg')
% 
% figure()
% f = fit(x,y,'exp2');
% plot(f,x,y)
% 
% figure()
% hist(events.angles,100)

% windowSize = 0.25;
% windowShift = 0.05;
% start_time = 0;
% end_time = 4;
% 
% angle_bins = linspace(0, 180, 37);
% angle_bin_centers = mean([angle_bins(1:end-1)' angle_bins(2:end)'],2);

% time_angle_hist_mean = nan(length(bin_edges)-1,length(angle_bins)-1);
% for i = 1:length(bin_edges)-1
%     for j = 1:length(angle_bins)-1
%      time_angle_hist_mean(i,j) = sum(events.time_difference >= bin_edges(i) & events.time_difference < bin_edges(i+1) & ... 
%          events.angles >= angle_bins(j) & events.angles < angle_bins(j+1));
%     end
% end
% 
% figure()
% imagesc(time_angle_hist_mean');
% set(gca,'YDir','normal')
% yticks(0.5:2:36.5)
% yticklabels(num2cell(angle_bins(1:2:end)'))
% xticks(linspace(0.5,length(time_bin_centers)-0.5,9))
% xticklabels(num2cell(linspace(bin_edges(1),bin_edges(end),9)))
% 
% 
% figure()
% imagesc((time_angle_hist_mean./(sum(time_angle_hist_mean,2)))')
% set(gca,'YDir','normal')
% yticks(0.5:2:36.5)
% yticklabels(num2cell(angle_bins(1:2:end)'))
% xticks(linspace(0.5,length(time_bin_centers)-0.5,9))
% xticklabels(num2cell(linspace(bin_edges(1),bin_edges(end),9)))
% caxis([0 0.15])
% colormap(parula)