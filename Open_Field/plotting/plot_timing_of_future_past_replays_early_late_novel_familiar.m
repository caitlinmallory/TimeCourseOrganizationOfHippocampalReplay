% Figure 2M: Analyzing and visualizing the timing of future and past replays
% in both early and late trials, as well as in novel and familiar conditions.

% Define trial ranges: early trials (1-10) and late trials (30-40)
trial_range = [1 10; 30 40];  

% Define colors for different conditions: 
% Black for all replays, Green for future replays, Purple for past replays
color_replay = [0 0 0];
color_future = [.4660 0.6740 0.1880];  % RGB for future replay
color_past = [0.4940 0.1840 0.5560];    % RGB for past replay
colors = [color_future; color_past];    % Combine future and past colors

% Extract data for early trials (1-10) for prospective (future) and 
% retrospective (past) replays, within a 10-second window after drink onset
prospective_early = replay.time_since_real_drink_onset(replay.drink_period_number >= trial_range(1,1) & replay.drink_period_number <= trial_range(1,2) & ...
    replay.future == 1 & replay.time_since_real_drink_onset >= 0 & replay.time_since_real_drink_onset <= 10);

retrospective_early = replay.time_since_real_drink_onset(replay.drink_period_number >= trial_range(1,1) & replay.drink_period_number <= trial_range(1,2) & ...
    replay.past == 1 & replay.time_since_real_drink_onset >= 0 & replay.time_since_real_drink_onset <= 10);

% Extract data for late trials (30-40) for prospective and retrospective replays
prospective_late = replay.time_since_real_drink_onset(replay.drink_period_number >= trial_range(2,1) & replay.drink_period_number <= trial_range(2,2) & ...
    replay.future == 1 & replay.time_since_real_drink_onset >= 0 & replay.time_since_real_drink_onset <= 10);

retrospective_late = replay.time_since_real_drink_onset(replay.drink_period_number >= trial_range(2,1) & replay.drink_period_number <= trial_range(2,2) & ...
    replay.past == 1 & replay.time_since_real_drink_onset >= 0 & replay.time_since_real_drink_onset <= 10);

% Combine the extracted data for early and late trials into one vector
data = [prospective_early; retrospective_early; prospective_late; retrospective_late];

% Create a table for the data with columns for the time since reward zone entry,
% direction (prospective or retrospective), and early/late trial categorization
tbl = table();
tbl.time_since_reward_zone_entry = data;  % Time since reward zone entry
tbl.direction = [ones(size(prospective_early)); 2*ones(size(retrospective_early)); ones(size(prospective_late)); 2*ones(size(retrospective_late))]; % 1=prospective, 2=retrospective
tbl.early_late = [ones(size(prospective_early)); ones(size(retrospective_early)); 2*ones(size(prospective_late)); 2*ones(size(retrospective_late))]; % 1=early, 2=late

% Create a boxplot for the time since reward zone entry for early/late trials, 
% grouping by prospective and retrospective replays (color by direction)
figure('Position', [1335 707 100 100])
b = boxchart(tbl.early_late, tbl.time_since_reward_zone_entry, 'GroupByColor', tbl.direction);
colororder(colors)  % Apply the custom colors for future/past replays
xticks([1 2])  % Set x-ticks for early and late trial groups
xticklabels({})  % Remove x-tick labels
xlim([0 3])  % Set x-axis limits
yticks([0:2:10])  % Set y-ticks from 0 to 10 seconds with step size 2
ylabel('Time since arrival (s)')  % Label for y-axis

% Perform ranksum tests to compare prospective and retrospective replays
% in both early and late trials
[p_early, h, z_early] = ranksum(prospective_early, retrospective_early);  % Early trial comparison
[p_late, h, z_late] = ranksum(prospective_late, retrospective_late);      % Late trial comparison

% Set figure properties and save as images (JPEG and PDF)
set(gcf, 'Color', 'white', 'Renderer', 'painters', 'PaperPositionMode', 'auto');
saveas(gcf, fullfile(fig_path, 'timing_of_future_past_replays_early_or_late_trials'), 'jpeg');
saveas(gcf, fullfile(fig_path, 'timing_of_future_past_replays_early_or_late_trials'), 'pdf');

% Extract data for novel and familiar trials based on session flags (10.1 represents novel trials)
prospective_novel = replay.time_since_real_drink_onset(cellfun(@(x) ismember(10.1, x), replay.session_flags) & ...
    replay.future == 1 & replay.time_since_real_drink_onset >= 0 & replay.time_since_real_drink_onset <= 10);

retrospective_novel = replay.time_since_real_drink_onset(cellfun(@(x) ismember(10.1, x), replay.session_flags) & ...
    replay.past == 1 & replay.time_since_real_drink_onset >= 0 & replay.time_since_real_drink_onset <= 10);

% Extract data for familiar trials (where session flag is not 10.1)
prospective_familiar = replay.time_since_real_drink_onset(~cellfun(@(x) ismember(10.1, x), replay.session_flags) & ...
    replay.future == 1 & replay.time_since_real_drink_onset >= 0 & replay.time_since_real_drink_onset <= 10);

retrospective_familiar = replay.time_since_real_drink_onset(~cellfun(@(x) ismember(10.1, x), replay.session_flags) & ...
    replay.past == 1 & replay.time_since_real_drink_onset >= 0 & replay.time_since_real_drink_onset <= 10);

% Combine the novel and familiar trial data into one vector
data = [prospective_novel; retrospective_novel; prospective_familiar; retrospective_familiar];

% Create labels for the trial categories: 1=novel prospective, 2=novel retrospective, 
% 3=familiar prospective, 4=familiar retrospective
labels = [ones(size(prospective_novel)); 2*ones(size(retrospective_novel)); 3*ones(size(prospective_familiar)); 4*ones(size(retrospective_familiar))];

% Create a new table for novel and familiar trials with columns for time since reward zone entry,
% direction (prospective or retrospective), and novel/familiar categorization
tbl = table();
tbl.time_since_reward_zone_entry = data;
tbl.direction = [ones(size(prospective_novel)); 2*ones(size(retrospective_novel)); ones(size(prospective_familiar)); 2*ones(size(retrospective_familiar))];
tbl.novel_familiar = [ones(size(prospective_novel)); ones(size(retrospective_novel)); 2*ones(size(prospective_familiar)); 2*ones(size(retrospective_familiar))];

% Create a boxplot for novel and familiar trials, grouping by prospective and retrospective replays
figure('Position', [1335 707 100 100])
b = boxchart(tbl.novel_familiar, tbl.time_since_reward_zone_entry, 'GroupByColor', tbl.direction);
colororder(colors)  % Apply the custom colors for future/past replays
xticks([1 2])  % Set x-ticks for novel and familiar trial groups
xticklabels({})  % Remove x-tick labels
xlim([0 3])  % Set x-axis limits
yticks([0:2:10])  % Set y-ticks from 0 to 10 seconds with step size 2
ylabel('Time since arrival (s)')  % Label for y-axis

% Perform ranksum tests to compare prospective and retrospective replays
% in both novel and familiar trials
[p_novel, ~, z_novel] = ranksum(prospective_novel, retrospective_novel);  % Novel trial comparison
[p_familiar, ~, z_familiar] = ranksum(prospective_familiar, retrospective_familiar);  % Familiar trial comparison

% Set figure properties and save as images (JPEG and PDF)
set(gcf, 'Color', 'white', 'Renderer', 'painters', 'PaperPositionMode', 'auto');
saveas(gcf, fullfile(fig_path, 'timing_of_future_or_past_replays_novel_or_familiar_trials'), 'jpeg');
saveas(gcf, fullfile(fig_path, 'timing_of_future_or_past_replays_novel_or_familiar_trials'), 'pdf');
