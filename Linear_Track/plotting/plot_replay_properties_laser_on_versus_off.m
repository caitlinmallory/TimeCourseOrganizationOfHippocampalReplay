%%
% setup_properties_table % contains a list of all properties you might want to look at/plot
load '/home/caitlin/Data/Processed_Data/properties_table.mat'

control_experimental_animals_to_include = [0]; % set to 1 to look at experimental (Jaws) animals, 0 to look at control (GFP) animals.
upper_limit_time_since_reward_zone_entry = inf; % inf is most permissive
lower_limit_time_since_reward_zone_entry = 0; % -inf is most permissive

% Forward congruent versus reverse congruent:
events = t_replay(t_replay.congruent_with_rat_location==1,:);
events.replay_class(events.direction==2 & events.congruent_with_rat_location==1) = 2;
events.replay_class(events.direction==1 & events.congruent_with_rat_location==1) = 1;

events.group(events.laser_state==0 & events.replay_class==2) = 1; % reverse 
events.group(events.laser_state==1 & events.replay_class==2) = 2; % reverse 
events.group(events.laser_state==0 & events.replay_class==1) = 3; % forward 
events.group(events.laser_state==1 & events.replay_class==1) = 4; % forward 

events = events(ismember(events.control_experimental_flag,control_experimental_animals_to_include) & events.time_since_reward_zone_entry > lower_limit_time_since_reward_zone_entry & events.time_since_reward_zone_entry < upper_limit_time_since_reward_zone_entry,:);

plot_individual_graphs = 0;
plot_boxplots=0;
plot_barplots=0;
plot_pdfs=0;
plot_cdfs=0;

fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures';

plot_supp_fig_8_DE

