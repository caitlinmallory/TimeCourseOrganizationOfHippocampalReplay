rats = [6 ] ;
fig_path = '/home/caitlin/Data/Processed_Data/Manuscript_Figures';
max_trial_number = inf;
max_latency = inf;
min_velocity = 0;
session_decoding_accuracy_thr = inf;
novel_range = [10.1 10.2];
thrForCategorization_include = 20;
flags_to_include = [11];
% 18 = dreadds; 19 = saline; 13 = laser off; 14 = laser on; what is 20?
flags_to_exclude = [100 20 19 13]; % MEC inactive
flags_to_exclude = [100 20 18 14]; % MEC active
flags_to_exclude = [100 20 101]; % excluded sessions with poor numbers of cells.
flags.flags_to_include = flags_to_include;
flags.flags_to_exclude = flags_to_exclude;
must_include_all_flags = 1; % 1: requires that all flags in Flags_to_include are present. 0: requires that at least one flag in Flags_to_include are present.

Rat_Names = {'Clover','Bo','MEC1','Bolt','Dash','CM1','Janni','Harpy','Imp','W18','W19','Billy3','Curly2','Goethe2'};
experimental_rats = [1 2 3 6]; % opto; with Jaws
control_rats1 = [4 5 ]; % opto, gfp only
control_rats2 = [7 8 9 10 11 12 13 14]; % no laser, no injections
binSize = 2;
data_tbl = table();

session_id = 0;
data_tbl = table();


for rat_num = 1:length(rats)
    rat = rats(rat_num);
    load_open_field_session_list_for_behavior

    for day = 1:length(dayFiles)
        display([Rat_Names{rat},' Day ', dayFiles{day}])
        cd(fullfile(directory,dayFiles{day}))

        load Experiment_Information
    
        run_segments = [];
        for i = 1:length(Experiment_Information.Segments)
            if ismember(11,Experiment_Information.Segments(i).Flags)
                run_segments = [run_segments; i]
            end
        end

        sessions_that_meet_criterion_day = check_for_relevant_sessions(flags,1,Experiment_Information,must_include_all_flags,1,novel_range);

        if ~isempty(sessions_that_meet_criterion_day)

            for session_count = 1:length(sessions_that_meet_criterion_day)
                %for session_count = 1
                session_id = session_id + 1;
                sessionNum = sessions_that_meet_criterion_day(session_count);
                sessionNum_decoder = Experiment_Information.Segments(sessionNum).Decoder;
                
                if ismember(14,Experiment_Information.Segments(sessionNum).Flags) | ismember(18,Experiment_Information.Segments(sessionNum).Flags)
                    laser_state = 1;
                else
                    laser_state = 0;
                end

                load Behavior_Analysis.mat
                row = find(run_segments == sessionNum);
                home_well = mode(Behavior_Analysis(row).goal_well_list(:,1));
                behavior_data = table();
                behavior_data.trial_number = [1:length(Behavior_Analysis(row).goal_well_list)]';
                
                behavior_data.home_trial = zeros(length(Behavior_Analysis(row).goal_well_list),1);
                behavior_data.home_trial(Behavior_Analysis(row).goal_well_list==home_well)=1;
                behavior_data.home_trial_number = cumsum(Behavior_Analysis(row).goal_well_list(:,1)==home_well);
                behavior_data.away_trial_number = cumsum(Behavior_Analysis(row).goal_well_list(:,1)~=home_well) ;
                behavior_data.latency = Behavior_Analysis(row).latencies;
                behavior_data.path_length = Behavior_Analysis(row).path_length;
                behavior_data.path_velocity = Behavior_Analysis(row).path_velocity;
                behavior_data.rat = rat*ones(height(behavior_data),1);
                behavior_data.session_id = session_id*ones(height(behavior_data),1);
                behavior_data.laser_state = laser_state.*ones(height(behavior_data),1);
                data_tbl = [data_tbl; behavior_data];
            end
        end
    end
end
%%

combined_behavior_data = data_tbl(data_tbl.trial_number <= max_trial_number & data_tbl.latency <= max_latency ...
    & data_tbl.path_velocity>min_velocity,:);

home_latency = combined_behavior_data.latency(combined_behavior_data.home_trial==1);
away_latency = combined_behavior_data.latency(combined_behavior_data.home_trial==0);
home_latency_0 = combined_behavior_data.latency(combined_behavior_data.home_trial==1 & combined_behavior_data.laser_state==0);
away_latency_0 = combined_behavior_data.latency(combined_behavior_data.home_trial==0 & combined_behavior_data.laser_state==0);
home_latency_1 = combined_behavior_data.latency(combined_behavior_data.home_trial==1 & combined_behavior_data.laser_state==1);
away_latency_1= combined_behavior_data.latency(combined_behavior_data.home_trial==0 & combined_behavior_data.laser_state==1);

home_path = combined_behavior_data.path_length(combined_behavior_data.home_trial==1);
away_path = combined_behavior_data.path_length(combined_behavior_data.home_trial==0);
home_path_0 = combined_behavior_data.path_length(combined_behavior_data.home_trial==1  & combined_behavior_data.laser_state==0);
away_path_0 = combined_behavior_data.path_length(combined_behavior_data.home_trial==0  & combined_behavior_data.laser_state==0);
home_path_1 = combined_behavior_data.path_length(combined_behavior_data.home_trial==1  & combined_behavior_data.laser_state==1);
away_path_1 = combined_behavior_data.path_length(combined_behavior_data.home_trial==0  & combined_behavior_data.laser_state==1);

home_path_velocity = combined_behavior_data.path_velocity(combined_behavior_data.home_trial==1);
away_path_velocity = combined_behavior_data.path_velocity(combined_behavior_data.home_trial==0);
home_path_velocity_0 = combined_behavior_data.path_velocity(combined_behavior_data.home_trial==1  & combined_behavior_data.laser_state==0);
away_path_velocity_0 = combined_behavior_data.path_velocity(combined_behavior_data.home_trial==0  & combined_behavior_data.laser_state==0);
home_path_velocity_1 = combined_behavior_data.path_velocity(combined_behavior_data.home_trial==1  & combined_behavior_data.laser_state==1);
away_path_velocity_1 = combined_behavior_data.path_velocity(combined_behavior_data.home_trial==0  & combined_behavior_data.laser_state==1);


colors = [0 0 1; 0 0 0];
transparency_pcnt = 0.5;
colors_2 = [[1 - transparency_pcnt*(1-colors(1,1)) 1 - transparency_pcnt*(1-colors(1,2)) 1 - transparency_pcnt*(1-colors(1,3))]; ...
    [1-transparency_pcnt*(1-colors(2,1)) 1-transparency_pcnt*(1-colors(2,2)) 1-transparency_pcnt*(1-colors(2,3))]];
%% Plot latency for laser off or on home and away trials:
data = [mean(home_latency_0) mean(home_latency_1) mean(away_latency_0) mean(away_latency_1)];
err = [std(home_latency_0)./sqrt(length(home_latency_0))... 
   std(home_latency_1)./sqrt(length(home_latency_1))... 
   std(away_latency_0)./sqrt(length(away_latency_0))... 
   std(away_latency_1)./sqrt(length(away_latency_1))]; 
plot_fig(data,err,colors,colors_2)
set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
saveas(gcf,fullfile(fig_path,['trial_latency_' num2str(rats)]),'jpeg')
saveas(gcf,fullfile(fig_path,['trial_latency_' num2str(rats)]),'pdf')
%% Plot pathlength for laser off or on home and away trials:
data = [mean(home_path_0) mean(home_path_1) mean(away_path_0) mean(away_path_1)];
err = [std(home_path_0)./sqrt(length(home_path_0))... 
   std(home_path_1)./sqrt(length(home_path_1))... 
   std(away_path_0)./sqrt(length(away_path_0))... 
   std(away_path_1)./sqrt(length(away_path_1))]; 
plot_fig(data,err,colors,colors_2)
set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
saveas(gcf,fullfile(fig_path,['trial_pathlength_' num2str(rats)]),'jpeg')
saveas(gcf,fullfile(fig_path,['trial_pathlength_' num2str(rats)]),'pdf')
%% Plot velocity for laser off or on home and away trials:
data = [mean(home_path_velocity_0) mean(home_path_velocity_1) mean(away_path_velocity_0) mean(away_path_velocity_1)];
err = [std(home_path_velocity_0)./sqrt(length(home_path_velocity_0))... 
   std(home_path_velocity_1)./sqrt(length(home_path_velocity_1))... 
   std(away_path_velocity_0)./sqrt(length(away_path_velocity_0))... 
   std(away_path_velocity_1)./sqrt(length(away_path_velocity_1))]; 
plot_fig(data,err,colors,colors_2)
%%

trialOrder = ["Home","Away"];
namedTrialType = categorical(combined_behavior_data.home_trial,[1,0],trialOrder);
groupOrder = ["Experimental","Control"];
namedGroup = categorical(combined_behavior_data.laser_state, [1,0],groupOrder);
% Anova effects:
p_anova = anovan(combined_behavior_data.latency,{namedTrialType,namedGroup},3,3,{'namedTrialType';'namedGroup'});
p_anova = anovan(combined_behavior_data.path_length,{namedTrialType,namedGroup},3,3,{'namedTrialType';'namedGroup'});
p_anova = anovan(combined_behavior_data.path_velocity,{namedTrialType,namedGroup},3,3,{'namedTrialType';'namedGroup'});


ranksum(home_latency_0,home_latency_1)
ranksum(away_latency_0,away_latency_1)

%%


sessions = unique(combined_behavior_data.session_id);

[C,ia, ic] = unique(combined_behavior_data(:,{'rat','laser_state','session_id'})); 

for i = 1:height(C)
    sub_data = combined_behavior_data(combined_behavior_data.session_id == C.session_id(i),:);
    C.latency_sub(i) = nanmean(sub_data.latency(sub_data.home_trial==0))-nanmean(sub_data.latency(sub_data.home_trial==1));
    C.latency_div(i) = nanmean(sub_data.latency(sub_data.home_trial==0))/nanmean(sub_data.latency(sub_data.home_trial==1));
    C.pathlength_sub(i) = nanmean(sub_data.path_length(sub_data.home_trial==0))-nanmean(sub_data.path_length(sub_data.home_trial==1));
    C.pathlength_div(i) = nanmean(sub_data.path_length(sub_data.home_trial==0))/nanmean(sub_data.path_length(sub_data.home_trial==1));
end

latency_sub_0 = C.latency_sub(C.laser_state==0);
latency_sub_1 = C.latency_sub(C.laser_state==1);
[p,h,z] = ranksum(latency_sub_0,latency_sub_1)
data = [nanmean(latency_sub_0) nanmean(latency_sub_1)]; 
err = [nanstd(latency_sub_0)./sqrt(length(latency_sub_0)) nanstd(latency_sub_1)./sqrt(length(latency_sub_1))];
plot_fig2(data,err,colors,colors_2)
set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
saveas(gcf,fullfile(fig_path,['session_latency_home_away_diff_sub_' num2str(rats)]),'jpeg')
saveas(gcf,fullfile(fig_path,['session_latency_home_away_diff_sub_' num2str(rats)]),'pdf')

% latency_div_0 = C.latency_div(C.laser_state==0);
% latency_div_1 = C.latency_div(C.laser_state==1);
% ranksum(latency_div_0,latency_div_1)
% data = [nanmean(latency_div_0) nanmean(latency_div_1)]; 
% err = [nanstd(latency_div_0)./sqrt(length(latency_div_0)) nanstd(latency_div_1)./sqrt(length(latency_div_1))];
% plot_fig2(data,err,colors,colors_2)

pathlength_sub_0 = C.pathlength_sub(C.laser_state==0);
pathlength_sub_1 = C.pathlength_sub(C.laser_state==1);
[p,h,z] = ranksum(pathlength_sub_0,pathlength_sub_1)
data = [nanmean(pathlength_sub_0) nanmean(pathlength_sub_1)]; 
err = [nanstd(pathlength_sub_0)./sqrt(length(pathlength_sub_0)) nanstd(pathlength_sub_1)./sqrt(length(pathlength_sub_1))];
plot_fig2(data,err,colors,colors_2)
set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
saveas(gcf,fullfile(fig_path,['session_pathlength_home_away_diff_sub_' num2str(rats)]),'jpeg')
saveas(gcf,fullfile(fig_path,['session_pathlength_home_away_diff_sub_' num2str(rats)]),'pdf')

% pathlength_div_0 = C.pathlength_div(C.laser_state==0);
% pathlength_div_1 = C.pathlength_div(C.laser_state==1);
% ranksum(pathlength_div_0,pathlength_div_1)
% data = [nanmean(pathlength_div_0) nanmean(pathlength_div_1)]; 
% err = [nanstd(pathlength_div_0)./sqrt(length(pathlength_div_0)) nanstd(pathlength_div_1)./sqrt(length(pathlength_div_1))];
% plot_fig2(data,err,colors,colors_2)



%%

function [] = plot_fig(data,err,colors,colors_2)
figure('Position',[1986 1051 100 100])
b1 = bar(1,data(1)); hold on;
e1 = errorbar(1,data(1),err(1),'k','linestyle','none')
e1.CapSize = 4;
hold on
b2 = bar(2,data(2)); hold on;
e2 = errorbar(2,data(2),err(2),'k','linestyle','none')
e2.CapSize = 4;
b3 = bar(4,data(3)); hold on;
e3 = errorbar(4,data(3),err(3),'k','linestyle','none')
e3.CapSize = 4;
b4 = bar(5,data(4)); hold on;
e4 = errorbar(5,data(4),err(4),'k','linestyle','none')
e4.CapSize = 4;
b1.FaceColor = colors(1,:);
b1.EdgeColor = 'none';
b2.FaceColor = colors_2(1,:);
b2.EdgeColor = 'none';
b3.FaceColor = colors(2,:);
b3.EdgeColor = 'none';
b4.FaceColor = colors_2(2,:);
b4.EdgeColor = 'none';
box off
xticks([1.5 4.5])
xticklabels({})
end


function [] = plot_fig2(data,err,colors,colors_2)
figure('Position',[600 600 100 100])
b1 = bar(1,data(1)); hold on;
e1 = errorbar(1,data(1),err(1),'k','linestyle','none');
e1.CapSize = 4;
hold on
b2 = bar(2,data(2)); hold on;
e2 = errorbar(2,data(2),err(2),'k','linestyle','none');
e2.CapSize = 4;
b1.FaceColor = colors(1,:);
b1.EdgeColor = 'none';
b2.FaceColor = colors_2(1,:);
b2.EdgeColor = 'none';

box off
xticks([1.5 4.5])
xticklabels({})
end