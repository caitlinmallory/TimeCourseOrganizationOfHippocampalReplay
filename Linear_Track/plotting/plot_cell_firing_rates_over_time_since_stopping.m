%%
define_by_place_fields=1;
define_by_place_fields_and_modulation=0;
define_by_modulation=0;
a = data_tbl.clusters_left_field_only & data_tbl.Well_Isolated==1; 
b = data_tbl.clusters_right_field_only & data_tbl.Well_Isolated==1; 
c = data_tbl.clusters_right_and_left_field & data_tbl.Well_Isolated==1; 

% define_by_place_fields=1;
% define_by_place_fields_and_modulation=0;
% define_by_modulation=0;
% a = data_tbl.clusters_left_field_only & data_tbl.left_com_spatial_firing_bin_location<0.33; 
% b = data_tbl.clusters_right_field_only & data_tbl.right_com_spatial_firing_bin_location>0.66;
% c = data_tbl.clusters_right_and_left_field;

% define_by_place_fields=0;
% define_by_place_fields_and_modulation=1;
% define_by_modulation=0;
% a = data_tbl.clusters_left_field_only & data_tbl.lr_firing_rate_modulation>0.5;
% b = data_tbl.clusters_right_field_only & data_tbl.lr_firing_rate_modulation>0.5;
% c = data_tbl.clusters_right_and_left_field;

% define_by_place_fields=0;
% define_by_place_fields_and_modulation=0;
% define_by_modulation=1;
% a = data_tbl.lr_firing_rate_modulation>0.5 & data_tbl.clusters_left_max_rate_higher;
% b = data_tbl.lr_firing_rate_modulation>0.5 & data_tbl.clusters_right_max_rate_higher;
% c = data_tbl.lr_firing_rate_modulation<=0.5;

bin_size = 0.025;
length_to_plot = 10;
bins_forward = length_to_plot/bin_size;
time_vec = bin_size/2:bin_size:(length_to_plot-bin_size/2);

%%
figure('Position',[147 649 800 250])
subplot(1,3,1)
plot(time_vec,nanmean([data_tbl.zscored_firing_rates_at_right(b,:); data_tbl.zscored_firing_rates_at_left(a,:)]),'b'); 
hold on;
plot(time_vec,nanmean([data_tbl.zscored_firing_rates_at_right(a,:); data_tbl.zscored_firing_rates_at_left(b,:)]),'r'); 
hold on;
plot(time_vec,nanmean([data_tbl.zscored_firing_rates_at_right(c,:); data_tbl.zscored_firing_rates_at_left(c,:)]),'m'); 
hold on;
title('firing rate during stopping period')
legend({'active on past lap','active on future lap','active on both laps'})

subplot(1,3,2)
plot(time_vec,nanmean([data_tbl.zscored_firing_rates_at_right_sde(b,:); data_tbl.zscored_firing_rates_at_left_sde(a,:)]),'b'); 
hold on;
plot(time_vec,nanmean([data_tbl.zscored_firing_rates_at_right_sde(a,:); data_tbl.zscored_firing_rates_at_left_sde(b,:)]),'r'); 
hold on;
plot(time_vec,nanmean([data_tbl.zscored_firing_rates_at_right_sde(c,:); data_tbl.zscored_firing_rates_at_left_sde(c,:)]),'m'); 
title('firing rate during SDEs')

subplot(1,3,3)
plot(time_vec,nanmean([data_tbl.zscored_firing_rates_at_right_replay(b,:); data_tbl.zscored_firing_rates_at_left_replay(a,:)]),'b'); 
hold on;
plot(time_vec,nanmean([data_tbl.zscored_firing_rates_at_right_replay(a,:); data_tbl.zscored_firing_rates_at_left_replay(b,:)]),'r'); 
hold on;
plot(time_vec,nanmean([data_tbl.zscored_firing_rates_at_right_replay(c,:); data_tbl.zscored_firing_rates_at_left_replay(c,:)]),'m'); 
hold on;
title('firing rate during replay')
set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
if define_by_place_fields==1
    saveas(gcf,'/home/caitlin/Data/Processed_Data/Manuscript_Figures/ZscoredFiringRateInReplays_UsingFields','png');
elseif define_by_modulation==1
saveas(gcf,'/home/caitlin/Data/Processed_Data/Manuscript_Figures/ZscoredFiringRateInReplays_UsingFiringRateModulation','png');
elseif define_by_place_fields_and_modulation==1
saveas(gcf,'/home/caitlin/Data/Processed_Data/Manuscript_Figures/ZscoredFiringRateInReplays_UsingFieldsAndFiringRateModulation','png');
end  
%%
figure('Position',[147 649 800 250])
subplot(1,3,1)
plot(time_vec,nanmean([data_tbl.firing_rates_at_right(b,:); data_tbl.firing_rates_at_left(a,:)]),'b'); 
hold on;
plot(time_vec,nanmean([data_tbl.firing_rates_at_right(a,:); data_tbl.firing_rates_at_left(b,:)]),'r'); 
hold on;
plot(time_vec,nanmean([data_tbl.firing_rates_at_right(c,:); data_tbl.firing_rates_at_left(c,:)]),'m'); 
hold on;
title('firing rate during stopping period')
xlim([0 10])

subplot(1,3,2)
plot(time_vec,nanmean([data_tbl.firing_rates_at_right_sde(b,:); data_tbl.firing_rates_at_left_sde(a,:)]),'b'); 
hold on;
plot(time_vec,nanmean([data_tbl.firing_rates_at_right_sde(a,:); data_tbl.firing_rates_at_left_sde(b,:)]),'r'); 
hold on;
plot(time_vec,nanmean([data_tbl.firing_rates_at_right_sde(c,:); data_tbl.firing_rates_at_left_sde(c,:)]),'m'); 
title('firing rate during SDEs')
xlim([0 10])

subplot(1,3,3)
plot(time_vec,nanmean([data_tbl.firing_rates_at_right_replay(b,:); data_tbl.firing_rates_at_left_replay(a,:)]),'b'); 
hold on;
plot(time_vec,nanmean([data_tbl.firing_rates_at_right_replay(a,:); data_tbl.firing_rates_at_left_replay(b,:)]),'r'); 
hold on;
plot(time_vec,nanmean([data_tbl.firing_rates_at_right_replay(c,:); data_tbl.firing_rates_at_left_replay(c,:)]),'m'); 
hold on;
title('firing rate during replay')
set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
if define_by_place_fields==1
    saveas(gcf,'/home/caitlin/Data/Processed_Data/Manuscript_Figures/FiringRateInReplays_UsingFields','png');
elseif define_by_modulation==1
saveas(gcf,'/home/caitlin/Data/Processed_Data/Manuscript_Figures/FiringRateInReplays_UsingFiringRateModulation','png');
elseif define_by_place_fields_and_modulation==1
saveas(gcf,'/home/caitlin/Data/Processed_Data/Manuscript_Figures/FiringRateInReplays_UsingFieldsAndFiringRateModulation','png');
end  
xlim([0 10])

%% For publication:
figure('Position',[100 100 300 300])
y1 = [data_tbl.firing_rates_at_right_replay(b,:); data_tbl.firing_rates_at_left_replay(a,:)];
y2 = [data_tbl.firing_rates_at_right_replay(a,:); data_tbl.firing_rates_at_left_replay(b,:)]; 
y3 = [data_tbl.firing_rates_at_right_replay(c,:); data_tbl.firing_rates_at_left_replay(c,:)]; 

y1_mean = nanmean(y1);
y1_sem = nanstd(y1)/sqrt(size(y1,1))
y2_mean = nanmean(y2);
y2_sem = nanstd(y2)/sqrt(size(y2,1))
y3_mean = nanmean(y3);
y3_sem = nanstd(y3)/sqrt(size(y3,1))

shadedErrorBar(time_vec,y1_mean,y1_sem,'lineProps','-b'); hold on;
shadedErrorBar(time_vec,y2_mean,y2_sem,'lineProps','-r'); hold on;
shadedErrorBar(time_vec,y3_mean,y3_sem,'lineProps','-m'); hold on;
ylim([0 4]);

[rho,p]=nancorr(time_vec',y1_mean')
[rho,p]=nancorr(time_vec',y2_mean')
[rho,p]=nancorr(time_vec',y3_mean')

early_past_cells=nanmean(y1(:,time_vec<=3),2)
late_past_cells=nanmean(y1(:,time_vec>3),2)
[p,h,z] = signrank(early_past_cells,late_past_cells)

early_future_cells=nanmean(y2(:,time_vec<=3),2)
late_future_cells=nanmean(y2(:,time_vec>3),2)
[p,h,z] = signrank(early_future_cells,late_future_cells)

early_bidirectional_cells=nanmean(y3(:,time_vec<=3),2)
late_bidirectional_cells=nanmean(y3(:,time_vec>3),2)
[p,h,z] = signrank(early_bidirectional_cells,late_bidirectional_cells)


title('firing rate during replay')
set(gcf, 'Color', 'white','Renderer','painters','PaperPositionMode','auto');
if define_by_place_fields==1
    saveas(gcf,'/home/caitlin/Data/Processed_Data/Manuscript_Figures/FiringRateInReplays_UsingFields','pdf');
elseif define_by_modulation==1
saveas(gcf,'/home/caitlin/Data/Processed_Data/Manuscript_Figures/FiringRateInReplays_UsingFiringRateModulation','pdf');
elseif define_by_place_fields_and_modulation==1
saveas(gcf,'/home/caitlin/Data/Processed_Data/Manuscript_Figures/FiringRateInReplays_UsingFieldsAndFiringRateModulation','pdf');
end  
xlim([0 10])

