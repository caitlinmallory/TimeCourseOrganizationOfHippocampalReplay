
names = {...

    'posterior_spread_full'; ...
    'posterior_spread_cropped'; ...
    'replay_spikeDensity_power'; ...% 1
    'replay_spikeDensity_power_pyr';...
    'mean_replay_spikeDensity_power';...
    'replay_ripple_power';...% 2
    'mean_replay_ripple_power';...
    'weighted_r';... % 3

    'mean_num_spikes_participating_excitatory_cells';... % 4
    'mean_fr_participating_excitatory_cells';... % 5
    'mean_num_spikes_all_excitatory_cells';... % 6
    'mean_fr_all_excitatory_cells';... % 7
    'fraction_of_excitatory_cells_participating';... % 8

    'mean_num_spikes_participating_inhibitory_cells';... % 9
    'mean_fr_participating_inhibitory_cells';... % 10
    'mean_num_spikes_all_inhibitory_cells';... % 11
    'mean_fr_all_inhibitory_cells';... % 12
    'fraction_of_inhibitory_cells_participating';... % 13

    'mean_num_spikes_participating_unimodal_cells';... % 14
    'mean_fr_participating_unimodal_cells';... % 15
    'mean_num_spikes_all_unimodal_cells';... % 16
    'mean_fr_all_unimodal_cells';... % 17
    'fraction_of_unimodal_cells_participating';... % 18

    'mean_num_spikes_participating_bimodal_cells';... % 19
    'mean_fr_participating_bimodal_cells';... % 20
    'mean_num_spikes_all_bimodal_cells';... % 21
    'mean_fr_all_bimodal_cells';... % 22
    'fraction_of_bimodal_cells_participating';... % 23

    'duration';... % 24
    'duration_og';
    'sharpness';... % 25
    'range';...
    'range_normalized';... % 26
    'range_over_time';... % 27
    'end_distance_from_rat';... % 28
    'time_in_session';... % 29
    'time_since_reward_zone_entry'; ... % 30
    'ratSpeed'; ...

    'ave_sharpness_at_peak'; ...
    'coverage_wu';...
    'coverage_wu_range';...
    'cum_coverage';...
    'max_jump_distance';...
    'ave_jump_distance';...
    'start_distance_from_rat'...
    
    };



titles = {...
    'Spread full';...
    'Spread cropped'; ...
    'Spike density'; ... % 1
    'Spike density pyr'; ...
    'Mean spike density'; ...
    'Ripple power'; ...% 2
    'Mean ripple power'; ...
    'Weighted R'; ... % 3

    'Active E-cells num spikes'; ... % 4
    'Active E-cells fr'; ... % 5
    'Ave. num of spikes'; ... % 6
    'All E-cells fr'; ... % 7
    '% E-cells'; ... % 8


    'Active I-cells num spikes'; ... % 9
    'Active I-cells fr'; ... % 10
    'All I-cells num spikes'; ... % 11
    'All I-cells fr'; ... % 12
    '% I-cells'; ... % 13

    'Active U-cells num spikes'; ... % 14
    'Active U-cells fr'; ... % 15
    'All U-cells num spikes'; ... % 16
    'All U-cells fr'; ... % 17
    '% U-cells'; ... % 18


    'Active B-cells num spikes'; ... % 19
    'Active B-cells fr'; ... % 20
    'All B-cells num spikes'; ... % 21
    'All B-cells fr'; ... % 22
    '% B-cells'; ... % 23


    'Duration'; ... % 24
    'Duration of C.E.'
    'Sharpness'; ... % 25
    'Range'; ... range
    'Range'; ... % range normalized
    'Slope'; ... % 27
    'End distance from rat'; ... % 28
    'Time in session'; ... % 29
    'Time since reward onset'; ... % 30
    'Rat speed'; ... %31

    'Ave sharpness at peak'; ...
    'Coverage Wu';...
    'Range Wu'; 
    'Cum. Coverage';...
    'Max jump distance';...
    'Ave jump distance';...
    'Start distance from rat'};


ylabels = {...

    'Spread'; ... % posterior spread full
    'Spread'; ... % posterior spread cropped
    'Z-score'; ... % 1 Spike density
    'Z-score'; ... % Spike density pyr
    'Z-score'; ... % Mean spike density
    'Z-score'; ... % 2 Ripple power
    'Z-score'; ... % Mean ripple power
    ''; ...        % 3 Weighted r

    'Num spikes'; ... % 4 active E cells. num spikes
    'Spikes/s'; ... % 5 active E cells fr
    'Num spikes'; ... % 6 all E cells num spikes
    'Spikes/s'; ... % 7 all E cells fr
    '%'; ... 8 percent participation of E cells


    'Num spikes'; ... % 9 active I cells num spikes
    'Spikes/s'; ... % 10 active I cells fr
    'Num spikes'; ... % 11 all I cells num spikes
    'Spikes/s'; ... % 12 all I cells fr
    '%'; % 13 percent participation of I cells

    'Num spikes'; ... % 14 active U cells num spikes
    'Spikes/s'; ... % 15 active U cells fr
    'Num spikes'; ... % 16 all U cells num spikes
    'Spikes/s'; ... % 17 all U cells fr.
    '%'; % 18 percent participation of U cells

    'Num spikes'; ... % 19 active B cells num spikes
    'Spikes/s'; ... % 20 active B cells fr
    'Num spikes'; ... % 21 all B cells num spikes
    'Spikes/s'; ... % 22 all B cells fr
    '%'; % 23 percent participation of B cells


    's'; ... % 24 duration
    's'; ... % duration of underlying CE
    ''; ... % 25 sharpness
    'cm'; ... range
    '% of track'; ... % 26 range normalized
    'm/s'; ... % 27 speed
    '% of track'; ... % 28 end distance from rat
    's'; ... % 22 time in session
    's'; ... % 30' time in stopping period
    'cm/s'; ... % 31 running speed
    '';...    'Ave sharpness at peak'; ...
    '%'; ... coverage_wu';...
    '%'; ... coverage wu;
    ''; ... 'cum_coverage';...
    ''; ... 'max jump distnace';
    ''; ... 'ave_jump_distance';...
    'bins'; ... 'start_distance_from_rat'};
    };

ylims = [...

    {[0.05 0.1; 0.05 0.1]}; % posterior spread full
    {[0.05 0.1; 0.05 0.1]}; % posterior spread cropped
    {[2.5 7.5; 2.5 7.5]};% 1  Spike density
    {[2.5 7.5; 2.5 7.5]};% Spike density pyr
    {[0 3.5; 0 3.5]}; % Mean spike density
    {[2 8; 2 8]};% 2 ripple power
    {[0 3; 0 3]}; % Mean ripple power
    {[0.7 0.9; 0.3 0.6]};% 3 weighted r

    {[7 11; 7 11]};% 4 active E cells num spikes
    {[40 90; 40 90]};% 5 active E cells fr
    {[0 5; 0 5]};% 6 all E cells num spikes
    {[10 30; 10 30]}; % 7 all E cells fr
    {[0.2 0.4; 0.2 0.4]}; %  E cells PP

    {[15 50; 15 50]};% 9 active I cells num spikes
    {[150 300; 150 300]}; % 10 active I cells fr
    {[15 50; 15 50]}; % 11 all I cells num spikes
    {[150 300; 150 300]}; % 12 all I cells fr
    {[0.8 1; 0.8 1]}; % 13 I E cells PP

    {[8 15; 8 15]};% 14 active U cells num spikes
    {[50 100; 50 100]};% 15 active U cells fr
    {[1 10; 1 10]};% 16 all U cells num spikes
    {[10 50; 10 50]}; % 17 all U cells fr
    {[0.1 0.4; 0.1 0.4]}; % 18 U cells PP

    {[6 15; 6 15]};% 19 active B cells num spikes
    {[40 100; 40 100]};% 20 active B cells fr
    {[2 12; 2 12]};% 21 all B cells num spikes
    {[10 40; 10 40]}; % 22 all B cells fr
    {[0.2 0.6; 0.2 0.6]}; % 23 B cells PP

    {[0.15 0.25; 0.15 0.25]}; % 24 duration
    {[0.15 0.25; 0.15 0.25]}; % 24 duration of underlying. C.e.

    {[0.2 0.5; 0.2 0.5]}; % 25  sharpness
    {[75 175; 75 175]}; % 26 range
    {[0.4 0.8; 0.2 0.4]}; % 26 range normalized    
    {[6 12; 4 8]}; % 27 slope
    {[0.4 1; 0.4 1]}; % 28 end distance from rat
    {[0 1000; 0 1000]}; % 29 time in session
    {[0 10; 0 10]}; % 30 time in stopping period
    {[0 5; 0 5]}; % 31 rat speed

    {[0 1; 0 1]}; %  'Ave sharpness at peak'; ...
    {[0.4 0.7; 0.2 0.4]}; % 'coverage_wu';...
    {[0.4 0.7; 0.2 0.4]}; % 'coverage_wu_range';...    
    {[0 260; 0 260]}; % 'cum_coverage';...
    {[0 0.4; 0 0.4]}; % 'max_jump_distance';...    
    {[0 0.4; 0 0.4]}; % 'ave_jump_distance';...
    {[0 30; 0 30]}]; %'start_distance_from_rat'};



hist_axis = [...

    {0:0.01:0.2}; % posterior spread full
    {0:0.01:0.2}; % posterior spread cropped
    {0:2:30}; %1 spike density
    {0:2:30}; %1 spike density pyr    
    {0:2:30}; % mean spike density
    {0:2:30}; % 2
    {0:2:30}; % mean ripple power
    {0.6:0.05:1};% 3

    {0:1:20};% 4
    {0:10:200}; % 5
    {0:1:5};% 6
    {0:10:100}; % 7
    {0:0.05:1};% 8

    {0:10:100};% 9
    {0:50:500}; % 10
    {0:10:100};% 11
    {0:50:500}; % 12
    {0:0.05:1};% 13

    {0:1:20}; % 14
    {0:10:200}; % 15
    {0:1:10};% 16
    {0:10:100}; %17
    {0:0.05:1};% 18

    {0:1:20}; % 19
    {0:10:200}; % 20
    {0:1:10};% 21
    {0:10:100}; %22
    {0:0.05:1};% 23


    {0:0.05:0.5};  % 24
    {0:0.05:0.5};  % 24

    {0:0.05:0.8}; % 25
    {0:10:100}
    {0.2:0.1:1}; % 26
    {0:2:30}; % 27
    {0:0.1:1}; % end distance from rat
    {0:300:2000}; % 29
    {0:2:60}; % 30
    {0:0.5:5};

    {0:0.1:1}; %  'Ave sharpness at peak'; ...
    {0:0.1:1}; % 'coverage_wu';...
    {0:0.1:1}; % 'coverage_wu_range';...    
    {0:10: 130}; % 'cum_coverage';...
     {0:0.05:0.4}; % 'ave_jump_distance';...   
    {0:0.05:0.4}; % 'ave_jump_distance';...
    {0:5:30}]; %'start_distance_from_rat'};


properties = table(names,titles,ylims,ylabels,hist_axis);