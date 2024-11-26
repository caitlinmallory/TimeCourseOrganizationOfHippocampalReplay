function [downsampleindices1, downsampleindices2] = downsample_match_time_in_stopping_period(x1,y1,x2,y2)



% matches two sessions based on time in session and time in stopping period

% x1 = t_replay.time_in_session(t_replay.laser_state==0);
% x2 = t_replay.time_in_session(t_replay.laser_state==1);
% 
% 
% y1 = t_replay.time_since_reward_zone_entry(t_replay.laser_state==0);
% y2 = t_replay.time_since_reward_zone_entry(t_replay.laser_state==1);

f = figure();
f.Position = [520 165 977 658];
subplot(2,2,1)
ecdf(x1);
hold on
ecdf(x2);
title('Original time into session')
legend({'laser off','laser on'})
text(200,0.8,['n = ' num2str(length(x1) + length(x2))])

subplot(2,2,2)
ecdf(y1);
hold on
ecdf(y2);
title('Original time into stopping period')

% num_session_time_bins = 10; %from Hardcastle et al (2017)
% num_time_into_stopping_period_bins = 30;
% xAxis = linspace(0,max([x1;x2]),num_session_time_bins+1);
% yAxis = linspace(0,max([y1;y2]),num_time_into_stopping_period_bins+1);
% 

% number.
% time_into_session_bin_size = 5*60; % 5 minute time bins
time_into_stopping_period_bin_size = 1; % 3 second time bins
num_time_into_stopping_period_bins = ceil(max([y1;y2])/time_into_stopping_period_bin_size);
% num_session_time_bins = ceil(max([x1;x2])/time_into_session_bin_size);
%num_session_time_bins = 1;


xAxis = linspace(0,max([y1;y2]),num_time_into_stopping_period_bins +1);


%initialize matrices
map1 = zeros(num_time_into_stopping_period_bins, 1);

downsampleindices1 = [];
downsampleindices2 = [];
%% Determine number of observations in each bin in session1
for i = 1:num_time_into_stopping_period_bins
   
        ind1 = find(y1 >= xAxis(i) & y1 < xAxis(i+1));
        map1(i)= numel(ind1);
        
        %same thing for map 2
        ind2 = find(y2 >= xAxis(i) & y2 < xAxis(i+1));

        % Number of observations to keep in each map
        numToKeep = min(numel(ind1),numel(ind2));
        
        downsampleindices1 = [downsampleindices1; datasample(ind1, numToKeep, 'Replace', false)];
        downsampleindices2 = [downsampleindices2; datasample(ind2, numToKeep, 'Replace', false)];
   
end


downsampleindices1 = sort(downsampleindices1);
downsampleindices2 = sort(downsampleindices2);


subplot(2,2,3)
ecdf(x1(downsampleindices1));
hold on
ecdf(x2(downsampleindices2));
title('Downsampled time into session')
legend({'laser off','laser on'})
text(200,0.8,['n = ' num2str(length(downsampleindices1) + length(downsampleindices2))])


subplot(2,2,4)
ecdf(y1(downsampleindices1));
hold on
ecdf(y2(downsampleindices2));
title('Downsampled time into stopping period')
text(25,0.4,[num2str(num_time_into_stopping_period_bins) ' stopping period bins'])




end