function [boundaries_merged,lengths_merged] = compute_allSequences_NaNseparated_merge_stopping_periods(Position_Data,stopping_period_inds,boundaries,delxThr,deltThr,binSize,spikeSampRate)

if size(boundaries,1)==1
    boundaries_merged = boundaries;
    lengths_merged = boundaries_merged(:,2)-boundaries_merged(:,1);
    return
end

%indicies of end of last sequence, beginning of next
ind_pre = boundaries(1:end-1,2);
ind_post = boundaries(2:end,1);
x_in_between_pre_post = {};
for i=1:length(ind_pre)
    %spatial positions
        x_in_between_pre_post{i} = Position_Data(stopping_period_inds(ind_pre(i)):stopping_period_inds(ind_post(i)),2)';
end

    
%times
t_pre = Position_Data(stopping_period_inds(ind_pre),1);
t_post = Position_Data(stopping_period_inds(ind_post),1);

delx = nan(length(ind_pre),1);
%sum the amount of movement between re-entries into the reward zone
for i = 1:length(ind_pre)
      delx(i) =  sum(abs(diff(x_in_between_pre_post{i})));
end

%temporal jump between re-entries into the reward zone
    delt = t_post - t_pre;
    delt = delt/spikeSampRate;

%find jumps with delx<delxThr and delt<deltThr
    ind_merge = find(delx<delxThr/binSize & delt<deltThr);
    
    boundaries_merged = boundaries;
    boundaries_merged(ind_merge,2) = 0;
    boundaries_merged(ind_merge+1,1) = 0;
    boundaries_merged(sum(boundaries_merged,2)==0,:) = [];
    
    ind = find(boundaries_merged(:,2)==0);
    boundaries_merged(ind,2) = boundaries_merged(ind+1,2);
    ind = find(boundaries_merged(:,1)==0);
    boundaries_merged(ind,:) = [];    
    
%compute new lengths of sequences
    lengths_merged = boundaries_merged(:,2)-boundaries_merged(:,1);
    
%     plot(boundaries(:,1),'.-')
%     hold on, plot(boundaries(:,2),'.-'), hold off
%     keyboard

