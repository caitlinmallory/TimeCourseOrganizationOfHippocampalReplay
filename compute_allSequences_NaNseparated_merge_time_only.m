function [boundaries_merged,lengths_merged] = compute_allSequences_NaNseparated_merge_time_only(x,boundaries,deltThr)



if size(boundaries,1)==1
    boundaries_merged = boundaries;
    lengths_merged = boundaries_merged(:,2)-boundaries_merged(:,1);
    return
end

%indicies of end of last sequence, beginning of next
ind_pre = boundaries(1:end-1,2);
ind_post = boundaries(2:end,1);

%times
t_pre = x(ind_pre,1);
t_post = x(ind_post,1);

%temporal jump between sequences
delt = t_post - t_pre;

%find jumps with delx<delxThr and delt<deltThr
ind_merge = find(delt<deltThr);

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

