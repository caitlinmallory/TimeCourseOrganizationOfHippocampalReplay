function [boundaries_merged,lengths_merged] = compute_allSequences_NaNseparated_merge(x,boundaries,delxThr,deltThr,binSize)



if size(boundaries,1)==1
    boundaries_merged = boundaries;
    lengths_merged = boundaries_merged(:,2)-boundaries_merged(:,1);
    return
end

%indicies of end of last sequence, beginning of next
ind_pre = boundaries(1:end-1,2);
ind_post = boundaries(2:end,1);

%spatial positions
if size(x,2)>2
    x_pre = x(ind_pre,2:3);
    x_post = x(ind_post,2:3);
else
    x_pre = x(ind_pre,2);
    x_post = x(ind_post,2);
end

%times
t_pre = x(ind_pre,1);
t_post = x(ind_post,1);

%spatial jump between sequences
if size(x,2)>2
    delx = sqrt((x_post(:,1)-x_pre(:,1)).^2 + (x_post(:,2)-x_pre(:,2)).^2);
else
    delx = sqrt((x_post(:,1)-x_pre(:,1)).^2);
end

%temporal jump between sequences
delt = t_post - t_pre;

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

