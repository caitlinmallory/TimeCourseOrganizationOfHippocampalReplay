function [boundaries_data,lengths_data] = compute_allSequences_NaNseparated_cm(x)

%function separates long 
% trajectory interrupted by NaNs into smaller, continuous trajectories

%replace data with 1's, NaN's with 0's
    y = ~isnan(x(:,1));

%function seqle finds chains composed of real data (separated by NaNs)
    z = seqle(y);
    identities = z{1}; %1 for data chain, 0 for NaN chain
    lengths = z{2};%length of each chain
    boundaries = z{3}; %chain boundaries

%Caitlin added the following: if you had a big jump followed by a small
%jump (i.e., nan followed by a one here), you want to keep the time bin following the big jump as part of the
%next sequence.
    


    for i = 2:length(z{1})
        if identities(i-1) == 0
            boundaries(1,i) = boundaries(1,i)-1;
            lengths(i) = lengths(i)+1;
        end
    end

%extract sequence
    lengths_data = lengths(identities==1)';
    boundaries_data = boundaries(:,identities==1)';    
%remove sequences that are 1 unit long
    %ind_remove = find(lengths_data==0); %I think this was a typo- changed
    %to 1. I was wrong. changed it back to 0.
    ind_remove = find(lengths_data==0);
    lengths_data(ind_remove,:) = [];
    boundaries_data(ind_remove,:) = [];
        
    
    
    
