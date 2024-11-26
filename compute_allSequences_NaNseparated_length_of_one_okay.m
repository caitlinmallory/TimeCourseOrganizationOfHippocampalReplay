function [boundaries_data,lengths_data] = compute_allSequences_NaNseparated_length_of_one_okay(x)

%function separates long trajectory interrupted by NaNs into smaller, continuous trajectories

%replace data with 1's, NaN's with 0's
    y = ~isnan(x(:,1));

%function seqle finds chains composed of real data (separated by NaNs)
    z = seqle(y);
    identities = z{1}; %1 for data chain, 0 for NaN chain
    lengths = z{2};%length of each chain
    boundaries = z{3}; %chain boundaries
    
%extract sequence
    lengths_data = lengths(identities==1)';
    boundaries_data = boundaries(:,identities==1)';
    
%remove sequences that are 1 unit long
    %ind_remove = find(lengths_data==0); %I think this was a typo- changed
    %to 1.
%     ind_remove = find(lengths_data==1);
%     lengths_data(ind_remove,:) = [];
%     boundaries_data(ind_remove,:) = [];
        
    
    
    
