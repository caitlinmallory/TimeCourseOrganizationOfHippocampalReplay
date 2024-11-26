function max_jump_distance = compute_max_jump_distance(posterior_com,bin_width,track_length) 

% ToDo: Figure out how to handle times when posterior was low.
if length(posterior_com) < 2
    max_jump_distance = nan;
    return
end
    
jump_distance = diff(posterior_com);
max_jump_distance = max(abs(jump_distance));

%convert from bins to cm:
max_jump_distance = max_jump_distance*bin_width;

%report as a fraction of the track length
max_jump_distance = max_jump_distance/track_length;