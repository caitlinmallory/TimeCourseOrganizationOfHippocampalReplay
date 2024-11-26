function [x_edges, y_edges, x_centers, y_centers, x_edges_coarse, y_edges_coarse, x_centers_coarse, y_centers_coarse]...
    = compute_bin_edges(min_variable_x,max_variable_x,min_variable_y,max_variable_y,bin_width,bin_width_coarse)

nx_steps = round((max_variable_x - min_variable_x)/bin_width);
ny_steps = round((max_variable_y - min_variable_y)/bin_width);
x_edges = linspace(min_variable_x,max_variable_x,nx_steps+1);
y_edges = linspace(min_variable_y,max_variable_y,ny_steps+1);
x_centers = (x_edges(1:end-1) + x_edges(2:end))/2;
y_centers = (y_edges(1:end-1) + y_edges(2:end))/2;


nx_steps_coarse = round((max_variable_x - min_variable_x)/bin_width_coarse);
ny_steps_coarse = round((max_variable_y - min_variable_y)/bin_width_coarse);
x_edges_coarse = linspace(min_variable_x,max_variable_x,nx_steps_coarse+1);
y_edges_coarse = linspace(min_variable_y,max_variable_y,ny_steps_coarse+1);
x_centers_coarse = (x_edges_coarse(1:end-1) + x_edges_coarse(2:end))/2;
y_centers_coarse = (y_edges_coarse(1:end-1) + y_edges_coarse(2:end))/2;