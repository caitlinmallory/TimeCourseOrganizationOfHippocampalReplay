load Position_Data
load Experiment_Information

num_segs = length(Experiment_Information.Segments);
for i = 1:num_segs
    Position_Data_sub = compute_dataTemporalConcatenation(Position_Data,Experiment_Information.Segments(i).Times);
    figure()
    plot(Position_Data_sub(:,2),Position_Data_sub(:,3))
end