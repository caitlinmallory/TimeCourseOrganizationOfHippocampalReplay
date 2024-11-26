all_rats = [1:11];

for rat_num = 1:length(all_rats)
    rats=all_rats(rat_num);
    combine_replay_by_stopping_period
    drawnow()
end