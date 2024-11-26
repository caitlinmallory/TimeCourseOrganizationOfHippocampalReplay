clearvars -except dayFiles day directory rat windows hand_clustered_only
%returns spikeDensity, measured as number of spikes in window with size
%equal to spikeDensityStepSize

%spikeDensity:
%1  2
%t  spikeDensity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restrict_to_excitatory_cells = 1;
restrict_to_well_isolated_cells = 0;

load Experiment_Information
load Analysis_Information
load clusters

spikeSampRate = Experiment_Information.spikeSampRate;
Run_Times = Experiment_Information.Run_Times;
Sleep_Times = Experiment_Information.Sleep_Times;

times_list = cat(2, (cat(2, Run_Times{:})), Sleep_Times{:});
Times_day = [min(times_list) max(times_list)];

%times to estimate spike density
timeEdges = (Times_day(1):spikeDensityStepSize*spikeSampRate:Times_day(2))';

%load clusters
if restrict_to_excitatory_cells == 1
    meets_excitatory_inhibitory_criterion = [clusters.Excitatory]';
else
    meets_excitatory_inhibitory_criterion = ones(height(clusters),1);
end

if restrict_to_well_isolated_cells == 1
    meets_well_isolated_criterion = [clusters.well_isolated]';
else
    meets_well_isolated_criterion = ones(height(clusters),1);
end

ind_cluster = find(meets_excitatory_inhibitory_criterion & meets_well_isolated_criterion);


%loop through clusters
spikeDensity = 0;
for i = 1:length(ind_cluster)
    spkTimes = clusters(ind_cluster(i)).spkTime;
    spkTimes(find(diff(spkTimes)==0),:) = [];
    if ~isempty(spkTimes)
        spikeDensity = spikeDensity + histc(spkTimes,timeEdges);
    end

    if length(find(diff(spkTimes)==0))>0
        keyboard
    end
end
spikeDensity = [timeEdges spikeDensity/spikeDensityStepSize];

%ToDO: why is this here?
ind = find(spikeDensity(:,2)>20000);

spikeDensity(ind,2) = spikeDensity(ind-1,2);


save('spikeDensity_pyr.mat','spikeDensity')

