function [speedTuningCurve, speed_bins] = compute_speedTuning(Run_Times,spkTime,spkSpeed,positions,posSampRate,speedThr)

speed_bins = linspace(speedThr,30,20)';

spkSpeed_filt = compute_dataTemporalConcatenation([spkTime,spkSpeed],Run_Times);

spikeDensity = histc(spkSpeed_filt(:,2),speed_bins); spikeDensity = spikeDensity(1:end-1);
posDensity = histc(positions(:,5),speed_bins); posDensity = posDensity(1:end-1);

speedTuningCurve = spikeDensity./posDensity/(1/posSampRate);
        