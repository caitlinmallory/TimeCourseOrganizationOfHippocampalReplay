function [LFP_filt_Phase,LFP_filt_Amp,LFP_filt] = compute_filteredLFP(freqRange,LFP,sampRate)

%bandpass filter
[FilterA,FilterB]=butter(2,freqRange/(sampRate/2));

%apply filter to LFP
LFP_filt = filtfilt(FilterA,FilterB,LFP);

%Hilbert Transform
LFP_filt_hilbert = hilbert(LFP_filt);
LFP_filt_Phase = wrapTo2Pi(angle(LFP_filt_hilbert));
LFP_filt_Amp = abs(LFP_filt_hilbert);

