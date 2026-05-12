% Boaz pipeline — what it does
% 
% Four MATLAB functions, all interactive:
% 
% Ca_EventDetector.m is single-trace event detection. Input is raw
% fluorescence, not dF/F. It computes dF/F internally as F / prctile ( F ,
% p ) F/prctile(F,p) with p = 20 p=20 for cyto and p = 5 p=5 for mito,
% applies msbackadj (a spline baseline from the bioinformatics toolbox)
% plus optional Hamming / LOESS / Savitzky-Golay smoothing, then detects
% with findpeaks (height, prominence, distance). It includes a
% lag-correction step that looks back in the unfiltered dF/F to undo
% filter-induced peak shifts, and finally lets the user accept or reject
% each event one-by-one in a figure.
% 
% Ca_EventCompare.m pairs events. For each mito peak, it finds the closest
% preceding cyto peak (last cyto whose time is below the mito time), opens
% a figure with the segment around the mito peak, and asks the user to
% approve or reject the pair. The output is a two-column matrix of approved
% (cyto idx, mito idx).
% 
% preComp.m is plumbing — it concatenates dFF and event tables into the
% structures Ca_EventCompare expects. The Leore-folder utilities (spks2ca,
% manAlign) pair spike times to calcium peaks with a frame tolerance and an
% interactive slider-based time-alignment tool.