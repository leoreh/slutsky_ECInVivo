function [bands, info] = calc_bands(varargin)

% calculates power in specific frequency bands from a psd
% 
% INPUT:
%   psdData         numeric. if mat than will calc bands for each row
%   freq            numeric. vec of frequency values, corresponding to
%                   psdData rows
%   flgNormBand     logical. normalize bands to broadband
% 
% OUTPUT
%   bands
%
% DEPENDENCIES
%   none
% 
% TO DO LIST
%
% 20 mar 24 LH      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'psdData', [], @isnumeric);
addOptional(p, 'freq', 1250, @isnumeric);
addOptional(p, 'flgNormBand', false, @islogical);

parse(p, varargin{:})
psdData         = p.Results.psdData;
freq            = p.Results.freq;
flgNormBand     = p.Results.flgNormBand;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ignore power-line interference freq
ignorePli = true; 
freq_pli = [47, 53];

% band params
bandNames = ["broad", "swa", "delta", "theta", "beta", "gamma"];
bandFreqs = [0.5, 100; 0.5, 1; 1, 4; 6, 12; 10, 30; 30, 80];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get indices to pli freq
if ignorePli
    lim_freq = (freq >= 0.5 & freq < freq_pli(1)) | (freq > freq_pli(2) & freq < Inf);
else
    lim_freq = (freq >= 0.5 & freq < Inf);
end

% normalize psd to broadband. running this before calculating bands is the
% same as dividing bands w/ broadband
if flgNormBand
    psdData = psdData ./ sum(psdData(:, lim_freq), 2);

    sbands = [2 : length(bandNames)];
    ytxt = 'Norm. PSD';
else
    sbands = [1 : length(bandNames)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calc power in band. bands is a 3d array of freqBand x state x session.
% psd is a 3d array of state x freq x session.
for iband = 1 : length(bandFreqs)
    
    bandIdx = InIntervals(freq, bandFreqs(iband, :));

    if InIntervals(freq_pli, bandFreqs(iband, :))
        bandIdx = bandIdx & ~InIntervals(freq, freq_pli);
    end
    
    bands(iband, :, :) = sum(psdData(:, bandIdx), 2);
end

% organize info
info.ignorePli = ignorePli;
info.freq_pli = freq_pli;
info.bandNames = bandNames;
info.bandFreqs = bandFreqs;
info.flgNormBand = flgNormBand;

end

% EOF