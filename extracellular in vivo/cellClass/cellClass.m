function  [CellClass] = cellClass (waves, varargin)

% classifies clusters to PYR \ INT according to waveform parameters
% (trough-to-peak and spike width, also calculates asymmetry).
% 
% INPUT
%   waves       matrix of sampels (rows) x units (columns). for example:
%               waves = cat(1, spikes.rawWaveform{spikes.su})'
%               waves = cat(1, spikes.rawWaveform{:})'
%   man         logical. manual selection or automatic via known values
%   fs          sampling frequency
%   basepath    recording session path {pwd}
%   graphics    plot figure {1}.
%   saveFig     save figure {1}.
%   saveVar     save variable {1}.
% 
% OUTPUT
%   CellClass   struct with fields:
%       pyr         logical vector where 1 = PYR and 0 = INT
%       tp          trough-to-peak times
%       spkw        spike width
%       asym        asymmetry
% 
% DEPENDENCIES
%   getWavelet      from buzcode
%   fft_upsample    from Kamran Diba
%
% 08 apr 19 LH. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'man', false, @isnumeric);
addOptional(p, 'fs', 24414.14, @isnumeric);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'saveFig', false, @islogical);
addOptional(p, 'saveVar', false, @islogical);

parse(p, varargin{:})
basepath = p.Results.basepath;
man = p.Results.man;
fs = p.Results.fs;
graphics = p.Results.graphics;
saveFig = p.Results.saveFig;
saveVar = p.Results.saveVar;

% params
nunits = size(waves, 2);
upsamp = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trough-to-peak time [ms]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : nunits
    w = fft_upsample(waves(:, i), upsamp);
    [minval, minpos] = min(w);
    [maxval, maxpos] = max(w(1 : minpos - 1));   
    [maxvalpost, maxpost] = max(w(minpos : end));               
    if ~isempty(maxpost)
        % trough-to-peak - Stark et al., 2013; Bartho et al., 2004
        tp(i) = maxpost;
        if ~isempty(maxval)
            % asymmetry - Sirota et al., 2008
            asym(i) = (maxvalpost - maxval) / (maxvalpost + maxval);
        end
    else
        warning('waveform may be corrupted')
        tp(i) = NaN;
        asym(i) = NaN;
    end
end
% samples to ms
tp = tp / fs * 1000 / upsamp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spike width by inverse of max frequency in spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1 : nunits
    w = waves(:, i);
%     w = fft_upsample(waves(:, i), upsamp);
    w = [w(1) * ones(1000, 1); w; w(end) * ones(1000, 1)];
    [wave f t] = getWavelet(w, fs, 500, 3000, 128);
    % wt = cwt(w, fs, 'amor', 'FrequencyLimits', [500 3000]);
    wave = wave(:, int16(length(t) / 4) : 3 * int16(length(t) / 4));
    
    % find maximum
    [maxPow, ix] = max(wave);
    [~, mix] = max(maxPow);
    ix = ix(mix);
    spkw(i) = 1000 / f(ix);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if graphics
    h = figure;
    scatter(tp, spkw)
    xlabel('trough-to-peak [ms]')
    ylabel('spike width [ms]')
    xlim([0 max(tp) + 0.1])
    ylim([0 max(spkw) + 0.1])
    hold on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select PYRs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate separatrix according to known values 
xx = [0 0.8];
yy = [2.4 0.4];
m = diff(yy) / diff(xx);
b = yy(1) - m * xx(1);  % y = ax + b

if ~man     % automatic selection according to separatrix
    pyr = spkw >= m * tp + b;
else        % manual selection of boundary, with separatrix as a guide   
    fprintf('\nDiscriminate pyr and int (select Pyramidal)\n\n');
    [pyr, boundary] = selectCluster([tp; spkw] ,[m b], h);
end

CellClass.pyr = pyr;
CellClass.tp = tp;
CellClass.spkw = spkw;
CellClass.asym = asym;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveVar   
    [~, filename] = fileparts(basepath);
    save([basepath, '\', filename, '.cellClass.mat'], 'CellClass')
end

end

% EOF