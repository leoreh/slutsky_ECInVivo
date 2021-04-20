function cc = cellclass(varargin)

% classifies clusters to PYR \ INT according to waveform parameters. note
% regarding upsampling: Peyrache uses resample followed by a gaussian
% filter. fft_upsample (by Stark) and interpft (Matlab) produce the same
% result with minimum edge effects.
% 
% INPUT
%   waves       matrix of sampels (rows) x units (columns). for example:
%               waves = cat(1, spikes.rawWaveform{spikes.su})'
%               waves = cat(1, spikes.rawWaveform{:})'
%   mfr         mean firing rate for each unit 
%   man         logical. manual selection or automatic via known values
%   fs          sampling frequency
%   basepath    recording session path {pwd}
%   graphics    plot figure {1}.
%   saveVar     save variable {1}.
% 
% OUTPUT
%   CellClass   struct with fields:
%       pyr         logical vector where 1 = PYR and 0 = INT
%       tp          trough-to-peak [ms] (bartho et al., 2004)
%       spkw        spike width [ms] (stark et al., 2013)
%       asym        asymmetry (Sirota et al., 2008)
%       hpk         half-peak width [ms] (Medrihan et al., 2017)
% 
% DEPENDENCIES
%   getWavelet      from buzcode
%   fft_upsample    from Kamran Diba
% 
% TO DO LIST
%   # add tail slope (Torrado Pacheco et al., Neuron, 2021). TP threshold
%   in that article was ~0.4 ms
%
% 08 apr 19 LH      updates: 
% 21 nov 19 LH      added spikes and FR in scatter
% 14 may 20 LH      added upsampling by fft

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'waves', []);
addOptional(p, 'mfr', [], @isnumeric);
addOptional(p, 'basepath', pwd);
addOptional(p, 'man', false, @islogical);
addOptional(p, 'fs', 20000, @isnumeric);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'saveVar', false, @islogical);

parse(p, varargin{:})
waves = p.Results.waves;
mfr = p.Results.mfr;
basepath = p.Results.basepath;
man = p.Results.man;
fs = p.Results.fs;
graphics = p.Results.graphics;
saveVar = p.Results.saveVar;

% params
if isempty(mfr)
    mfr = ones(1, size(waves, 2)) * 20;
else
    mfr = rescale(mfr, 10, 50);
end

nunits = size(waves, 2);

upsamp = 10;
nsamp = 10 * size(waves, 1);
nfs = fs * upsamp;

% for ALT 2 of spkw
fb = cwtfilterbank('SignalLength', nsamp, 'VoicesPerOctave', 32,...
    'SamplingFrequency', nfs, 'FrequencyLimits', [1 fs / 4]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc metrics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\ncalculating waveform parameters\n\n')
tic

% initialize
tp = nan(1, nunits);
spkw = nan(1, nunits);
hpk = nan(1, nunits);
asym = nan(1, nunits);

for i = 1 : nunits
    
    w = waves(:, i);
    w = interpft(w, nsamp); % upsample in the frequency domain
    
    % ---------------------------------------------------------------------
    % exclude positive or distorted spikes
    
    % since waveforms were detrended in getSPKfromDat, there mean is close
    % to 0 and thus the baseline does not need to be calculated. An
    % alternative is to calculate the mean of the waveform before the
    % spike. if spike detection was done with spktimesWh this should not
    % happen
    %     bline = 0;
    %     [~, abspk] = max(abs(w));
    %     if w(abspk) > bline || abspk > length(w) * 0.75
    %         fprintf('\nclu %d exhibits positive spike, skipping\n', i)
    %         tp(i) = NaN;
    %         spkw(i) = NaN;
    %         hpk(i) = NaN;
    %         asym(i) = NaN;
    %         continue
    %     end
  
    % ---------------------------------------------------------------------
    % trough-to-peak time [ms]
    
    [~, minpos] = min(w);
    [maxval, ~] = max(w(1 : minpos - 1));   
    [maxvalpost, maxpost] = max(w(minpos + 1 : end));               
    if ~isempty(maxpost)
        % trough-to-peak - Bartho et al., 2004
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
    
    % ---------------------------------------------------------------------
    % half peak width
    wu = w / maxvalpost;
    th1 = find(wu(minpos : maxpost + minpos) < 0.5);
    if any(th1)
        th1 = maxpost - th1(end);
        th2 = find(wu(maxpost + minpos : end) < 0.5);        
        if any(th2)
            th2 = th2(1);
        else
            th2 = th1;
        end
    else
        th1 = NaN;
        th2 = NaN;
    end
    hpk(i) = (th1 + th2) * 1000 / nfs;
    
    % ---------------------------------------------------------------------
    %  spike width by inverse of max frequency in spectrum

    % ALT 1: getWavelet
    %       boundary value replication (symmetrization)
    w = [w(1) * ones(1000, 1); w; w(end) * ones(1000, 1)];
    [cfs, f, ~] = getWavelet(w, nfs, 400, 3000, 128);
    [maxPow, ix] = max(cfs);
    [~, mix] = max(maxPow);
    ix = ix(mix);
    spkw(i) = 1000 / f(ix);
    
    % ALT 2: MATLAB cwt
%         [cfs, f, ~] = cwt(w, 'FilterBank', fb);
%         [~, ifreq] = max(abs(squeeze(cfs)), [], 1);
%         maxf = f(ifreq(round(length(w) / 2)));
%         spkwnew(i) = 1000 / maxf;
end

% samples to ms
tp = tp * 1000 / nfs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate separatrix according to known values 
xx = [0 0.8];
yy = [2.4 0.4];
m = diff(yy) / diff(xx);
b = yy(1) - m * xx(1);  % y = ax + b
sep = [m b];

if graphics  
    s = scatter(tp, spkw, mfr, 'filled');
    hold on
    xlabel('trough-to-peak [ms]')
    ylabel('spike width [ms]')
    xb = get(gca, 'XLim');
    yb = get(gca, 'YLim');
    plot(xb, [sep(1) * xb(1) + sep(2), sep(1) * xb(2) + sep(2)])
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

cc.pyr = pyr;
cc.tp = tp;
cc.spkw = spkw;
cc.asym = asym;
cc.hpk = hpk;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveVar   
    [~, basename] = fileparts(basepath);
    save([basepath, '\', basename, '.cellClass.mat'], 'cc')
    
    % cell metrics
    cmName = [basename, '.cell_metrics.cellinfo.mat'];
    if exist(cmName, 'file')
        load(cmName)
        cell_metrics.cc_spkw = cc.spkw;
        cell_metrics.cc_tp = cc.tp;
        cell_metrics.cc_asym = cc.asym;
        cell_metrics.cc_hpk = cc.hpk;
        save(cmName, 'cell_metrics')
    end
end

fprintf('\nthat took %.1f minutes\n', toc / 60)

return

% EOF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paths used to compare CE results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basepaths{1} = 'G:\RA\hDLX_Gq_WT2\200820_bslDay1';
basepaths{2} = 'G:\RA\hDLX_Gq_Tg\210820_bslDay2Raw2';
basepaths{3} = 'D:\Data\lh86\lh86_210301_072600';
basepaths{4} = 'G:\lh81\lh81_210207_045300';
cell_metrics = CellExplorer('basepaths', basepaths);

