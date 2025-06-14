function st = spktimes_metrics(varargin)

%  calculate spike timing metrics from ACG and ISI histogram. based in part
%  on calc_ACG_metrics from cell explorer. metrics include:
%  busrtIndex_Royer2012, Mizuseki2012, Doublets. calculates two ACGs:
%  narrow (100ms, 0.5ms bins) and wide (1s, 1ms bins)


% INPUT:
%   spktimes        cell of spike times in s. 
%   fs              numeric. sampling frequency. 
%   sunits          numeric vec. indices of selected units for calculation
%                   {[]}.
%   bins            cell array of n x 2 mats of intervals.
%                   metrices will be calculated for each cell by limiting
%                   spktimes to the intervals. can be for example
%                   ss.boutTimes. must be the same units as spikes.times
%                   (e.g. [s])
%   basepath        path to recording
%   flg_save         logical. save variables (update spikes and save su)
%   flg_force       logical. force analysis even if struct file exists
%                   {false}
%   flg_all         logical. analyse all parameters (slow) {true}
%
% OUTPUT:
%   st              struct
%
% DEPENDENCIES:
%   CCG
%
% TO DO LIST:
%   rmv dependency on spikes struct and cell explorer
%   add metric from Miura et al., j. neurosci., 2007
%
% 24 nov 21 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'spktimes', []);
addOptional(p, 'fs', [], @isnumeric);
addOptional(p, 'sunits', []);
addOptional(p, 'bins', {[0 Inf]});
addOptional(p, 'basepath', pwd, @ischar);
addOptional(p, 'flg_save', true, @islogical);
addOptional(p, 'flg_force', false, @islogical);
addOptional(p, 'flg_all', true, @islogical);

parse(p, varargin{:})
spktimes    = p.Results.spktimes;
fs          = p.Results.fs;
sunits      = p.Results.sunits;
bins        = p.Results.bins;
basepath    = p.Results.basepath;
flg_save     = p.Results.flg_save;
flg_force   = p.Results.flg_force;
flg_all     = p.Results.flg_all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% file names
[~, basename] = fileparts(basepath);
stFile = fullfile(basepath, [basename, '.st_metrics.mat']);
spkFile = [basename '.spikes.cellinfo.mat'];
sessionFile = [basename, '.session.mat'];

% check if already analyzed 
if exist(stFile, 'file') && ~flg_force
    load(stFile)
    return
end

% load spikes if empty
if isempty(spktimes)
    if exist(spkFile, 'file')
        load(spkFile)
    else
        error('spikes file does not exist, input spktimes')
    end
    spktimes = spikes.times;
end

if isempty(bins)
    bins = [0 Inf];
end
if ~iscell(bins)
    bins = {bins};
end

% selected untis
if isempty(sunits)
    sunits = 1 : length(spktimes);
end

% load session info for fs
if exist(sessionFile, 'file')
    load(sessionFile)
end
if isempty(fs)
    if exist('session', 'var')
        fs = session.extracellular.sr;
    else
        fs = 20000;
    end
end

% acg params
st.info.runtime = datetime("now");
st.info.bins = bins;
st.info.acg_wide_bins = 500;
st.info.acg_wide_bnsz = 0.001;
st.info.acg_wide_dur = 1;
st.info.acg_narrow_bins = 100;
st.info.acg_narrow_bnsz = 0.0005;
st.info.acg_narrow_dur = 0.1;

% spk params
nunits = length(sunits);
nbins = length(bins);
minSpkThr = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
st.acg_wide     = nan(st.info.acg_wide_dur / st.info.acg_wide_bnsz + 1, nbins, nunits);
st.acg_narrow   = nan(st.info.acg_narrow_dur / st.info.acg_narrow_bnsz + 1, nbins, nunits);
st.doublets     = nan(nbins, nunits);
st.royer        = nan(nbins, nunits);
st.royer2       = nan(nbins, nunits);
st.lidor        = nan(nbins, nunits);
st.mizuseki     = nan(nbins, nunits);
st.cv           = nan(nbins, nunits);
st.cv2          = nan(nbins, nunits);
st.lv           = nan(nbins, nunits);
st.lvr          = nan(nbins, nunits);

for iunit = 1 : length(sunits)
    for ibin = 1 : length(bins)
            
        % limit spktimes to window
        spkIdx = InIntervals(spktimes{sunits(iunit)}, bins{ibin});
        st_unit = spktimes{sunits(iunit)}(spkIdx);
        nspks = length(st_unit);
        isi = diff(st_unit);
        nisi = length(isi);
        
        if nspks < minSpkThr
            continue
        end
        
        % acg
        [st.acg_wide(:, ibin, iunit), st.info.acg_wide_tstamps] = CCG(st_unit,...
            ones(size(st_unit)), 'binSize', st.info.acg_wide_bnsz,...
            'duration', st.info.acg_wide_dur, 'norm', 'rate', 'Fs', 1 / fs);
        
        [st.acg_narrow(:, ibin, iunit), st.info.acg_narrow_tstamps] = CCG(st_unit,...
            ones(size(st_unit)), 'binSize', st.info.acg_narrow_bnsz,...
            'duration', st.info.acg_narrow_dur, 'norm', 'rate', 'Fs', 1 / fs);
        
        % burstiness ------------------------------------------------------
        
        % doublets: max bin count from 2.5-8ms normalized by the average number
        % of spikes in the 8-11.5ms bins
        st.doublets(ibin, iunit) = max(st.acg_narrow(st.info.acg_narrow_bins + 1 + 5 : st.info.acg_narrow_bins + 1 + 16, ibin, iunit)) /...
            mean(st.acg_narrow(st.info.acg_narrow_bins + 1 + 16 : st.info.acg_narrow_bins + 1 + 23, ibin, iunit));
        
        % royer 2012: average number of spikes in the 3-5 ms bins divided by the
        % average number of spikes in the 200-300 ms bins
        st.royer(ibin, iunit) = mean(st.acg_wide(st.info.acg_wide_bins + 1 + 3 : st.info.acg_wide_bins + 1 + 5, ibin, iunit)) /...
            mean(st.acg_wide(st.info.acg_wide_bins + 1 + 200 : st.info.acg_wide_bins + 1 + 300, ibin, iunit));
        
        % royer 2012: max(0:10) / mean(40:50) normalized
        peakbrst = max(st.acg_wide(st.info.acg_wide_bins + 1 :...
            st.info.acg_wide_bins + 1 + 10, ibin, iunit));
        basebrst = mean(st.acg_wide(st.info.acg_wide_bins + 1 + 40 :...
            st.info.acg_wide_bins + 1 + 50, ibin, iunit));
        if peakbrst > basebrst
            st.royer2(ibin, iunit) = (peakbrst - basebrst) / peakbrst;
        else
            st.royer2(ibin, iunit) = (peakbrst - basebrst) / basebrst;
        end        

        % lidor: sum of spikes in 2-10 ms normalized to sum in 35-50
        t1 = find(st.info.acg_narrow_tstamps > 0.002 & st.info.acg_narrow_tstamps < 0.01);
        t2 = find(st.info.acg_narrow_tstamps > 0.035 & st.info.acg_narrow_tstamps < 0.05);
        burst_temp = sum(st.acg_narrow(t1, ibin, iunit));
        bl_temp = sum(st.acg_narrow(t2, ibin, iunit));
        st.lidor(ibin, iunit) = (burst_temp - bl_temp) ./ (burst_temp + bl_temp);
        
        % Mizuseki 2011: fraction of spikes with a ISI for following or preceding
        % spikes < 0.006
        burst_temp = zeros(1, length(st_unit) - 1);
        for ispk = 2 : length(st_unit) - 1
            burst_temp(ispk) = any(diff(st_unit(ispk - 1 : ispk + 1)) < 0.006);
        end
        st.mizuseki(ibin, iunit) = sum(burst_temp > 0) / length(burst_temp);
        
        % firing irregularity ---------------------------------------------
        
        % Cv (coefficient of variation): Shinomoto 2003. note that when
        % calculated for boutTimes, this is terribly biased due to large
        % ISI's between bouts
        st.cv(ibin, iunit) = std(isi) / mean(isi);
        
        % Cv2 (local cv): Holt 1996, taken from CE
        cv2_temp = 2 * abs(isi(1 : end - 1) - isi(2 : end)) ./...
            (isi(1 : end - 1) + isi(2 : end));
        st.cv2(ibin, iunit) = mean(cv2_temp(cv2_temp < 1.95));
        
        % Lv: Shinomoto 2003 and Kobayashi 2019.
        lv_term = 0;
        for ispk = 1 : nisi - 1
            diff_term = 3 * (isi(ispk) - isi(ispk + 1))^2;
            sum_term = (isi(ispk) + isi(ispk + 1))^2;
            lv_term = lv_term + diff_term / sum_term;
        end
        st.lv(ibin, iunit) = lv_term / (nisi - 1);
        
        % LvR: Shinomoto 2009
        ref = 0.005;                    % refractory constant [s]
        lv_term = 0;
        for ispk = 1 : nisi - 1
            ccorr_term = 4 * isi(ispk) * isi(ispk + 1);
            sum_term = isi(ispk) + isi(ispk + 1);
            left_term = 1 - ccorr_term / (sum_term ^ 2);
            right_term = 1 + (4 * ref) / sum_term;
            lv_term = lv_term + left_term * right_term;
        end
        st.lvr(ibin, iunit) = 3 / (nisi - 1) * lv_term;
        
        if flg_all
            % fit triple exponential to acg. adapted from CE (fit_ACG.m).
            % requires the Curve Fitting Toolbox. no idea whats going on here
            g = fittype('max(c*(exp(-(x-f)/a)-d*exp(-(x-f)/b))+h*exp(-(x-f)/g)+e,0)',...
                'dependent',{'y'},'independent',{'x'},...
                'coefficients',{'a','b','c','d','e','f','g','h'});
            a0 = [20, 1, 30, 2, 0.5, 5, 1.5, 2];
            lb = [1, 0.1, 0, 0, -30, 0, 0.1, 0];
            ub = [500, 50, 500, 15, 50, 20, 5, 100];
            offset = 101;
            x = ([1 : 100] / 2)';
            [f0, ~] = fit(x, st.acg_narrow(x * 2 + offset, ibin, iunit),...
                g, 'StartPoint', a0, 'Lower', lb, 'Upper', ub);
            fit_params = coeffvalues(f0);
            st.tau_rise(ibin, iunit) = fit_params(2);
            st.tau_burst(ibin, iunit) = fit_params(7);
            st.acg_refrac(ibin, iunit) = fit_params(6);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% plot_spktimesMetrics()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flg_save
    
    save(stFile, 'st')

    % update cell metrics
    cmName = [basename, '.cell_metrics.cellinfo.mat'];
    if exist(cmName, 'file')
        load(cmName)
        cell_metrics.st_doublets    = st.doublets;
        cell_metrics.st_royer       = st.royer;
        cell_metrics.st_royer2      = st.royer2;
        cell_metrics.st_lidor       = st.lidor;
        cell_metrics.st_mizuseki    = st.mizuseki;
        cell_metrics.st_cv          = st.cv;
        cell_metrics.st_cv2         = st.cv2;
        cell_metrics.st_lv          = st.lv;
        cell_metrics.st_lvr         = st.lvr;
        save(cmName, 'cell_metrics')
    end
end

end

% EOF