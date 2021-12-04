function st = spktimesMetrics(varargin)

%  calculate metrics of st_metrics from ACG. based on calc_ACG_metrics from
%  cell explorer. metrics include: busrtIndex_Royer2012, Mizuseki2012,
%  Doublets. calculates two ACGs: narrow (100ms, 0.5ms bins) and wide (1s,
%  1ms bins)


% INPUT:
%   spikes          struct (see getSpikes)
%   sunits          numeric vec. indices of selected units for calculation
%                   {[]}.
%   winCalc         cell array where each cell is an n x 2 mat of indices.
%                   metrices will be calculated for each cell by limiting
%                   spktimes to the within the times specified by the n x 2
%                   mat. can be for example ss.stateEpochs. must
%                   be the same units as spikes.times (e.g. [s])
%   basepath        path to recording
%   graphics        logical. plot graphics {true} or not (false)
%   saveVar         logical. save variables (update spikes and save su)
%
% OUTPUT:
%   
%
% DEPENDENCIES:
%   CCG
%
% TO DO LIST:
%
% 24 nov 21 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'spikes', []);
addOptional(p, 'sunits', []);
addOptional(p, 'winCalc', {[0 Inf]}, @iscell);
addOptional(p, 'basepath', pwd, @ischar);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'saveVar', true, @islogical);

parse(p, varargin{:})
spikes      = p.Results.spikes;
sunits      = p.Results.sunits;
winCalc     = p.Results.winCalc;
basepath    = p.Results.basepath;
graphics    = p.Results.graphics;
saveVar     = p.Results.saveVar;

% load spikes if empty
if isempty(spikes)
    [~, filename] = fileparts(basepath);
    spkname = [filename '.spikes.cellinfo.mat'];
    if exist(spkname, 'file')
        load(spkname)
    else
        error('%s not found', spkname)
    end
end

% make sure spikes has required fields
if ~all(isfield(spikes, {'shankID', 'cluID', 'times'}))
    error('spikes missing required fields')
end

% selected untis
if isempty(sunits)
    sunits = 1 : length(spikes.times);
end

% load session info
[~, basename] = fileparts(basepath);
sessionName = [basename, '.session.mat'];
if ~exist(sessionName, 'file')
    session = CE_sessionTemplate(pwd, 'viaGUI', false,...
        'force', true, 'saveVar', true);
else
    load(sessionName)
end

% acg params
st.info.runtime = datetime(now, 'ConvertFrom', 'datenum');
st.info.winCalc = winCalc;
st.info.acg_wide_bins = 500;
st.info.acg_wide_bnsz = 0.001;
st.info.acg_wide_dur = 1;
st.info.acg_narrow_bins = 100;
st.info.acg_narrow_bnsz = 0.0005;
st.info.acg_narrow_dur = 0.1;

% spk params
fs = session.extracellular.sr;
nunits = length(sunits);
nwin = length(winCalc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
st.acg_wide     = zeros(st.info.acg_wide_dur / st.info.acg_wide_bnsz + 1, nunits, nwin);
st.acg_narrow   = zeros(st.info.acg_narrow_dur / st.info.acg_narrow_bnsz + 1, nunits, nwin);
st.doublets     = zeros(nunits, nwin);
st.royer        = zeros(nunits, nwin);
st.lidor        = zeros(nunits, nwin);
st.mizuseki     = zeros(nunits, nwin);
st.cv           = zeros(nunits, nwin);
st.cv2          = zeros(nunits, nwin);
st.lv           = zeros(nunits, nwin);
st.lvr          = zeros(nunits, nwin);

for iunit = sunits
    for iwin = 1 : length(winCalc)
            
        % limit spktimes to window
        spkIdx = InIntervals(spikes.times{iunit}, winCalc{iwin});
        spktimes = spikes.times{iunit}(spkIdx);
        nspks = length(spktimes);
        isi = diff(spktimes);
        nisi = length(isi);
        
        % acg
        [st.acg_wide(:, iunit, iwin), st.info.acg_wide_tstamps] = CCG(spktimes,...
            ones(size(spktimes)), 'binSize', st.info.acg_wide_bnsz,...
            'duration', st.info.acg_wide_dur, 'norm', 'rate', 'Fs', 1 / fs);
        
        [st.acg_narrow(:, iunit, iwin), st.info.acg_narrow_tstamps] = CCG(spktimes,...
            ones(size(spktimes)), 'binSize', st.info.acg_narrow_bnsz,...
            'duration', st.info.acg_narrow_dur, 'norm', 'rate', 'Fs', 1 / fs);
        
        % burstiness ------------------------------------------------------
        
        % doublets: max bin count from 2.5-8ms normalized by the average number
        % of spikes in the 8-11.5ms bins
        st.doublets(iunit, iwin) = max(st.acg_narrow(st.info.acg_narrow_bins + 1 + 5 : st.info.acg_narrow_bins + 1 + 16, iunit, iwin)) /...
            mean(st.acg_narrow(st.info.acg_narrow_bins + 1 + 16 : st.info.acg_narrow_bins + 1 + 23, iunit, iwin));
        
        % royer 2012: average number of spikes in the 3-5 ms bins divided by the
        % average number of spikes in the 200-300 ms bins
        st.royer(iunit, iwin) = mean(st.acg_wide(st.info.acg_wide_bins + 1 + 3 : st.info.acg_wide_bins + 1 + 5, iunit, iwin)) /...
            mean(st.acg_wide(st.info.acg_wide_bins + 1 + 200 : st.info.acg_wide_bins + 1 + 300, iunit, iwin));
        
        % lidor: sum of spikes in 2-10 ms normalized to sum in 35-50
        t1 = find(st.info.acg_narrow_tstamps > 0.002 & st.info.acg_narrow_tstamps < 0.01);
        t2 = find(st.info.acg_narrow_tstamps > 0.035 & st.info.acg_narrow_tstamps < 0.05);
        burst_temp = sum(st.acg_narrow(t1, iunit, iwin));
        bl_temp = sum(st.acg_narrow(t2, iunit, iwin));
        st.lidor(iunit, iwin) = (burst_temp - bl_temp) ./ (burst_temp + bl_temp);
        
        % Mizuseki 2011: fraction of spikes with a ISI for following or preceding
        % spikes < 0.006
        burst_temp = zeros(1, length(spktimes) - 1);
        for ispk = 2 : length(spktimes) - 1
            burst_temp(ispk) = any(diff(spktimes(ispk - 1 : ispk + 1)) < 0.006);
        end
        st.mizuseki(iunit, iwin) = sum(burst_temp > 0) / length(burst_temp);
        
        % firing irregularity ---------------------------------------------
        
        % Cv (coefficient of variation): Shinomoto 2003
        st.cv(iunit, iwin) = std(isi) / mean(isi);
        
        % Cv2 (local cv): Holt 1996, taken from CE
        cv2_temp = 2 * abs(isi(1 : end - 1) - isi(2 : end)) ./...
            (isi(1 : end - 1) + isi(2 : end));
        st.cv2(iunit, iwin) = mean(cv2_temp(cv2_temp < 1.95));
        
        % Lv: Shinomoto 2003 and Kobayashi 2019.
        lv_term = 0;
        for ispk = 1 : nisi - 1
            diff_term = 3 * (isi(ispk) - isi(ispk + 1))^2;
            sum_term = (isi(ispk) + isi(ispk + 1))^2;
            lv_term = lv_term + diff_term / sum_term;
        end
        st.lv(iunit, iwin) = lv_term / (nisi - 1);
        
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
        st.lvr(iunit, iwin) = 3 / (nisi - 1) * lv_term;
        
        % Fano factor: variability in fr relative to mfr
        % brustiness.ff = var(fr) / mean(fr)
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% if graphics
%     plot_spktimesMetrics()
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if saveVar
       
    [~, basename] = fileparts(basepath);
    save([basepath, filesep, basename, '.st_metrics.mat'], 'st')

    % update cell metrics
    cmName = [basename, '.cell_metrics.cellinfo.mat'];
    if exist(cmName, 'file')
        load(cmName)
        cell_metrics.st_doublets    = st.doublets;
        cell_metrics.st_royer       = st.royer;
        cell_metrics.st_mizuseki    = st.mizuseki;
        cell_metrics.st_lv          = st.lv;
        cell_metrics.st_cv          = st.cv;
        cell_metrics.st_cv2         = st.cv2;
        save(cmName, 'cell_metrics')
    end
end

end

% EOF