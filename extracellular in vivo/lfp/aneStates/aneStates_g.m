function as = aneStates_g(varargin)

% analyze anesthesia states of one group (e.g. WT)
%
% INPUT
%   ch          vector. channel/s to analyze.
%   basepath    string. recording session path {pwd}
%   basename    string. mouse basename. if empty extracted from basepath
%   graphics    logical. plot figure {1}
%   saveVar     logical. save variable {1}
%   saveFig     logical. save figure {1}
%   forceA      logical {0}. force analysis even if .mat exists
%
% OUTPUT
%   as          struct
%
% TO DO LIST
%       #
%
% CALLS
%       aneStates.m
%       stdshade
%
% EXAMPLE
%
% 21 jan 20 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'basepath', pwd, @isstr);
addParameter(p, 'rm', [], @isnumeric);
addParameter(p, 'graphics', true, @islogical)
addParameter(p, 'saveVar', true, @islogical);
addParameter(p, 'saveFig', true, @islogical);
addParameter(p, 'forceA', false, @islogical);

parse(p, varargin{:})
basepath = p.Results.basepath;
rm = p.Results.rm;
graphics = p.Results.graphics;
saveVar = p.Results.saveVar;
saveFig = p.Results.saveFig;
forceA = p.Results.forceA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = 1250;
[~, grpname] = fileparts(basepath);
cd(basepath)
filename = dir('*.abf');
files = natsort({filename.name});
nfiles = 1 : length(files);
nfiles(rm) = [];
% if exist([grpname '_as.mat']) && ~forceA
%     load([grpname '_as.mat'])
%     return
% end
for i = 1 : length(nfiles)
    [~, basename] = fileparts(files{nfiles(i)});
    % load individual data  
    [bs, iis, ep] = aneStates_m('ch', 1, 'basepath', basepath,...
        'basename', basename, 'graphics', graphics, 'saveVar', saveVar,...
        'saveFig', saveFig, 'forceA', false, 'binsize', 30, 'smf', 7);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % arrange population data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    as.grp = grpname;
    as.mouse{i} = basename;
    as.recDur(i) = ep.recDur;
    as.deepDur(i) = ep.deepDur;
    as.surDur(i) = ep.surDur;
    as.iis{i} = iis.rate;
    as.bsr{i} = bs.bsr;
    as.dband{i} = ep.dband;
    as.t{i} = bs.cents / fs / 60;
    as.nspks(i) = ep.nspks;
    as.deep_nspks(i) = ep.deep_nspks;
    as.sur_nspks(i) = ep.sur_nspks;
    as.thr(i) = iis.thr(2);
    as.sur_delta(i) = ep.sur_delta;
    as.deep_delta(i) = ep.deep_delta;
    
    % after last recording is loaded, arrange in mat and save
    if nfiles(i) == nfiles(end)
        
        [~, idx] = max(as.recDur);      % index to longest recording
        maxdur = length(as.bsr{idx});
        
        % bsr, iis and delta
        mat = cellfun(@(x)[x(:); NaN(maxdur-length(x), 1)], as.bsr,...
            'UniformOutput', false);
        as.bsr = cell2mat(mat)';
        mat = cellfun(@(x)[x(:); NaN(maxdur-length(x), 1)], as.iis,...
            'UniformOutput', false);
        as.iis = cell2mat(mat)';
        mat = cellfun(@(x)[x(:); NaN(maxdur-length(x), 1)], as.dband,...
            'UniformOutput', false);
        as.dband = cell2mat(mat)';
        % timestamps
        mat = cellfun(@(x)[x(:); NaN(maxdur-length(x), 1)], as.t,...
            'UniformOutput', false);
        as.t = cell2mat(mat)';
        as.t = as.t(idx, :);
        
        if saveVar
            save([grpname '_as.mat'], 'as')
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
graphics = false;
if graphics
    
    % data for correlation
    cdata(:, 1) = as.iis(~isnan(as.iis));
    cdata(:, 2) = as.bsr(~isnan(as.bsr));
    cdata(:, 3) = as.dband(~isnan(as.dband));
    % for histograms     
    nbins = floor(length(as.mouse) / 2);
    % for xlim and time correlation
    xt = prctile(as.recDur, 50) / fs / 60;
    [~, idx] = min(abs(as.t - xt * fs * 60));
    
    figure;
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    suptitle(grpname)
    
    % iis
    sp1 = subplot(3, 4, [1 : 2]);
    plot(as.t / fs / 60, as.iis, 'LineWidth', 1)
    hold on
    stdshade(as.iis, 0.5, 'k', as.t / fs / 60)
    axis tight
    x = xlim;
    xlim([x(1) xt])
    ylabel('IIS rate [spikes / bin]')
    set(gca, 'TickLength', [0 0])
    box off
    title('IIS rate')
    [r, pval] = corr(mean(as.iis(:, 1 : idx), 'omitnan')',...
        as.t(1 : idx)', 'Type', 'Spearman');
    y = ylim;
    text(x(1), y(2), sprintf('r = %.2f, p = %.3f', r, pval))
    
    % bs
    sp2 = subplot(3, 4, [5 : 6]);
    plot(as.t / fs / 60, as.bsr, 'LineWidth', 1)
    hold on
    stdshade(as.bsr, 0.5, 'k', as.t / fs / 60)
    axis tight
    x = xlim;
    ylim([0 1])
    xlim([x(1) xt])
    ylabel('BSR')
    set(gca, 'TickLength', [0 0])
    box off
    title('BSR')
    [r, pval] = corr(mean(as.bsr(:, 1 : idx), 'omitnan')',...
        as.t(1 : idx)', 'Type', 'Spearman');
    y = ylim;
    text(x(1), y(2), sprintf('r = %.2f, p = %.3f', r, pval))
    
    % delta
    sp3 = subplot(3, 4, [9 : 10]);
    plot(as.t / fs / 60, as.dband, 'LineWidth', 1);
    hold on
    stdshade(as.dband, 0.5, 'k', as.t / fs / 60)
    axis tight
    x = xlim;
%     ylim([0 1])
    xlim([x(1) xt])
    ylabel('norm. delta power')
    xlabel('Time [m]')
    set(gca, 'TickLength', [0 0])
    box off
    title('Delta')
    [r, pval] = corr(mean(as.dband(:, 1 : idx), 'omitnan')',...
        as.t(1 : idx)', 'Type', 'Spearman');
    y = ylim;
    text(x(1), y(2), sprintf('r = %.2f, p = %.3f', r, pval))
    
    % iis vs bsr
    subplot(3, 4, 3)
    hold on
    for i = 1 : length(nfiles)
        scatter(as.bsr(i, :), as.iis(i, :), 50, '.')
    end
    xlabel('BSR')
    ylabel('IIS rate')
    set(gca, 'TickLength', [0 0])
    box off
    title('BSR vs. IIS')
    [r, pval] = corr(cdata(:, 2), cdata(:, 1), 'Type', 'Spearman');
    y = ylim;
    x = xlim;
    text(x(1), y(2), sprintf('r = %.2f, p = %.3f', r, pval))
    
    % iis vs delta
    subplot(3, 4, 7)
    hold on
    for i = 1 : length(nfiles)
        scatter(as.dband(i, :), as.iis(i, :), 50, '.')
    end
    xlabel('delta power')
    ylabel('IIS rate')
    set(gca, 'TickLength', [0 0])
    box off
    title('Delta vs. IIS')
    [r, pval] = corr(cdata(:, 3), cdata(:, 1), 'Type', 'Spearman');
    y = ylim;
    x = xlim;
    text(x(1), y(2), sprintf('r = %.2f, p = %.3f', r, pval))
    
    % bsr vs delta
    subplot(3, 4, 11)
    hold on
    for i = 1 : length(nfiles)
        scatter(as.dband(i, :), as.bsr(i, :), 50, '.')
    end
    xlabel('delta power')
    ylabel('BSR')
    set(gca, 'TickLength', [0 0])
    box off
    title('Delta vs. BSR')
    [r, pval] = corr(cdata(:, 3), cdata(:, 2), 'Type', 'Spearman');
    y = ylim;
    x = xlim;
    text(x(1), y(2), sprintf('r = %.2f, p = %.3f', r, pval))
    
    %  nspks histogram
    subplot(3, 4, 4)
    h = histogram(as.nspks, nbins, 'Normalization', 'count');
    h.EdgeColor = 'none';
    h.FaceColor = 'k';
    h.FaceAlpha = 0.4;
    hold on
    h = histogram(as.nspksBS, nbins, 'Normalization', 'count');
    h.EdgeColor = 'none';
    h.FaceColor = 'b';
    h.FaceAlpha = 0.4;
    legend({'All', 'in BS'})
    xlabel('Number of Spikes [#)]')
    ylabel('Mice [#]')
    set(gca, 'TickLength', [0 0])
    box off
    title('IIS')
    
    %  duration histogram
    subplot(3, 4, 8)
    h = histogram(as.recDur / fs / 60, nbins, 'Normalization', 'count');
    h.EdgeColor = 'none';
    h.FaceColor = 'k';
    h.FaceAlpha = 0.4;
    hold on
    h = histogram(as.deepDur / fs / 60, nbins, 'Normalization', 'count');
    h.EdgeColor = 'none';
    h.FaceColor = 'b';
    h.FaceAlpha = 0.4;
    legend({'All', 'BS'})
    xlabel('Duration [m)]')
    ylabel('Mice [#]')
    set(gca, 'TickLength', [0 0])
    box off
    title('Duration')
    
    %  threshold histogram
    subplot(3, 4, 12)
    h = histogram(as.thr, nbins, 'Normalization', 'count');
    h.EdgeColor = 'none';
    h.FaceColor = 'k';
    h.FaceAlpha = 0.4;
    xlabel('Threshold [mV)]')
    ylabel('Mice [#]')
    set(gca, 'TickLength', [0 0])
    box off
    title('Threshold')
   
    linkaxes([sp1, sp2, sp3], 'x');
    
    if saveFig
        export_fig(grpname, '-tif', '-transparent')
    end
end
end

