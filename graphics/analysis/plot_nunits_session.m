function plot_nunits_session(varargin)

% plot a stacked bar plot of the number of units per spike group
%
% INPUT:
%   basepath        char. path of recording
%   frBoundries     2 x 2 mat. include only units with mfr within
%                   frBoundries. 1st row refers to rs, 2nd to fs
%   saveFig         logical
% 
% OUTPUT
%
% DEPENDENCIES
% 
% TO DO LIST
%
% 12 jan 22 LH      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'frBoundries', [0 Inf; 0 Inf], @isnumeric);
addParameter(p, 'saveFig', true, @islogical);

parse(p, varargin{:})
basepath        = p.Results.basepath;
frBoundries     = p.Results.frBoundries;
saveFig         = p.Results.saveFig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, basename] = fileparts(basepath);

% load vars from each session
varsFile = ["fr"; "spikes"; "cell_metrics"; "datInfo"; "session"];
varsName = ["fr"; "spikes"; "cm"; "datInfo"; "session"];
if ~exist('v', 'var')
    v = getSessionVars('basepaths', {basepath}, 'varsFile', varsFile,...
        'varsName', varsName);
end

% plot only su or all units
suFlag = true;                 

% spike groups
grp = unique(v.spikes.shankID);

% select units
units_grp = zeros(length(grp), 3);
for igrp = 1 : length(grp)
    units_grp(igrp, 1 : 2) = sum(selectUnits('basepath', basepath, 'spikes', v.spikes,...
        'cm', v.cm, 'fr', v.fr, 'suFlag', suFlag, 'grp', igrp,...
        'frBoundries', frBoundries, 'forceA', true, 'saveVar', false), 2)';
    units_grp(igrp, 3) = sum(v.spikes.shankID == igrp) - sum(units_grp(igrp, :));
end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMatlabGraphics(false)
clr = [0.4 0.7 1; 0.95 0.4 0.4; 0.66 0.4 0.66];
fh = figure;

subplot(1, 2, 1)
bh = bar(units_grp, 'stacked');
bh(1).FaceColor = clr(1, :);
bh(2).FaceColor = clr(2, :);
bh(3).FaceColor = clr(3, :);
legend({'RS', 'FS', 'Other'})
xticks(1 : length(grp))
xlabel('Spike Group')
ylabel('No. Units')
box off

subplot(1, 2, 2)
bh = bar(sum(units_grp, 1), 'FaceColor', 'flat');
bh.CData = clr;
ylabel('No. Units')
xticklabels({'RS', 'FS', 'Other'})

sgtitle(basename)

if saveFig
    figpath = fullfile(basepath, 'graphics');
    mkdir(figpath)
    figname = fullfile(figpath, [basename, '_nunits']);
    export_fig(figname, '-tif', '-transparent', '-r300')
end

end