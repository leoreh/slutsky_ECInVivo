function [emg, spec, boutTimes] = recRep_wrapper(varargin)

% wrapper for creating a recording representative (see
% recRep_plot).

% creates a figure where each panel (row) is a different plot across time.
% the relation between panels can be determined such that, for example, the
% buttom panel is a zoom-in view of the top panel.
%
% INPUT:
%   basepath    char. path to session folder {pwd}
%   panels2plot string array. determines which plots to include and in what
%               order. can be "spec", "emg", "hypnogram", "raster", "raw" 
%   tlayout     vector of 2 elements descrbigin the tiled layout
%   panelTiles  cell describing the spread of panels across tiles
%   th          handle of tile layout
%   xDur        numeric vector where each element describes the x-range for
%               the corrsponding panel
%   xStart      numeric vector where each element describes time point
%               where the panel starts, ie xlim(1)
%   xRelations  cell of zoom-in relations (see recRep_wrapper)
%   tLine       numeric. time point in seconds to add a dashed line
%   rawCh       numeric vector describing the channels to load for raw data
%   saveFig     logical
%
% OUTPUT
% 
% CALLS
%
% TO DO LIST
%
% 05 mar 24 LH  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'panels2plot', []);
addParameter(p, 'tlayout', []);
addParameter(p, 'panelTiles', []);
addParameter(p, 'th', []);
addParameter(p, 'xDur', []);
addParameter(p, 'xStart', []);
addParameter(p, 'xRelations', []);
addParameter(p, 'tLine', []);
addParameter(p, 'rawCh', []);
addParameter(p, 'saveFig', false, @islogical);

parse(p, varargin{:})
basepath        = p.Results.basepath;
panels2plot     = p.Results.panels2plot;
tlayout         = p.Results.tlayout;
panelTiles      = p.Results.panelTiles;
th              = p.Results.th;
xDur            = p.Results.xDur;
xStart          = p.Results.xStart;
xRelations      = p.Results.xRelations;
tLine           = p.Results.tLine;
rawCh           = p.Results.rawCh;
saveFig         = p.Results.saveFig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine panel arrangement
npanels = length(panels2plot);

% downsample emg / eeg signals
fs_emg = 250;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[emg, spec, boutTimes] = recRep_prep('basepath', basepath,...
    'panels2plot', panels2plot, 'fs_emg', fs_emg);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

recRep_plot('basepath', basepath,...
    'xDur', xDur, 'xStart', xStart, 'panels2plot', panels2plot,...
    'spec', spec, 'emg', emg, 'boutTimes', boutTimes,...
    'saveFig', saveFig, 'fs_emg', fs_emg, 'rawCh', rawCh,...
    'xRelations', xRelations, 'tlayout', tlayout, 'panelTiles', panelTiles,...
    'tLine', tLine, 'th', th);

end

% EOF

