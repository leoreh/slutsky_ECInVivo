function otl = get_otlSpec(varargin)

% gets spectrogram outliers, typically to exclude them as artifacts.
% must recalculate spectrogram to make sure window isn't smoothed.
% uses eeg signal form accusleep

% INPUT
%   basepath    char. fullpath to recording folder {pwd} ch
%   ch          numeric vector depicting channels over which
%               the spectrogram should be calculated
%   thrFactor   threshold in stds  / mad for defining outliers {6}
%   graphics    logical. plot figure form get_otl {false}
%   saveVar     logical. organize and save struct {true}
%   flgForce    logical. force analysis or load if exists already {false}
%   flgInspect  logical. manually inspect outliers with accusleep.
%               this offers the chance to manually curate them.
%
% OUTPUT
%   otl         struct
%
% TO DO LIST
%
% 26 mar 24 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'basepath', pwd, @ischar)
addParameter(p, 'flgForce', false, @islogical)
addParameter(p, 'flgInspect', false, @islogical)
addParameter(p, 'thrFactor', 6, @isnumeric)
addParameter(p, 'graphics', false, @islogical)
addParameter(p, 'saveVar', true, @islogical)

parse(p, varargin{:})
basepath        = p.Results.basepath;
flgForce        = p.Results.flgForce;
flgInspect      = p.Results.flgInspect;
thrFactor       = p.Results.thrFactor;
saveVar         = p.Results.saveVar;
graphics        = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% files
[~, basename] = fileparts(basepath);
fname_sig = fullfile(basepath, [basename, '.sleep_sig.mat']);
fname_states = fullfile(basepath, [basename, '.sleep_states.mat']);
fname_otl = fullfile(basepath, [basename, '.otl_spec.mat']);
fname_labels = fullfile(basepath, [basename, '.otl_labels.mat']);

% load if exists
if ~flgForce && exist(fname_otl, 'file')
    load(fname_otl, 'otl')
    return
end

% load / re-calc spec
sig = load(fname_sig, 'eeg');
sig = sig.eeg;

spec = calc_spec('sig', sig, 'fs', 1250, 'graphics', false, 'saveVar', false,...
    'padfft', 0, 'winstep', 1, 'logfreq', false, 'ftarget', [],...
    'ch', {[1]}, 'force', true, 'window', 1);
s = spec.s;
freq = spec.freq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detection of spec outliers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grab only a specific frequency range to consider
[~, fidx1] = min(abs(freq - 15));
[~, fidx2] = min(abs(freq - 47));
fidx = fidx1 : fidx2;

data = s;
otl = get_otl(data, 'thrFactor', thrFactor, 'graphics', graphics);

if saveVar
    save(fname_otl, 'otl')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flgInspect

    % load data
    sSig = load(fname_sig);
    load(fname_states, 'ss')
    cfg = as_loadConfig();

    % insert outliers to "bin" state
    binIdx = ss.labels == cfg.nstates + 1;
    if any(binIdx)
        warning('bin state already exists')
    end
    ss.labels(otl.boolean) = cfg.nstates + 1;
    AccuSleep_viewer(sSig, ss.labels, fname_labels)
    waitfor(gcf)

    % load manual labels and update otl struct
    if exist(fname_labels, 'file')
        load(fname_labels)

        otl.idx = find(labels == cfg.nstates + 1);
        otl.boolean = false(size(data, 1), 1);
        otl.boolean(otl.idx) = true;
        otl.bouts = binary2bouts('vec', otl.boolean, 'minDur', [], 'maxDur', [],...
            'interDur', [], 'exclude', false, 'flgPrnt', false);

        if saveVar
            save(fname_otl, 'otl')
        end
    end
end


% EOF

