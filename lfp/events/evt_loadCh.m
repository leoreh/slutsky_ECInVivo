function [lfp, fs] = evt_loadCh(basepath, basename, session, ch, win, bit2uv)
% EVT_LOADCH Load one averaged LFP channel from the binary .lfp file.
%
%   [lfp, fs] = EVT_LOADCH(basepath, basename, session, ch, win, bit2uv)
%
%   SUMMARY:
%       Shared raw-channel loader for both signal-loading paths (ripp_sigLoad
%       and the 'lfp' source of ed_sigLoad). Reads the requested channel(s) from
%       <basename>.lfp over the analysis window, averages them when more than
%       one is given, and applies the bit-to-microvolt conversion. The
%       conversion auto-detects the acquisition system from the wideband rate
%       (TDT vs Intan) unless supplied. Returns a column vector in microvolts.
%
%   INPUTS:
%       basepath - (Char)   Session directory holding <basename>.lfp.
%       basename - (Char)   File stem.
%       session  - (Struct) Session metadata; reads extracellular.srLfp,
%                           .nChannels, and .sr (for the conversion default).
%       ch       - (Num)    Zero-indexed channel id(s) to load.
%       win      - (Vec)    Window [start end] (s); Inf end spans to file end.
%       bit2uv   - (Num)    Conversion factor. Empty -> auto (TDT/Intan).
%
%   OUTPUTS:
%       lfp      - (Vec)    [n x 1] Windowed signal (microvolts).
%       fs       - (Num)    LFP sampling frequency [Hz] (extracellular.srLfp).
%
%   DEPENDENCIES:
%       binary_load.
%
%   HISTORY:
%       Created: 260706 (absorbs the duplicated wrapper channel-load blocks).

fs     = session.extracellular.srLfp;
nChans = session.extracellular.nChannels;

% Voltage conversion (auto: TDT vs Intan) unless supplied
if isempty(bit2uv)
    if round(session.extracellular.sr) == 24414
        bit2uv = 1;         % TDT / Tucker
    else
        bit2uv = 0.195;     % Intan
    end
end

% Window -> load duration (Inf spans to the file end)
if isinf(win(2))
    sigDur = Inf;
else
    sigDur = win(2) - win(1);
end

% Load the requested channel(s); average when more than one
fname = fullfile(basepath, [basename, '.lfp']);
lfp = double(binary_load(fname, 'duration', sigDur, 'fs', fs, ...
    'nCh', nChans, 'start', win(1), 'ch', ch, ...
    'downsample', 1, 'bit2uv', bit2uv));
if size(lfp, 2) > 1
    lfp = mean(lfp, 2);
end
lfp = lfp(:);

end     % EOF
