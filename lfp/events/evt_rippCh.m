function ch = evt_rippCh(basepath, basename, session)
% EVT_RIPPCH Resolve the ripple/detection channel from the ripple output.
%
%   ch = EVT_RIPPCH(basepath, basename, session)
%
%   SUMMARY:
%       The consumer-side channel resolver. Returns the channel the ripple
%       pipeline detected on, read from <basename>.ripp.mat (ripp.info.rippCh) -
%       the single source of truth - so downstream tools (ED, the curation GUI)
%       follow the ripple analysis instead of the session's Ripple tag. Falls
%       back to the session tag (evt_pickCh), then channel 1, so it stays
%       runnable before ripples have been detected.
%
%   INPUTS:
%       basepath - (Char)   Session directory.
%       basename - (Char)   File stem. {folder name}
%       session  - (Struct) Session metadata, used only for the fallback.
%                           {loaded on demand}
%
%   OUTPUTS:
%       ch       - (Num)    Zero-indexed detection channel(s).
%
%   DEPENDENCIES:
%       evt_pickCh; basepaths2vars (only for the fallback when session is empty).
%
%   HISTORY:
%       Created: 260717 (decouple consumers from channelTags.Ripple - read the
%                channel from ripp.info.rippCh instead).

if nargin < 2 || isempty(basename), [~, basename] = fileparts(basepath); end
if nargin < 3, session = []; end

% Primary: the channel the ripple pipeline recorded
f = fullfile(basepath, [basename, '.ripp.mat']);
if isfile(f)
    S = load(f, 'ripp');
    if isfield(S, 'ripp') && isfield(S.ripp, 'info') ...
            && isfield(S.ripp.info, 'rippCh') && ~isempty(S.ripp.info.rippCh)
        ch = S.ripp.info.rippCh;
        return;
    end
end

% Fallback: the session Ripple tag (load the session only if not supplied)
if isempty(session)
    v = basepaths2vars('basepaths', {basepath}, 'vars', {'session'}, ...
        'flgPrnt', false);
    if isfield(v, 'session'), session = v.session; end
end
ch = evt_pickCh(session);

% channelTags.Ripple may be the CellExplorer struct form ({.channels}); unwrap
% it so a fallback still returns a plain channel vector.
if isstruct(ch) && isfield(ch, 'channels'), ch = ch.channels; end

end     % EOF
