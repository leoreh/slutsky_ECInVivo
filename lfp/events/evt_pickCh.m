function ch = evt_pickCh(session)
% EVT_PICKCH Resolve the detection channel from the session's ripple tag.
%
%   ch = EVT_PICKCH(session)
%
%   SUMMARY:
%       Shared channel resolution for both event pipelines. Returns the
%       ripple-tagged channel so ED and ripples detect on the same LFP. Falls
%       back to channel 1 (with a warning) when the session, its channelTags,
%       or the Ripple tag are absent - keeping detection runnable on minimal
%       session layouts.
%
%   INPUTS:
%       session - (Struct) Session metadata; reads channelTags.Ripple. May be
%                          empty or lack the field.
%
%   OUTPUTS:
%       ch      - (Num)    Zero-indexed detection channel.
%
%   DEPENDENCIES:
%       None.
%
%   HISTORY:
%       Created: 260706 (absorbs the duplicated wrapper channel-pick blocks).

ch = 1;
if ~isempty(session) && isfield(session, 'channelTags') ...
        && isfield(session.channelTags, 'Ripple') ...
        && ~isempty(session.channelTags.Ripple)
    ch = session.channelTags.Ripple;
else
    warning('evt_pickCh:noRippleTag', ...
        'no channelTags.Ripple; detection uses channel %d.', ch);
end

end     % EOF
