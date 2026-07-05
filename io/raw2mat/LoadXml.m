function [xml, rxml] = LoadXml(fbasename, varargin)
% LOADXML  Parse a NeuroScope / ndManager .xml into a LoadPar-compatible struct.
%
% Native replacement for the previous xmltools-based parser: uses MATLAB's
% built-in xmlread (DOM) instead of the hand-rolled xmltools, so it is no
% longer sensitive to tag formatting (multi-attribute tags, etc.). The output
% struct is field-for-field identical to the legacy LoadXml and remains
% backwards compatible with LoadPar-style downstream code.
%
% INPUT
%   fbasename   path to the .xml (with or without the .xml extension)
%
% OUTPUT
%   xml         struct with:
%                 .FileName, .Date
%                 .nBits .nChannels .SampleRate .SampleTime .VoltageRange
%                 .Amplification .Offset .lfpSampleRate
%                 .AnatGrps(g).Channels/.Skip     (anatomical groups)
%                 .SpkGrps(g).Channels/.nSamples/.PeakSample/.nFeatures
%                 .nElecGps .ElecGp{g}            (spike groups, if present)
%                 .HiPassFreq                     (if a process_mhipass program exists)
%   rxml        the raw parsed DOM document (kept for signature compatibility)
%
% Channel numbers are 0-based, matching the xml (callers add 1 where needed).

xml = struct;

if ~contains(fbasename, '.xml')
    fbasename = [fbasename '.xml'];
end

rxml = xmlread(fbasename);
root = rxml.getDocumentElement;          % <parameters>

xml.FileName = fbasename;

% general info
gi = getChild(root, 'generalInfo');
if ~isempty(gi)
    d = getChild(gi, 'date');
    if ~isempty(d), xml.Date = nodeText(d); end
end

% acquisition system
acq = getChild(root, 'acquisitionSystem');
if ~isempty(acq)
    xml.nBits         = getNum(acq, 'nBits');
    xml.nChannels     = getNum(acq, 'nChannels');
    xml.SampleRate    = getNum(acq, 'samplingRate');
    xml.SampleTime    = 1e6 / xml.SampleRate;   % backwards compatible
    xml.VoltageRange  = getNum(acq, 'voltageRange');
    xml.Amplification = getNum(acq, 'amplification');
    xml.Offset        = getNum(acq, 'offset');
end

% field potentials (lfp sampling rate)
fp = getChild(root, 'fieldPotentials');
if ~isempty(fp)
    xml.lfpSampleRate = getNum(fp, 'lfpSamplingRate');
end

% anatomical groups (with per-channel skip flag)
ad = getChild(root, 'anatomicalDescription');
if ~isempty(ad)
    groups = getChildren(getChild(ad, 'channelGroups'), 'group');
    for g = 1:numel(groups)
        chans = getChildren(groups{g}, 'channel');
        for c = 1:numel(chans)
            xml.AnatGrps(g).Channels(c) = str2double(nodeText(chans{c}));
            xml.AnatGrps(g).Skip(c)     = str2double(char(chans{c}.getAttribute('skip')));
        end
    end
end

% spike groups
sd = getChild(root, 'spikeDetection');
if ~isempty(sd)
    groups = getChildren(getChild(sd, 'channelGroups'), 'group');
    if isempty(groups)
        xml.nElecGps = 0;
    else
        for g = 1:numel(groups)
            chans = getChildren(getChild(groups{g}, 'channels'), 'channel');
            for c = 1:numel(chans)
                xml.SpkGrps(g).Channels(c) = str2double(nodeText(chans{c}));
            end
            if ~isempty(getChild(groups{g}, 'nSamples'))
                xml.SpkGrps(g).nSamples   = getNum(groups{g}, 'nSamples');
                xml.SpkGrps(g).PeakSample = getNum(groups{g}, 'peakSampleIndex');
                xml.SpkGrps(g).nFeatures  = getNum(groups{g}, 'nFeatures');
            end
            xml.nElecGps  = numel(groups);          % backwards compatibility
            xml.ElecGp{g} = xml.SpkGrps(g).Channels;
        end
    end
end

% high-pass frequency from a process_mhipass program
pr = getChild(root, 'programs');
if ~isempty(pr)
    programs = getChildren(pr, 'program');
    for p = 1:numel(programs)
        nm = getChild(programs{p}, 'name');
        if ~isempty(nm) && strcmp(nodeText(nm), 'process_mhipass')
            params = getChildren(getChild(programs{p}, 'parameters'), 'parameter');
            for q = 1:numel(params)
                pn = getChild(params{q}, 'name');
                if ~isempty(pn) && strcmp(nodeText(pn), 'frequency')
                    xml.HiPassFreq = str2double(nodeText(getChild(params{q}, 'value')));
                    break
                end
            end
        end
    end
end

end

% ------------------------------------------------------------------------
function e = getChild(node, tag)
% first direct child element named tag, or [] if none
e = [];
if isempty(node), return; end
ch = node.getChildNodes;
for k = 0:ch.getLength - 1
    c = ch.item(k);
    if c.getNodeType == 1 && strcmp(char(c.getNodeName), tag)   % 1 = ELEMENT_NODE
        e = c; return
    end
end
end

function es = getChildren(node, tag)
% all direct child elements named tag, as a cell array
es = {};
if isempty(node), return; end
ch = node.getChildNodes;
for k = 0:ch.getLength - 1
    c = ch.item(k);
    if c.getNodeType == 1 && strcmp(char(c.getNodeName), tag)
        es{end+1} = c; %#ok<AGROW>
    end
end
end

function t = nodeText(e)
% trimmed text content of a leaf element
t = strtrim(char(e.getTextContent));
end

function v = getNum(parent, tag)
% numeric value of the first child element named tag ([] if absent)
e = getChild(parent, tag);
if isempty(e), v = []; return; end
v = str2double(nodeText(e));
end
