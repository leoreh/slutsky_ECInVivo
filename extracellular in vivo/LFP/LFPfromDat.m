function LFPfromDat


% loads lfp data. can specify channels, intervals, average across channels,
% resample, invert and more.
%  
% INPUT
%   basename    string. filename of lfp file. if empty retrieved from
%               basepath. if .lfp should not include extension, if .wcp
%               should include extension
%   basepath    string. path to load filename and save output {pwd}
%   extension   load from {'lfp'} (neurosuite), 'abf', 'wcp', or 'dat'.
%   forceL      logical. force reload {false}.
%   fs          numeric. requested sampling frequency {1250}
%   interval    numeric mat. list of intervals to read from lfp file [s]
%               can also be an interval of traces from wcp
%   ch          vec. channels to load
%   pli         logical. filter power line interferance (1) or not {0}
%   dc          logical. remove DC component (1) or not {0}
%   invertSig   logical. invert signal s.t. max is positive
%   saveVar     save variable {1}.
%   chavg       cell. each row contain the lfp channels you want to average
%   
% DEPENDENCIES
%   import_wcp
%   ce_LFPfromDat (if extension = 'dat')
% 
% OUTPUT
%   lfp         structure with the following fields:
%   fs
%   fs_orig
%   extension
%   interval    
%   duration    
%   chans
%   timestamps 
%   data  
% 
% 01 apr 19 LH & RA
% 19 nov 19 LH          load mat if exists  
% 14 jan 19 LH          adapted for wcp and abf 
%                       resampling
%
% TO DO LIST
%       # lfp from dat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'basename', '');
addOptional(p, 'extension', 'lfp');
addOptional(p, 'forceL', false, @islogical);
addOptional(p, 'fs', 1250, @isnumeric);
addOptional(p, 'interval', [0 inf], @isnumeric);
addOptional(p, 'ch', [1 : 16], @isnumeric);
addOptional(p, 'pli', false, @islogical);
addOptional(p, 'dc', false, @islogical);
addOptional(p, 'invertSig', false, @islogical);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'chavg', {}, @iscell);

parse(p,varargin{:})
basepath = p.Results.basepath;
basename = p.Results.basename;
extension = p.Results.extension;
forceL = p.Results.forceL;
fs = p.Results.fs;
interval = p.Results.interval;
ch = p.Results.ch;
pli = p.Results.pli;
dc = p.Results.dc;
invertSig = p.Results.invertSig;
saveVar = p.Results.saveVar;
chavg = p.Results.chavg;

fsOut = 1250;
[~, basename] = fileparts(basepath);

import iosr.dsp.*


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end



% EOF