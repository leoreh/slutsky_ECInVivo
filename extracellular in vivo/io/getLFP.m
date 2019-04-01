function lfp = getLFP(varargin)

% gets video struct from .xls of ToxTrac (Rodriguez et al., 2018)
%  
% INPUT
%   filename    string. filename of lfp file
%   basepath    string. path to load filename and save output {pwd}
%   sheet       string. sheet to extract data from {'Seq_0_Arena_1_1'}
%   graphics    plot figure {1}.
%   saveFig     save figure {1}.
%   saveVar     save variable {1}.
% 
% OUTPUT
%   vid         structure with the following fields:
%       fps         frames per second of video file
%       nframes     number of frames analyzed
%       drate       detection rate of object (mouse)
%       t           frame timestamps
%       x           x coordinates (pixels)
%       y           y coordinates (pixels)
%       spd         instantaneous speed. starts at third frame
%       acc         instantaneous acceleration. starts at fifth frame
% 
% 
% 01 apr 19 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'filename');
addOptional(p, 'sheet', 'Seq_0_Arena_1_1', @ischar);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'saveFig', false, @islogical);
addOptional(p, 'saveVar', true, @islogical);

parse(p,varargin{:})
basepath = p.Results.basepath;
sheet = p.Results.sheet;
graphics = p.Results.graphics;
saveFig = p.Results.saveFig;
saveVar = p.Results.saveVar;

if isempty(intervals)
    intervals = [0 inf];
end

cd(basepath)
if isempty(filename)
    [~, filename] = fileparts(basepath);
    filename = [filename '.lfp'];
end

% check if file exists
if ~exist(filename)
error('file %s does not exist', filename)
end


end
