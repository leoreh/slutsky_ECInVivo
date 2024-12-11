function ied = curate(ied,varargin)
% open manual curation window & wait for user to close it.
% save IED upon request.
%
%   INPUT (in 1st position):
%       IED.data object, after detection.
%   INPUT (optional, name value):
%       saveVar     logical {true}. save variable
%       basepath    recording session path {pwd}
%       basename    string. if empty extracted from basepath
%   **IMPORTANT:** if saving & basepath / basename is given, they take
%                  precedence. File will be save as they define,
%                  ied.file_loc will be overwritten.
% OUTPUT
%   IED.data object

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isa(ied,'IED.data')
    error("first input must be IED.data obj")
end

p = inputParser;
addParameter(p, 'saveVar', true, @islogical);
addParameter(p, 'basepath', pwd, @isstr);
addParameter(p, 'basename', [], @isstr);
parse(p, varargin{:})

saveVar = p.Results.saveVar;
basepath = p.Results.basepath;
basename = p.Results.basename;

if any(~ismember(["basepath","basename"],p.UsingDefaults)) && saveVar
    % if user gave any info & want to save, overwrite existing
    if isempty(basename)
        [~,basename] = fileparts(basepath);
    end
    ied.file_loc = fullfile(basepath,join([basename "ied.mat"],"."));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open windows, let user mark
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cur_win = IED_curation_GUI("ied",ied);
waitfor(cur_win.UIFigure)

% decide analysis stage - "during_curation" or "curated"
if ied.last_mark == numel(ied.pos)
    ied.status = "curated";
else
    ied.status = "during_curation";
end

% update git info
try
    ied.git_last_step = get_gits_status(["slutsky_ECInVivo", "CellExplorer"]);
catch err
    ied.git_last_step = join(["Error: " err.message]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save result if needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if saveVar
    save(ied.file_loc,"ied")
    fprintf("\n****** Save in %s ******\n",ied.file_loc)
else
    ied.file_loc = '';
end

% EOF