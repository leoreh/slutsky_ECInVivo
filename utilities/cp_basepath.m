function cp_basepath(varargin)
% CP_BASEPATH Copies files from a session directory to a new location.
%
% This function copies files from a single session directory to a new location
% while maintaining the mouse/session directory structure. For example, files
% from 'E:\Data\lh140\lh140_230624_090049\*.clu.1' will be copied to
% 'D:\Data\lh140\lh140_230624_090049\*.clu.1'. Uses PowerShell's robocopy
% for faster file operations.
%
% INPUT (Optional Key-Value Pairs):
%   fNames       (cell array) Cell array of file patterns to copy. Can include
%                wildcards (e.g., '*.clu.1', '*.session.mat*'). Each pattern
%                will be expanded to match all files in the session directory.
%   basepath     (char) Full path to the session directory (e.g.,
%                'E:\Data\lh140\lh140_230624_090049')
%   newpath      (char) New base path where files will be copied to. If empty,
%                uses 'D:\Data'. {[]}
%   overwrite    (logical) Flag to overwrite existing files. {false}
%   verbose      (logical) Flag to display progress information. {false}
%
% OUTPUT:
%   None. Files are copied to the new location.
%
% EXAMPLE:      See end of function
%
% DEPENDENCIES:
%   None
%
% HISTORY:
%   Aug 2024 (AI Assisted) - Created function with input parsing and
%                            enhanced documentation
%             - Added support for wildcard patterns in filenames
%             - Added verbose option for progress display
%             - Switched to PowerShell robocopy for faster file operations
%             - Optimized to use single robocopy call for all files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments using inputParser
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'fNames', @(x) iscell(x) && all(cellfun(@ischar, x)));
addOptional(p, 'basepath', @ischar);
addOptional(p, 'newpath', 'D:\Data', @ischar);
addOptional(p, 'overwrite', false, @islogical);
addOptional(p, 'verbose', false, @islogical);

parse(p, varargin{:});
fNames = p.Results.fNames;
basepath = p.Results.basepath;
newpath = p.Results.newpath;
overwrite = p.Results.overwrite;
verbose = p.Results.verbose;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract mouse and session names from basepath
[~, sessionName] = fileparts(basepath);
[~, mouseName] = fileparts(fileparts(basepath));

% Create new session directory
sessionPath = fullfile(newpath, mouseName, sessionName);
if ~exist(sessionPath, 'dir')
    mkdir(sessionPath);
    if verbose
        fprintf('Created: %s\n', sessionPath);
    end
end

% Prepare robocopy command
robocopyOpts = '/NP /NJH /NJS';  % No progress, no job header, no job summary
if ~overwrite
    robocopyOpts = [robocopyOpts ' /XC /XN /XO'];  % Exclude existing files if they are older/newer/changed (effectively skip if exists and not identical)
else
    robocopyOpts = [robocopyOpts ' /IS /IT']; % Include same files, include tweaked files (effectively overwrite)
end
if verbose
    robocopyOpts = [robocopyOpts ' /V'];  % Verbose output
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET FILE NAMES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Expand each file pattern and collect all matching files
cpFiles = {};
for iPattern = 1:length(fNames)
    % Get all files matching the pattern
    pattern = fullfile(basepath, fNames{iPattern});
    files = dir(pattern);

    % Add full paths to the list
    if ~isempty(files)
        cpFiles = [cpFiles; {files.name}'];
    end
end

% Remove duplicates (in case multiple patterns match the same file)
cpFiles = unique(cpFiles)';

% Join file names with spaces, ensuring each is quoted
quotedFileNames = cellfun(@(f) ['"' f '"'], cpFiles, 'UniformOutput', false);
filesArgString = strjoin(quotedFileNames, ' ');

if verbose
    fprintf('Found %d files to copy\n', length(cpFiles));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPY FILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Construct robocopy command
% The syntax is: robocopy <SourceDir> <DestDir> [File1 File2 ...] [Options]
cmd = sprintf('robocopy "%s" "%s" %s %s', ...
    basepath, ...           % source directory
    sessionPath, ...        % destination directory
    filesArgString, ...     % list of files to copy
    robocopyOpts);         % options

% Execute command
[status, cmdout] = system(cmd);

if verbose && ~isempty(cmdout)
    fprintf('Robocopy output:\n%s\n', cmdout);

    % Robocopy considers status codes 0-7 as success (or success with non-critical issues)
    % 0 = No files copied (all skipped or source empty)
    % 1 = All files copied successfully.
    % 2 = Extra files/dirs were detected. Examine log. No files were copied.
    % 3 = Some files were copied. Additional files were present. No failure was encountered.
    % etc. up to 7.
    if status < 8
        fprintf('Copy operation overall status: %d\n', status);
    else
        fprintf('Copy operation encountered errors. Status: %d\n', status);
    end
end

end     % EOF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE CALL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fNames = {...
%     '*.acceleration.mat',...
%     '*.datInfo*',...
%     '*.cell_metrics*',...
%     '*.fr.mat',...
%     '*.frEmg.mat',...
%     '*.nrs',...
%     '*.ripp.*',...
%     '*.session.mat*',...
%     '*.spikes.*',...
%     '*.spktimes.mat',...
%     '*.sr.mat',...
%     '*.st_*',...
%     '*.swv_',...
%     '*.units.mat',...
%     '*.xml',...
%     '*.sleep_states*',...
%     '*.sleep_labelsMan*',...
%     % '*.lfp',...
%     % '*.sleep_*',...
%     };
%
%
% % get group basepaths
% basepaths = mcu_basepaths('wt');
% nPaths = length(basepaths);
%
% newpath = 'D:\OneDrive - Tel-Aviv University\Data\Baclofen\Control';
%
% for iPath = 1 : nPaths
%     basepath = basepaths{iPath};
%     cp_basepath('fNames', fNames, 'basepath', basepath, 'newpath', newpath,...
%         'overwrite', true, 'verbose', true)
% end