function call_powershell_cp_basepath(varargin)
% CALL_POWERSHELL_CP_BASEPATH Calls the Copy-BasePath.ps1 PowerShell script.
%
% SUMMARY:
% This function serves as a MATLAB wrapper to call the PowerShell script
% 'Copy-BasePath.ps1', passing the necessary arguments.
%
% INPUT (Key-Value Pairs, same as the original cp_basepath and PowerShell script):
%   fNames       (cell array) Cell array of file patterns.
%   basepath     (char) Full path to the session directory.
%   newpath      (char, optional) New base path. Defaults to 'D:\Data' in PowerShell.
%   overwrite    (logical, optional) Flag to overwrite. Defaults to false.
%   verbose      (logical, optional) Flag for verbose output. Defaults to false.
%
% EXAMPLE:
%   fNames = {'*.clu.1', '*.res.1'};
%   basepath = 'E:\Data\lh140\lh140_230624_090049';
%   call_powershell_cp_basepath('fNames', fNames, 'basepath', basepath, 'verbose', true);
%
%   call_powershell_cp_basepath('fNames', {'*.session.mat'}, ...
%                               'basepath', 'C:\Test\mouse1\mouse1_sess1', ...
%                               'newpath', 'D:\BackupData', ...
%                               'overwrite', true, 'verbose', true);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % arguments using inputParser
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p = inputParser;
    % Match the parameter names used in your PowerShell script
    % For FilePatterns, PowerShell expects an array.
    addParameter(p, 'fNames', {}, @(x) iscell(x) && all(cellfun(@ischar, x)));
    addParameter(p, 'basepath', '', @ischar);
    addParameter(p, 'newpath', '', @ischar); % Let PowerShell handle default if empty
    addParameter(p, 'overwrite', false, @islogical);
    addParameter(p, 'verbose', false, @islogical);

    parse(p, varargin{:});

    matlab_fNames    = p.Results.fNames;
    matlab_basepath  = p.Results.basepath;
    matlab_newpath   = p.Results.newpath;
    matlab_overwrite = p.Results.overwrite;
    matlab_verbose   = p.Results.verbose;

    % --- Construct the PowerShell command ---

    % Path to your PowerShell script.
    % Adjust this path if the script is not in the current MATLAB directory
    % or a directory on the system PATH.
    % You can use mfilename('fullpath') if this wrapper is in the same dir
    % as the .ps1 file to make it more robust.
    % For example:
    % scriptFolder = fileparts(mfilename('fullpath'));
    % psScriptPath = fullfile(scriptFolder, 'Copy-BasePath.ps1');
    psScriptPath = 'D:\Code\Personal\powershell\Copy-BasePath.ps1'; % Assuming it's in the current directory

    % Base PowerShell command
    % -NoProfile: Speeds up startup by not loading user profiles.
    % -ExecutionPolicy Bypass: Allows running unsigned scripts for this session.
    % -File: Specifies the script file to run.
    command = sprintf('powershell.exe -ExecutionPolicy Bypass -NoProfile -File "%s"', psScriptPath);

    % Add FilePatterns argument
    % PowerShell string arrays from command line can be passed as comma-separated,
    % single-quoted strings: e.g., -FilePatterns '*.txt','*.log'
    if isempty(matlab_fNames)
        error('fNames parameter is mandatory and cannot be empty.');
    end
    % Escape single quotes within patterns themselves (e.g. if a pattern is 'file''s.txt')
    % and then wrap each pattern in single quotes.
    fNamesFormatted = cellfun(@(x) ['''' strrep(x, '''', '''''') ''''], matlab_fNames, 'UniformOutput', false);
    fNamesStr = strjoin(fNamesFormatted, ',');
    command = [command, sprintf(' -FilePatterns %s', fNamesStr)]; % No outer quotes for the value itself

    % Add BasePath argument
    if isempty(matlab_basepath)
        error('basepath parameter is mandatory and cannot be empty.');
    end
    % Paths with spaces must be quoted
    command = [command, sprintf(' -BasePath "%s"', matlab_basepath)];

    % Add NewPath argument (optional)
    if ~isempty(matlab_newpath)
        command = [command, sprintf(' -NewPath "%s"', matlab_newpath)];
    end

    % Add Overwrite switch (optional)
    if matlab_overwrite
        command = [command, ' -Overwrite'];
    end

    % Add Verbose common parameter (optional)
    % The PowerShell script uses [CmdletBinding()] which enables -Verbose
    if matlab_verbose
        command = [command, ' -Verbose'];
        disp(['Executing PowerShell: ', command]); % Display command if MATLAB verbose is on
    end

    % --- Execute the PowerShell command ---
    [status, cmdout] = system(command);

    % --- Handle output and status ---
    if matlab_verbose || status ~= 0 % Show output if MATLAB verbose or if there was an error
        disp('-------------------- PowerShell Output Begin --------------------');
        disp(cmdout);
        disp('--------------------- PowerShell Output End ---------------------');
    end

    if status ~= 0
        warning('MATLAB:call_powershell_cp_basepath:ExecutionError', ...
                'PowerShell script execution failed with status %d. See output above for details.', status);
    elseif matlab_verbose
        disp('PowerShell script executed successfully.');
    end

end