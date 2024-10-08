function [stat_tab, git_repos] = get_gits_status(git_dirs)
    % GET_GITS_STATUS - Collects and summarizes the update status of specified Git repositories.
    %   Having git executable install and on path is recommended: https://www.git-scm.com/downloads
    % INPUTS:
    %   git_dirs    - (Optional) Text array specifying Git repository
    %                 full paths or folder names on path. 
    %                 Default: ["slutsky_ECInVivo", "CellExplorer"].
    %
    % OUTPUTS:
    %   stat_tab    - Table with Git repository status: full path, last commit ID, and modification.
    %                 
    %   git_repos   - Array of Git repository objects.
    %
    % EXAMPLES:
    %   dirs = ["path/to/repo1", "path/to/repo2"];
    %   [status_table, repositories] = get_gits_status(dirs);
    %
    % DETAILS:
    %   - If full path is not given, top folders are chosen in case of shadowing.
    %
    %   - Modification Status:
    %     Determined by the presence of modified or untracked files. If neither, the repository is
    %     considered unmodified (false).
    %
    % NOTES:
    %   - MATLAB 2023b is necessary to collect Git repository objects or when
    %     the Git executable is not installed externally.
    %
    %   - Use of External Git:
    %     The function attempts to use the external Git executable for faster and more informative results.
    %     If external Git is unavailable, it falls back to slower MATLAB-based Git commands.
    %
    %   - Modification Result Format:
    %     The `modified_res` column attempts to mimic the output of `git status --porcelain`.
    %     Note that when not using external Git this is not an exact match.
    %
    %   - If no output is requested, the function prints the result table to the
    %     command window without assigning it to a variable ('ans').
    %
    % See also: gitrepo
    %
    % Author: LdM https://github.com/Liordemarcas
    % Published: 240114

    arguments
        git_dirs (1,:) string {mustBeText} = ["slutsky_ECInVivo", "CellExplorer", "my_MATLAB_code"]; % never should actually use pwd, nargin will take precedence
    end

    %%%%%%%% find folders to work on
    % Ech input folder should be either a folder name name (than search on
    % path) or folder full path.
    % when using path, in case of shadowing, take the top folder.
    % If you can't find a folder, skip it (instead of error).

    % find path to folders
    current_subfolders = dir();
    current_subfolders(~[current_subfolders.isdir]) = [];
    current_subfolders = {current_subfolders.name};
    for iFold = numel(git_dirs):-1:1

        if isfolder(git_dirs(iFold)) && ~ismember(git_dirs(iFold), current_subfolders)
            % what was given is a full folder path, no need to search.
            % note that creating full path is needed if the name is a
            % current subfolder without full path.
            continue
        else
            % collect folder info from path
            fold_file_list = what(git_dirs(iFold));

            % collect path only if found
            if ~isempty(fold_file_list)
                % take only top folder in case of shadowing
                git_dirs(iFold) = string(fold_file_list(1).path);
            else
                % place missing for delete later
                git_dirs(iFold) = missing;
            end
        end
    end
    git_dirs(ismissing(git_dirs)) = [];

    % init status table
    [~, dirs_names] = fileparts(git_dirs);
    stat_tab = table('Size',[numel(git_dirs), 4],'VariableTypes', ["string", "string", "logical", "string"],...
        'VariableNames',["full_path","last_commit","modified","modified_res"],'RowNames',dirs_names);

    % create cleanup object, to return to base folder upon completion / error
    return_dir = pwd;
    cl = onCleanup(@() cd(return_dir));

    % collect repos and table summary of thier status
    for iFold = numel(git_dirs):-1:1
        % write path
        stat_tab.full_path(iFold) = git_dirs(iFold);

        % move to the folder, to use fast git commands
        cd(git_dirs(iFold))

        % collect last commit - require access to external git
        [stat, result] = system('git rev-parse --short HEAD');
        if stat ~= 0
            % for some reason git didn't work - try by MATLAB (way slower).
            % Requires MATLAB 2023b.
            git_repos{iFold} = gitrepo(git_dirs(iFold)); % get repo object.
            stat_tab.last_commit(iFold) = git_repos{iFold}.LastCommit.ID{1}(1:7); % save last commit

            % find if anything was modified
            mod_files       = git_repos{iFold}.ModifiedFiles;
            untracked_files = git_repos{iFold}.UntrackedFiles;
            if isempty(mod_files) && isempty(untracked_files)
                stat_tab.modified(iFold) = false;
            else
                stat_tab.modified(iFold) = true;

                % write modification, try to make format consistant with git status --porcelain.
                % Note that this is not really the same, more work is needed to match.
                mod_files = extractAfter(mod_files,git_dirs(iFold) + filesep);
                untracked_files = extractAfter(untracked_files,git_dirs(iFold) + filesep);
                stat_tab.modified_res(iFold) = join(...
                    [" M " + mod_files; ...
                    "?? " + untracked_files], ...
                    newline);
                stat_tab.modified_res(iFold) = stat_tab.modified_res(iFold) + newline;
            end

        else
            % git is working! faster solution, and more informative in
            % modified_res.

            % save last commit info
            stat_tab.last_commit(iFold) = result(1:7);

            % collect if repository is modified
            [~, result] = system('git status --porcelain');
            % if porcelain result is empty, no change was made
            if isempty(result)
                stat_tab.modified(iFold) = false;
            else
                stat_tab.modified(iFold) = true;
                stat_tab.modified_res(iFold) = string(result);
            end
        end

        % write down repos data if requested
        if nargout == 2 && ( ~exist('git_repos','var') || ~isa(git_repos{iFold},'matlab.git.GitRepository') )
            % get repo object if it was requested & not collected yet
            git_repos{iFold} = gitrepo(git_dirs(iFold)); % get repo object
        end

    end

    if nargout == 2
        git_repos = vertcat(git_repos{:});
    elseif nargout == 0
        disp(stat_tab)
        clearvars stat_tab git_repos
    end
end