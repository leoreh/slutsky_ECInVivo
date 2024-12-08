function path_shadows_solver(name2solve)

    all_options = what(name2solve);
    all_options = {all_options.path};
    if numel(all_options) > 1
        [indx,tf] = listdlg('ListString',all_options,...
            'PromptString',["This package has multiple optional on your path." "Choose path to Keep:" ""],...
            "CancelString","Abort",...
            "ListSize",[500 250],"SelectionMode","single","Name","Path Conflict Problem");
        if tf == 0
            error("path_shadows_solver:Abourt",'Unsolved Path Conflict.')
        else
            all_options(indx) = [];
            folds2rmv = fileparts(all_options);
            for iFold = 1:numel(folds2rmv)
                warning('off', 'MATLAB:rmpath:DirNotFound');
                rmpath(genpath(folds2rmv{iFold}))
                warning('on', 'MATLAB:rmpath:DirNotFound');
            end
        end
    end

end