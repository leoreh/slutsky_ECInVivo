classdef IED_general_controlers < general_controlers

    methods
        function save_button_pushed(app)
            % instead of saving the app, save only the ied object
            % block app while saving
            progress_dlg = uiprogressdlg(app.app_base.UIFigure,...
                "Indeterminate","on","Message","Saving IED","Cancelable","off","Title","ied saving");
            cl = onCleanup(@() close(progress_dlg)); % remove progress_dlg block on any form of exit (i.e. error)
            
            ied = app.model.ied;

            % find where to save if needed
            if isempty(ied.file_loc) || ~isfolder(fileparts(ied.file_loc))
                progress_dlg.Message = "Choose where to save";
                [save_file,save_path] = uiputfile('*.ied.mat');
                if isnumeric(save_file)
                    % user cancled, return
                    return
                else
                    % place location in ied
                    ied.file_loc = fullfile(save_path,save_file);
                end
            end

            progress_dlg.Message = "Saving (might take a while)";
            save(ied.file_loc,"ied","-v7.3")
        end

        function UIFigure_close_req(app)
            % raise exit validation
            
            answer = uiconfirm(app.app_base.UIFigure,...
                {'Ready to leave?' 'Changes are already updated in the IED.data obj'}, 'Exit Curation',...
                "Options",["Yes","No"],...
                'CancelOption',"No","DefaultOption","No");
            if answer == "Yes"
                % delete graphics
                delete(app.app_base.UIFigure)
                app.app_base.delete_graphics_componants()
            end
            
        end
    end

end