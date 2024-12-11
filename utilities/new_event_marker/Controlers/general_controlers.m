classdef general_controlers < handle

    properties
        app_base
        model data2disp_model = data2disp_general_model()
    end
    
    properties
        save_button
        save_loc_menu
        add_data_menu
        export_data_menu
    end
    
    properties
        save_loc
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% Setup, Adding Componants %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods
        function app = general_controlers(app_base, model)
            app.app_base = app_base;
            app.model = model;
        end

        function set_controlers(app, options)
            arguments
                app
                options.save_button
                options.save_loc_menu
                options.add_data_menu
                options.export_data_menu
                options.add_close_req
            end

            if isfield(options, "save_button")
                app.create_save_button(options.save_button{:});
            end

            if isfield(options, "save_loc_menu")
                app.create_save_loc_menu(options.save_loc_menu);
            end

            if isfield(options, "add_data_menu")
                app.create_add_data_menu(options.add_data_menu);
            end

            if isfield(options, "export_data_menu")
                app.create_export_data_menu(options.export_data_menu);
            end

            if isfield(options, "add_close_req") && options.add_close_req
                app.app_base.UIFigure.CloseRequestFcn = @(~,~) app.UIFigure_close_req;
            end
        end

    end

    methods (Hidden)
        function create_save_button(app, button_p, button_row, button_col)
            app.save_button = uibutton(button_p, 'push');
            app.save_button.ButtonPushedFcn = @(~,~) app.save_button_pushed;
            app.save_button.Icon = fullfile( matlabroot, "toolbox", "matlab", "icons", "file_save.png");
            app.save_button.IconAlignment = 'center';
            app.save_button.Text = '';
            app.save_button.Layout.Row = button_row;
            app.save_button.Layout.Column = button_col;
        end

        function create_save_loc_menu(app, menu_p)
            app.save_loc_menu = uimenu(menu_p);
            app.save_loc_menu.Text = "&Change Save Loc";
            app.save_loc_menu.Tooltip = "Change app saving locaction";
            app.save_loc_menu.MenuSelectedFcn = @(~,~) app.save_loc_clicked;
        end

        function create_add_data_menu(app, menu_p)
            app.add_data_menu = uimenu(menu_p);
            app.add_data_menu.Text = "&Add Data";
            app.add_data_menu.Tooltip = 'Add new dataset into display';
            app.add_data_menu.MenuSelectedFcn = @(~,~) app.add_data_clicked;
        end

        function create_export_data_menu(app, menu_p)
            app.export_data_menu = uimenu(menu_p);
            app.export_data_menu.Text = "&Export Data";
            app.export_data_menu.Tooltip = "Export events data to a file";
            app.export_data_menu.MenuSelectedFcn = @(~,~) app.export_data_clicked;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% Componants Responce %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods
        % Button pushed function: save_button
        function save_button_pushed(app)
            % save the app for loading later
            
            % find where to save
            if isempty(app.save_loc)
                [file_name, file_path] = uiputfile('*.mat','Choose where to save',...
                    'event_marker_' + string(datetime("now","Format","uuuuMMdd_HHmmss")));
                if isnumeric(file_name)
                    % user aborted
                    return
                else
                    app.save_loc = fullfile(file_path,file_name);
                end
            end
            
            progress_dlg = uiprogressdlg(app.app_base.UIFigure,...
                "Indeterminate","on","Message","Saving App","Cancelable","off","Title","Saving...");
            cl = onCleanup(@() close(progress_dlg)); % remove progress_dlg block on any form of exit (i.e. error)

            % save the app
            s.app = app.app_base;
            save(app.save_loc,'-struct',"s")
        end

        % Menu selected function:add_data_menu
        function add_data_clicked(app)
            % load variables from workspace & files, add the to app

            %prevent interaction while working
            progress_dlg = uiprogressdlg(app.app_base.UIFigure,"Indeterminate","on","Message","Building Display","Cancelable","off","Title","Building Display");
            cl1 = onCleanup(@() close(progress_dlg)); % remove progress_dlg block on any form of exit (i.e. error)

            app.model.add_data()
        end

        % Menu selected function: export_data_menu
        function export_data_clicked(app)
            % let user export events data to file
            
            % find where to save - note that this is only loosly attached
            % to app.save_loc, and lways asking. This is to let user
            % seperate data & app.
            [file_name, file_path] = uiputfile('*.mat','Choose where to save',app.save_loc);
            if isnumeric(file_name)
                % user aborted
                return
            else
                save_at = fullfile(file_path,file_name);
            end

            % show thinking
            progress_dlg = uiprogressdlg(app.app_base.UIFigure,...
                "Indeterminate","on","Message","Saving Data","Cancelable","off","Title","Saving...");
            cl = onCleanup(@() close(progress_dlg)); % remove progress_dlg block on any form of exit (i.e. error)

            app.model.save_only_evt_data(save_at)
        end

        % Menu Selected function: save_loc_menu
        function save_loc_clicked(app)
            % let user define where they want to save app, when saving.

            % create default for uiputfile - use current save path if exist
            if isempty(app.save_loc)
                def_name =  'event_marker_' + string(datetime("now","Format","uuuuMMdd_HHmmss"));
            else
                def_name = app.save_loc;
            end
            
            % open uiput file, save new path
            [file_name, file_path] = uiputfile('*.mat','Choose where to save',def_name);
            if isnumeric(file_name)
                % user aborted, do not change save loc
                app.save_loc = app.save_loc;
            else
                app.save_loc = fullfile(file_path,file_name);
            end

            % notify user about change
            if isempty(app.save_loc)
                msg = "Current save path is unset.";
            else
                msg = "Current save path is:" + newline + app.save_loc;
            end
            uialert(app.app_base.UIFigure,msg,'Save Path Change','Icon','info')
        end

        % Close request function: app.app_base.UIFigure
        function UIFigure_close_req(app)
            % ask to save before exit
            
            %%%%%%%%%% ask user stuff
            % ask if to save
            answer = uiconfirm(app.app_base.UIFigure,'Save app before closing?','Save On Quit',...
                'Options',["Yes","Only Data","No","Cancle"],'CancelOption','Cancle');
            if answer == "Cancle"
                % this give user the chance to regret
                return
            end
            
            % find where to save if needed
            if answer ~= "No" && isempty(app.save_loc)
                [file_name, file_path] = uiputfile('*.mat','Choose where to save',...
                    'event_marker_' + string(datetime("now","Format","uuuuMMdd_HHmmss")));
                if isnumeric(file_name)
                    % user aborted
                    return
                else
                    app.save_loc = fullfile(file_path,file_name);
                end
            end

            %%%%%%%%%% save as requested
            if answer ~= "No"
                % show thinking
                progress_dlg = uiprogressdlg(app.app_base.UIFigure,...
                    "Indeterminate","on","Message","Saving","Cancelable","off","Title","Saving...");
                cl = onCleanup(@() close(progress_dlg)); % remove progress_dlg block on any form of exit (i.e. error)

                if answer == "Yes"
                    % save the whole app
                    progress_dlg.Message = "Saving App";
                    s.app = app.app_base;
                    save(app.save_loc,'-struct',"s")

                elseif answer == "Only Data"
                    % save only model data
                    progress_dlg.Message = "Saving Data";
                    app.model.save_only_evt_data(app.save_loc)

                end
            end

            % delete graphics
            delete(app.app_base.UIFigure)
            app.app_base.delete_graphics_componants()
        end
    
    end

end