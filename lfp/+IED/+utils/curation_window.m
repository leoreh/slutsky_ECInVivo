classdef curation_window < handle
% open a window to curate IED.data - 
% user decide for each discharged if it is good or not.
% Note that since IED.data is a handle class, curation_window update it
% directly.
%
%   INPUT:
%       ied - IED.data obj, data with detected discharges you wish to
%             collect.
%   OUTPUT:
%       app - curation_window obj, graphic manager to govern the app
%             behavior, NOT deleted after uifigure closed.
%
% Based on getIIS by LH (see +IED/legacy folder)
% By: LdM 
% Published: 230827
%
%   see also IED.data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%% Properties that correspond to app components %%%%%%%%%%%%%%
    properties (Access = public)
        IED_curation_UIFigure  matlab.ui.Figure
        GridLayout             matlab.ui.container.GridLayout
        title_instruction      matlab.ui.control.Label
        zoom_in_axe            matlab.ui.control.UIAxes
        zoom_out_axe           matlab.ui.control.UIAxes
        zoom_in_axe_emg        matlab.ui.control.UIAxes
        save_button            matlab.ui.control.Button
    end

    
    %%%%%%%%%%%%%% graphic properties %%%%%%%%%%%%%%%%
    properties (Access = public)
        %%% zoom out props %%%
        zoom_out_line     % line in the zoom out graph
        zoom_out_peak     % marker showing current discharge
        zoom_out_thr      % lines showing threshold values
        zoom_out_marg = 5 % zoom out margins (one side), in seconds

        %%% zoom in props %%%
        zoom_in_line       % line in the zoom in graph
        zoom_in_peak       % marker showing current discharge
        zoom_in_thr        % lines showing threshold values
        zoom_in_marg = 0.5 % zoom in margins (one side), in seconds
        
        %%% zoom in emg props %%%
        zoom_in_emg_line       % line in the zoom in graph
        zoom_in_emg_peak       % marker showing current discharge

        %%% general %%%
        block_keypress = false % when true, window-press give precedence to ax toolbars
    end
    
    %%%%%%%%%%%%%%  data properties  %%%%%%%%%%%%%%%%
    properties (Access = public)
%         pos
%         sig
%         fs
%         accepted
%         thr
%         thrDir
        ied         % handle to IED.data to curate
        tstamps     % time stamps to match ied.sig
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    %%%%%%%%%%%%%%  app management %%%%%%%%%%%%%%%%
    methods (Hidden)
        function [x_data, y_data] = collect_data(app,nDischarge,marg,src)
            % collect what x_data & y_data are relevant to discharge.
            % INPUT:
            %   app - base obj.
            %   nDischarge - numeric scalar, idx of discharge in detection
            %   marg       - how much to collect around the discharge [in sec]
            %   src        - text scalar, from where in the ied to collect the data ("sig" / "emg")
            % OUTPUT:
            %   x_data, y_data - numeric vectors, data to display.
            
            % get pos to work on
            Discharge_pos = app.ied.pos(nDischarge);

            % convert marg to samples
            marg = app.ied.fs * marg;

            % collect window start & end
            win_start = max([1 Discharge_pos - marg]);
            win_end   =  min([length(app.ied.(src)) Discharge_pos + marg]);

            % transform to x & y vectors
            x_data = app.tstamps(win_start : win_end);
            y_data = app.ied.(src)(win_start : win_end);
        end

        function change_discharge(app)
            % change the display to show the requested discharge
            discharge_num = app.ied.last_mark;

            % collect discharge info
            [discharge_loc_x, discharge_loc_y] = app.collect_data(discharge_num, 0, "sig");

            % update zoom out
            [app.zoom_out_line.XData, app.zoom_out_line.YData] = app.collect_data(discharge_num, app.zoom_out_marg, "sig");
            [app.zoom_out_peak.XData, app.zoom_out_peak.YData] = deal(discharge_loc_x,discharge_loc_y);

            % update zoom in
            [app.zoom_in_line.XData, app.zoom_in_line.YData] = app.collect_data(discharge_num, app.zoom_in_marg, "sig");
            [app.zoom_in_peak.XData, app.zoom_in_peak.YData] = deal(discharge_loc_x,discharge_loc_y);
            title(app.zoom_in_axe, sprintf('spk %d/%d', app.ied.last_mark, numel(app.ied.pos)) )
            ylim(app.zoom_in_axe, [min(app.zoom_in_line.YData) max(app.zoom_in_line.YData)])

            % update zoom in emg
            if isvalid(app.zoom_in_axe_emg)
                [app.zoom_in_emg_line.XData, app.zoom_in_emg_line.YData] = app.collect_data(discharge_num, app.zoom_in_marg, "emg");
                [discharge_loc_x, discharge_loc_y] = app.collect_data(discharge_num, 0, "emg");
                [app.zoom_in_emg_peak.XData, app.zoom_in_emg_peak.YData] = deal(discharge_loc_x,discharge_loc_y);
            end
        end
        
    end
    
    %%%%%%%%%%%%%%  Callbacks that handle component events %%%%%%%%%%%%%%%%
    methods (Hidden)

        % Close request function: IED_curation_UIFigure
        function IED_curation_UIFigureCloseRequest(app)
            % raise exit validation
            
            answer = uiconfirm(app.IED_curation_UIFigure,...
                {'Ready to leave?' 'Changes are already updated in the IED.data obj'}, 'Exit Curation',...
                "Options",["Yes","No"],...
                'CancelOption',"No","DefaultOption","No");
            if answer == "Yes"
                delete(app.IED_curation_UIFigure)
            end
            
        end

        % Window key press function: IED_curation_UIFigure
        function IED_curation_UIFigureWindowKeyPress(app, event)
            % Manage key response - what each key means.
            % Note this manage both keyboard & mouse clicks, by being
            % called from IED_curation_UIFigureWindowButtonDown. Normal
            % callback only detect keyboard clicks, otherwise.
            
            % validate that press wasn't on save_button
            % keyboard
            
            key = event.Key;
            switch key
                case {'1','normal','numpad1'} % normal = left click
                    % accept discharge

                    % flash green
                    app.zoom_in_axe.Title.Color  = 'g';
                    pause(0.05)
                    app.zoom_in_axe.Title.Color  = 'k';

                    % mark accepted (unnecessary, as this is default)
                    app.ied.accepted(app.ied.last_mark) = true;

                    % move to next discharge
                    if numel(app.ied.pos) < app.ied.last_mark+1
                        close(app.IED_curation_UIFigure)
                        return
                    else
                        app.ied.last_mark = app.ied.last_mark+1;
                    end
                case {'2','return','alt','open','numpad2'} % alt = right click & alt key, open = double click
                    % decline discharge

                    % flash red
                    app.zoom_in_axe.Title.Color  = 'r';
                    pause(0.05)
                    app.zoom_in_axe.Title.Color  = 'k';

                    % mark denied
                    app.ied.accepted(app.ied.last_mark) = false;

                    % move to next discharge
                    if numel(app.ied.pos) < app.ied.last_mark+1
                        close(app.IED_curation_UIFigure)
                        return
                    else
                        app.ied.last_mark = app.ied.last_mark+1;
                    end
                    
                case {'3','extend','numpad3'} % extend = middle click
                    % return to previous discharge

                    % validate that it is not the first one
                    if 1 == app.ied.last_mark
                        uialert(app.IED_curation_UIFigure,...
                            'Can''t go back from first IED!','invalid return request')
                    else
                        app.ied.last_mark = app.ied.last_mark-1;
                    end
                
                otherwise
                    % any other key, just skip it
                    return
            end
            
            % update graphs
            app.change_discharge()
        end

        % Window button up function: IED_curation_UIFigure
        function IED_curation_UIFigureWindowButtonUp(app)
            % run when mouse key is pressed.
            % make sure axtoolbar wasn't released, and pass mouse key to
            % key-pressed function
            pause(0.01)
            if app.block_keypress
                % Clicked object was released.
                % Remove block, but do not move discharge
                app.block_keypress = false;
            else
                % collect what mouse key was pressed, pass it to
                % keypress-manager
                key_pass.Key = app.IED_curation_UIFigure.SelectionType;
                app.IED_curation_UIFigureWindowKeyPress(key_pass)
            end
        end

        % Selection changed function: zoom_in_axe's & zoom_out_axe's axtoolbar
        function axtoolbar_pressed(app)
            % make sure mouse click behavior won't happen when user try to
            % click axtoolbar.
            % Changing the property here is evaluated before
            % IED_curation_UIFigureWindowButtonUp, and inside the next
            % function app.toolbar_pressed block behavior.
            app.block_keypress = true;
        end
        
        % Button pushed function: zoom_in_axe's & zoom_out_axe's toolbarpushbutton
        function pushbutton_axtoolbar_pressed(app,src,evt)
            % axtoolbar SelectionChangedFcn do not respond to pushbuttons.
            % this wrapper make sure to conserve buttons functionality,
            % while still preventing response to user calls.
            app.axtoolbar_pressed()
            switch src.Tag
                case 'restoreview'
                    matlab.graphics.controls.internal.resetHelper(evt.Axes,true)
                case 'Save As'
                    matlab.graphics.internal.export.exportCallback(matlab.graphics.controls.internal.exportHelper(evt.Axes))
                case 'Copy as Image'
                    matlab.graphics.internal.export.exportTo(matlab.graphics.controls.internal.exportHelper(evt.Axes),'target','clipboard','format','image')
                case 'Copy as Vector Graphic'
                    matlab.graphics.internal.export.exportTo(matlab.graphics.controls.internal.exportHelper(evt.Axes),'target','clipboard','format','vector','vector',true)
            end
            
        end
   
        % Button pushed function: save_button
        function save_ied(app)
            % save ied object at its current state
            
            % block IED_curation_UIFigureWindowKeyPress for this press.
            % this will be restore in IED_curation_UIFigureWindowButtonUp,
            % following mouse release.
            app.block_keypress = true;

            % block app while saving
            progress_dlg = uiprogressdlg(app.IED_curation_UIFigure,...
                "Indeterminate","on","Message","Saving IED","Cancelable","off","Title","ied saving");
            cl = onCleanup(@() close(progress_dlg)); % remove progress_dlg block on any form of exit (i.e. error)
            
            % find where to save if needed
            if isempty(app.ied.file_loc) || ~isfolder(fileparts(app.ied.file_loc))
                progress_dlg.Message = "Choose where to save";
                [save_file,save_path] = uiputfile('*.ied.mat');
                if isnumeric(save_file)
                    % user cancled, return
                    return
                else
                    % place location in ied
                    app.ied.file_loc = fullfile(save_path,save_file);
                end
            end

            progress_dlg.Message = "Saving (might take a while)";
            ied = app.ied; %#ok<PROP> this is so saved name will be ied as expected
            save(ied.file_loc,"ied","-v7.3") %#ok<PROP> this is so saved name will be ied as expected
        end
    end

    
    %%%%%%%%%%%%%%  Component initialization %%%%%%%%%%%%%%%%
    methods (Hidden)

        % Create UIFigure and components
        function createComponents(app)
            
            % collect matlab relese
            mat_ver = ver('matlab');
            mat_ver = mat_ver.Release;
            
            % Create IED_curation_UIFigure and hide until all components are created
            app.IED_curation_UIFigure = uifigure('Visible', 'off');
            app.IED_curation_UIFigure.Color = [1 1 1];
            app.IED_curation_UIFigure.Position = [100 100 640 480];
            app.IED_curation_UIFigure.Name = 'IED Curation';
            app.IED_curation_UIFigure.CloseRequestFcn = @(~,~) app.IED_curation_UIFigureCloseRequest;
            app.IED_curation_UIFigure.WindowButtonUpFcn =  @(~,~) app.IED_curation_UIFigureWindowButtonUp;
            app.IED_curation_UIFigure.WindowKeyPressFcn = @(~,evt) app.IED_curation_UIFigureWindowKeyPress(evt);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.IED_curation_UIFigure);
            app.GridLayout.ColumnWidth = {'0.2x','1x'};
            app.GridLayout.RowHeight = {'0.2x', '1x', '1x', '1x'};
            if (str2double(mat_ver(3:6)) > 2020) || ( (str2double(mat_ver(3:6))==2020) && (mat_ver(7) == "b") )
                % only from 2020b and later
                app.GridLayout.BackgroundColor = [1 1 1];
            end

            % Create zoom_out_axe
            app.zoom_out_axe = uiaxes(app.GridLayout);
            xlabel(app.zoom_out_axe, 'Time [s]')
            ylabel(app.zoom_out_axe, 'Voltage [mV]')
            zlabel(app.zoom_out_axe, 'Z')
            app.zoom_out_axe.YLim = [0 1];
            if str2double(mat_ver(3:6)) >= 2021
                % only from 2021a and later
                app.zoom_out_axe.XLimitMethod = 'tight';
            else
                axis(app.zoom_out_axe,'tight')
            end
            app.zoom_out_axe.TickLength = [0 0];
            app.zoom_out_axe.NextPlot = 'add';
            app.zoom_out_axe.Layout.Row = 2;
            app.zoom_out_axe.Layout.Column = [1 2];
            app.zoom_out_axe.XAxis.Exponent = 0;
            
            % manage zoom_out_axe toolbar
            buttons2disp = ["export","brush","pan","zoomin","zoomout","restoreview"];
            bt_strip = axtoolbar(app.zoom_out_axe,buttons2disp,"SelectionChangedFcn",@(~,~) app.axtoolbar_pressed);
            % add custom data-tips button - it refused to be added to uiaxe otherwise
            axtoolbarbtn(bt_strip,"state","Icon","datacursor","Tooltip","Data Tips",...
                "ValueChangedFcn",@(e,d)datacursormode(ancestor(d.Source,'figure'),d.Value));
            % make sure push buttons are also clickable without moving discharged
            push_buttons = findobj(bt_strip,'Type','toolbarpushbutton');
            [push_buttons.ButtonPushedFcn] = deal(@(src,evt) app.pushbutton_axtoolbar_pressed(src,evt));

            % Create zoom_in_axe
            app.zoom_in_axe = uiaxes(app.GridLayout);
            ylabel(app.zoom_in_axe, 'Voltage [mV]')
            zlabel(app.zoom_in_axe, 'Z')
            app.zoom_in_axe.YLim = [0 1];
            if str2double(mat_ver(3:6)) >= 2021
                % only from 2021a and later
                app.zoom_in_axe.XLimitMethod = 'tight';
            else
                axis(app.zoom_in_axe,'tight')
            end
            app.zoom_in_axe.TickLength = [0 0];
            app.zoom_in_axe.NextPlot = 'add';
            app.zoom_in_axe.Layout.Row = 3;
            app.zoom_in_axe.Layout.Column = [1 2];
            app.zoom_in_axe.XAxis.Exponent = 0;

            % manage zoom_in_axe toolbar
            bt_strip = axtoolbar(app.zoom_in_axe,buttons2disp,"SelectionChangedFcn",@(~,~) app.axtoolbar_pressed);
            % add custom data-tips button - it refused to be added to uiaxe otherwise
            axtoolbarbtn(bt_strip,"state","Icon","datacursor","Tooltip","Data Tips",...
                "ValueChangedFcn",@(e,d)datacursormode(ancestor(d.Source,'figure'),d.Value));
             % make sure push buttons are also clickable without moving discharged
            push_buttons = findobj(bt_strip,'Type','toolbarpushbutton');
            [push_buttons.ButtonPushedFcn] = deal(@(src,evt) app.pushbutton_axtoolbar_pressed(src,evt));
            
            % Create zoom_in_axe_emg
            app.zoom_in_axe_emg = uiaxes(app.GridLayout);
            ylabel(app.zoom_in_axe_emg, 'Voltage [mV]')
            zlabel(app.zoom_in_axe_emg, 'Z')
            app.zoom_in_axe_emg.YLim = [0 1];
            if str2double(mat_ver(3:6)) >= 2021
                % only from 2021a and later
                app.zoom_in_axe_emg.XLimitMethod = 'tight';
            else
                axis(app.zoom_in_axe_emg,'tight')
            end
            app.zoom_in_axe_emg.TickLength = [0 0];
            app.zoom_in_axe_emg.NextPlot = 'add';
            app.zoom_in_axe_emg.Layout.Row = 4;
            app.zoom_in_axe_emg.Layout.Column = [1 2];
            app.zoom_in_axe_emg.XAxis.Exponent = 0;
            title(app.zoom_in_axe_emg, 'EMG zoomed in')

            % manage zoom_in_axe_emg toolbar
            bt_strip = axtoolbar(app.zoom_in_axe_emg,buttons2disp,"SelectionChangedFcn",@(~,~) app.axtoolbar_pressed);
            % add custom data-tips button - it refused to be added to uiaxe otherwise
            axtoolbarbtn(bt_strip,"state","Icon","datacursor","Tooltip","Data Tips",...
                "ValueChangedFcn",@(e,d)datacursormode(ancestor(d.Source,'figure'),d.Value));
             % make sure push buttons are also clickable without moving discharged
            push_buttons = findobj(bt_strip,'Type','toolbarpushbutton');
            [push_buttons.ButtonPushedFcn] = deal(@(src,evt) app.pushbutton_axtoolbar_pressed(src,evt));


            % Create title_instruction
            app.title_instruction = uilabel(app.GridLayout);
            app.title_instruction.HorizontalAlignment = 'center';
            app.title_instruction.FontName = 'Arial';
            app.title_instruction.FontSize = 14;
            app.title_instruction.FontWeight = 'bold';
            app.title_instruction.Layout.Row = 1;
            app.title_instruction.Layout.Column = [1 2];
            app.title_instruction.Text = {'Inspect IIS:'; 'left/1 = accept; right/2/enter/alt = decline; middle/3 = previous'};

            % create save button
            app.save_button = uibutton(app.GridLayout,"push");
            app.save_button.Text = "Save ied object";
            app.save_button.WordWrap = "on";
            app.save_button.HorizontalAlignment = "center";
            app.save_button.ButtonPushedFcn = @(~,~) app.save_ied();
            app.save_button.Layout.Row = 1;
            app.save_button.Layout.Column = 1;

            % Show the figure after all components are created
            app.IED_curation_UIFigure.Visible = 'on';
        end
    end

    %%%%%%%%%%%%%%  App creation and deletion %%%%%%%%%%%%%%%%
    methods (Access = public)

        % Construct app
        function app = curation_window(ied)
            % place graphic where it should be, insert properites into the
            % object

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % arguments
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if ~isa(ied,'IED.data')
                error("first input must be IED.data obj")
            end

            % point to ied from app
            app.ied = ied;
            
            % params
            app.tstamps = ( 1:length(app.ied.sig) )' / app.ied.fs;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % init graphics
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Create UIFigure and components
            createComponents(app)

            % collect discharge info
            [discharge_loc_x, discharge_loc_y] = app.collect_data(app.ied.last_mark, 0, "sig");
            
            % prep zoom out graph
            [x_data, y_data] = app.collect_data(app.ied.last_mark, app.zoom_out_marg, "sig");
            app.zoom_out_line = plot(app.zoom_out_axe ,x_data, y_data);
            app.zoom_out_peak = scatter(app.zoom_out_axe, discharge_loc_x, discharge_loc_y, 'r*');
            sig_mean = mean(app.ied.sig);
            sig_std = std(app.ied.sig);
            diplay_limit = 1.5*app.ied.thr(1)*sig_std;
            ylim(app.zoom_out_axe,[sig_mean-diplay_limit, sig_mean+diplay_limit])
            
            % prep zoom in graph
            [x_data, y_data] = app.collect_data(app.ied.last_mark, app.zoom_in_marg, "sig");
            app.zoom_in_line = plot(app.zoom_in_axe ,x_data, y_data);
            app.zoom_in_peak = scatter(app.zoom_in_axe, discharge_loc_x, discharge_loc_y, 'r*');
            ylim(app.zoom_in_axe, [min(y_data) max(y_data)])
            title(app.zoom_in_axe, sprintf('spk %d/%d', app.ied.last_mark, numel(app.ied.pos)));

            % prep zoom in emg graph
            if ~isempty(app.ied.emg)
                % ied have emg data, use it
                
                % create lines
                [x_data, y_data] = app.collect_data(app.ied.last_mark, app.zoom_in_marg, "emg");
                [discharge_loc_x, discharge_loc_y] = app.collect_data(app.ied.last_mark, 0, "emg");
                app.zoom_in_emg_line = plot(app.zoom_in_axe_emg ,x_data, y_data);
                app.zoom_in_emg_peak = scatter(app.zoom_in_axe_emg, discharge_loc_x, discharge_loc_y, 'r*');

                % choose ylimit - same as zoom_out
                ylim(app.zoom_in_axe_emg,[sig_mean-diplay_limit, sig_mean+diplay_limit])
            else
                % ied do not have emg data, remove unused axe
                row2rmv = app.zoom_in_axe_emg.Layout.Row;
                delete(app.zoom_in_axe_emg)
                app.GridLayout.RowHeight(row2rmv) = [];
            end

            % add threshold markers
            if ismember(app.ied.thrDir,["positive","both"])
                app.zoom_out_thr(1) = yline(app.zoom_out_axe, app.ied.thr(2), '--r');
                app.zoom_in_thr(1)  = yline(app.zoom_in_axe,  app.ied.thr(2), '--r');
            end
            if ismember(app.ied.thrDir,["negative","both"])
                app.zoom_out_thr(2) = yline(app.zoom_out_axe, -app.ied.thr(2), '--r');
                app.zoom_in_thr(2)  = yline(app.zoom_in_axe,  -app.ied.thr(2), '--r');
            end

            
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.IED_curation_UIFigure)
        end
    end
end

% EOF