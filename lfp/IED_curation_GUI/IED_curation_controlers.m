classdef IED_curation_controlers < handle
    properties
        app_base
        model IED_data2disp
    end

    properties
        decline_button matlab.ui.control.Button
        accept_button  matlab.ui.control.Button
        return_button  matlab.ui.control.Button
    end
    
    
    methods
        function app = IED_curation_controlers(app_base, model)
            app.app_base = app_base;
            app.model = model;
        end

        function set_controlers(app, options)
            arguments
                app

                options.decline_button
                options.accept_button
                options.return_button
                options.jump2detection
            end

            if isfield(options, "decline_button")
                app.create_decline_button(options.decline_button{:})
            end

            if isfield(options, "accept_button")
                app.create_accept_button(options.accept_button{:})
            end

            if isfield(options, "return_button")
                app.create_return_button(options.return_button{:})
            end

            if isfield(options, "jump2detection")
                app.create_jump2detection_button(options.jump2detection{:})
            end
        end
    
    end

    methods (Hidden)
        function create_decline_button(app, button_parent, button_row, button_col)
            button = matlab.ui.control.Button();
            button.ButtonPushedFcn = @(~,~) app.decline_button_pushed;
            button.Tooltip = {'Decline IED'};
            button.Text = 'Decline';
            button.Parent = button_parent;
            button.Layout.Row = button_row;
            button.Layout.Column = button_col;
            app.decline_button = button;
        end

        function create_accept_button(app, button_parent, button_row, button_col)
            button = matlab.ui.control.Button();
            button.ButtonPushedFcn = @(~,~) app.accept_button_pushed;
            button.Tooltip = {'Accept IED'};
            button.Text = 'Accept';
            button.Parent = button_parent;
            button.Layout.Row = button_row;
            button.Layout.Column = button_col;
            app.accept_button = button;
        end
        
        function create_return_button(app, button_parent, button_row, button_col)
            button = matlab.ui.control.Button();
            button.ButtonPushedFcn = @(~,~) app.return_button_pushed;
            button.Tooltip = {'Return to Previous IED'};
            button.Text = 'Return';
            button.Parent = button_parent;
            button.Layout.Row = button_row;
            button.Layout.Column = button_col;
            app.return_button = button;
        end
    
        function create_jump2detection_button(app, button_parent, button_row, button_col)
            button = matlab.ui.control.Button();
            button.ButtonPushedFcn = @(~,~) app.jump2detection_pushed;
            button.Tooltip = {'Jump view to current detection'};
            button.Text = 'Jump 2 detection';
            button.WordWrap = "on";
            button.Parent = button_parent;
            button.Layout.Row = button_row;
            button.Layout.Column = button_col;
            app.return_button = button;
        end
    end


    methods
        function decline_button_pushed(app)
            try
                app.model.change_detection_status(false)
                if ~isempty(app.decline_button) && isvalid(app.decline_button)
                    % flash button red
                    app.decline_button.FontColor = 'r';
                    pause(0.1)
                    app.decline_button.FontColor = 'k';
                end
                app.model.move_IED_detection("Next")
            catch err
                if strcmp(err.identifier,'IED_data2disp:move_IED_detection:out_of_bounds')
                    if strcmp(err.message, 'Cant Move after last detection')
                        close(app.app_base.UIFigure)
                    else
                        uialert(app.app_base.UIFigure, err.message, 'Out of detections','Icon','error', ...
                            'CloseFcn',@(~,~) focus(app.app_base.UIFigure))
                    end
                else
                    rethrow(err)
                end
            end
            % refocus after click - only if app is still there
            if isvalid(app.app_base) && isvalid(app.app_base.UIFigure)
                focus(app.app_base.UIFigure)
            end
        end

        function accept_button_pushed(app)
            try
                app.model.change_detection_status(true)
                if ~isempty(app.accept_button) && isvalid(app.accept_button)
                    % flash button green
                    app.accept_button.FontColor = 'g';
                    pause(0.1)
                    app.accept_button.FontColor = 'k';
                end
                app.model.move_IED_detection("Next")
            catch err
                if strcmp(err.identifier,'IED_data2disp:move_IED_detection:out_of_bounds')
                    if strcmp(err.message, 'Cant Move after last detection')
                        close(app.app_base.UIFigure)
                        return
                    else
                        uialert(app.app_base.UIFigure, err.message, 'Out of detections','Icon','error',...
                            'CloseFcn',@(~,~) focus(app.app_base.UIFigure))
                    end
                else
                    rethrow(err)
                end
            end
            focus(app.app_base.UIFigure)
        end

        function return_button_pushed(app)
            try
                app.model.move_IED_detection("Previous")
            catch err
                if strcmp(err.identifier,'IED_data2disp:move_IED_detection:out_of_bounds')
                    uialert(app.app_base.UIFigure, err.message, 'Out of detections','Icon','error',...
                        'CloseFcn',@(~,~) focus(app.app_base.UIFigure))
                else
                    rethrow(err)
                end
            end
            focus(app.app_base.UIFigure)
        end
        
        function jump2detection_pushed(app)
            ied = app.model.ied;
            detect_time = ied.pos(ied.last_mark)./ied.fs;
            app.model.move_time(detect_time);
            focus(app.app_base.UIFigure)
        end

    end

    methods
        function [help_text] = key_responce(app, key_pressed)
            % create keyboard shortcuts

            % create help text
            if nargout > 0
                help_text = ...
                    "1 -> accept detection" + newline + ...
                    "2 / return -> decline detection" + newline + ...
                    "3 -> retrun to previous";
                return
            end

            % create key responce
            switch key_pressed
                case {'1','numpad1'}
                    app.accept_button_pushed()
                case  {'2','return', 'numpad2'}
                    app.decline_button_pushed()
                case {'3','extend','numpad3'}
                    app.return_button_pushed()
            end
        end
    
    end
end