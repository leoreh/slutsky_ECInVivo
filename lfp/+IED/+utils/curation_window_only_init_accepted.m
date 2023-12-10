classdef curation_window_only_init_accepted < IED.utils.curation_window
    % same as curation_window, but only run on the values that were accepted
    % initaly
    properties
        pos2use (1,:) double {mustBeVector,mustBeNumeric}
    end

    % constructor
    methods
        function app = curation_window_only_init_accepted(ied)

            % collect positions marked as accepted
            pos2use = find(ied.accepted);
            if ~ismember(ied.last_mark,pos2use)
                ied.last_mark = pos2use(1);
            end

            % build app
            app@IED.utils.curation_window(ied)
            app.pos2use = pos2use;
        end
    end

    % override superclass
    methods

        % Window key press function: IED_curation_UIFigure
        function IED_curation_UIFigureWindowKeyPress(app, event)
            % Manage key response - what each key means.
            % Note this manage both keyboard & mouse clicks, by being
            % called from IED_curation_UIFigureWindowButtonDown.
            
            % find event location in pos2use
            [~,pos2use_idx_now] = ismember(app.ied.last_mark, app.pos2use);

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
                    if numel(app.pos2use) < pos2use_idx_now+1
                        close(app.IED_curation_UIFigure)
                        return
                    else
                        app.ied.last_mark = app.pos2use(pos2use_idx_now+1);
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
                    if numel(app.pos2use) < pos2use_idx_now+1
                        close(app.IED_curation_UIFigure)
                        return
                    else
                        app.ied.last_mark = app.pos2use(pos2use_idx_now+1);
                    end

                case {'3','extend','numpad3'} % extend = middle click
                    % return to previous discharge

                    % validate that it is not the first one
                    if 1 == pos2use_idx_now
                        uialert(app.IED_curation_UIFigure,...
                            'Can''t go back from first IED!','invalid return request')
                    else
                        app.ied.last_mark = app.pos2use(pos2use_idx_now-1);
                    end

                otherwise
                    % any other key, just skip it
                    return
            end

            % update graphs
            app.change_discharge()
        end

    end

end