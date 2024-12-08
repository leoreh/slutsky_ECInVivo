classdef event_map < handle

    properties
        app_base
        model data2disp_model = data2disp_general_model()
    end

    properties (SetObservable=true)
        % hold the data to display.
        % Is a table, with variables "local" & "global", and rows are event
        % types. Variables contain lines.
        % Therefore, indexsable by event_lines.local(event_type).
        event_lines table = table('Size',[0 2],...
            'VariableTypes', ["matlab.graphics.chart.primitive.Line" "matlab.graphics.chart.primitive.Line"], ...
            'VariableNames',["local" "global"]);
    end

    properties
        % parts of map displaying the events only inside current
        % model.events_time2disp (localy)
        local_map_ax matlab.ui.control.UIAxes
        global_map_ax matlab.ui.control.UIAxes
        
        ROI_time images.roi.Line

        % collect legend for all lines
        event_lg    matlab.graphics.illustration.Legend
    end
    
    properties (Dependent)
        local_exists
        global_exists

        all_lines
        lines_types
    end

    properties
        update_local_in_win_change_lst  event.listener
        update_local_event_changed_lst  event.listener
        update_global_in_win_change_lst event.listener
        update_global_event_changed_lst event.listener

        events_types_changed_listener   event.listener

        global_lim_changed_lst          event.listener
    end

    events
        line_pressed
    end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%% Setup, Adding Componants %%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods
        function app = event_map(app_base, model)
            % connect view to model & base
            
            arguments
                app_base
                model data2disp_model
            end

            % save handles to main objects
            app.app_base = app_base;
            app.model = model;
        end

        function set_view(app, options)
            % create where local map & global map
            arguments
                app
                options.local_map_loc
                options.global_map_loc
            end

            % create local map when requested
            if isfield(options, "local_map_loc")
                app.create_local_map(options.local_map_loc{:})
            end

            % create global map & ROI when requested
            if isfield(options, "global_map_loc")
                app.create_global_map(options.global_map_loc{:})
            end
            
            % add listener to model
            app.events_types_changed_listener = listener(app.model, "events_types_changed", ...
                @(~,~) app.remake_evt_lines);
        end
    end

    methods (Hidden)
        function create_global_map(app, map_ax_parent, map_ax_row, map_ax_col)
            % create global_map_ax
            app.global_map_ax = uiaxes(map_ax_parent);
            app.global_map_ax.YLim = [-1.1 1.1];
            app.global_map_ax.XAxisLocation = 'origin';
            app.global_map_ax.XMinorTick = 'on';
            app.global_map_ax.YTick = [];
            app.global_map_ax.TickDir = 'both';
            app.global_map_ax.XAxis.Exponent = 0;
            app.global_map_ax.Layout.Row = map_ax_row;
            app.global_map_ax.Layout.Column = map_ax_col;
            app.global_map_ax.Interactions = dataTipInteraction;
            app.global_map_ax.InteractionOptions.LimitsDimensions = "x";

            % rescale event_map parts so you can see / interact with them
            % in the current XLims
            % app.global_lim_changed_lst = listener(app.global_map_ax, 'XLim', 'PostSet', ...
            %     @(~,evt) app.map_lims_changed());
            app.global_map_ax.XAxis.LimitsChangedFcn =  @(~,~) app.map_lims_changed;

            % create ROI_time
            app.ROI_time = images.roi.Line(app.global_map_ax,'Deletable',false, ...
                'Position',[0 0; 0 0],'DrawingArea', 'unlimited');
            app.ROI_time.addlistener('MovingROI',@(~,evt) app.ROI_moving(evt));
            
            % change events as needed
            app.update_global_event_changed_lst = listener(app.model, "events_changed", @(~,~) app.update_global_event_map());

            % update ROI location based on model
            app.update_global_in_win_change_lst = listener(app.model, "in_win_change", @(~,~) app.move_ROI());

            % add event_legend if one do not exist
            if isempty(app.event_lg) || ~isvalid(app.event_lg)
                app.event_lg = legend(app.global_map_ax);
                app.event_lg.ItemHitFcn = @(~,evt) app.event_ax_legend_on_click(evt);
            end
        end
    
        function create_local_map(app, map_ax_parent, map_ax_row, map_ax_col)
            % create local_map_ax
            app.local_map_ax = uiaxes(map_ax_parent);
            xlabel(app.local_map_ax, 'Time [Sec]')
            app.local_map_ax.Toolbar.Visible = 'off';
            app.local_map_ax.YLim = [-1.1 1.1];
            app.local_map_ax.YTick = [];
            app.local_map_ax.YTickLabel = [];
            app.local_map_ax.Layout.Row = map_ax_row;
            app.local_map_ax.Layout.Column = map_ax_col;
            app.local_map_ax.XAxis.Exponent = 0;
            app.local_map_ax.PositionConstraint = "innerposition";
            app.local_map_ax.InteractionOptions.LimitsDimensions = "x";
            app.local_map_ax.Toolbar.Visible = "on";
            
            % add event_legend if one do not exist
            if isempty(app.event_lg) || ~isvalid(app.event_lg)
                app.event_lg = legend(app.local_map_ax);
                app.event_lg.ItemHitFcn = @(~,evt) app.event_ax_legend_on_click(evt);
            end

            % change view when model changed
            app.update_local_in_win_change_lst = listener(app.model, "in_win_change", @(~,~) app.update_local_map_lines);
            app.update_local_event_changed_lst = listener(app.model, "events_changed", @(~,~) app.update_local_map_lines);
        end
    end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%% Control ROI & Map Behavior %%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   methods

       function ROI_moving(app,evt)
           % Run while ROI is moving.
           % Does:
           %    1) Prevent ROI movment on the Y axis.
           %    2) If ROI reaches axes limits, move to match %% CANT,
           %    ROI LIMITED BY INIT XLIMS

           % restore previous Y loc
           app.ROI_time.Position(:,2) = evt.PreviousPosition(:,2);
           
           % prevent crossing data limits
           current_pos = sort(app.ROI_time.Position(:,1));
           ROI_size = diff(current_pos);
           if current_pos(1) < 0
               current_pos = [0 ROI_size];
           end
           if current_pos(2) > app.model.longest_data
               current_pos = [app.model.longest_data-ROI_size app.model.longest_data];
           end
           app.ROI_time.Position(:,1) = current_pos;

           % check if ROI pos passed Xlims - 
           % note that this will behave unexpectadly if ROI size > XLims
           current_lims = app.global_map_ax.XLim;
           win_size = diff(current_lims);
           win_mov_val = win_size*0.005;
           
           if (current_pos(1) <= current_lims(1)) && ( (current_pos(1)-win_mov_val) >= 0)
               current_lims = current_pos(1)-win_mov_val + [0 win_size];
           end
           if (current_pos(2) >= current_lims(2)) && ( (current_pos(2)+win_mov_val) <= app.model.longest_data)
               current_lims = current_pos(2)+win_mov_val - [win_size 0];
           end
           
           % reset current lims to match
           app.global_map_ax.XLim = current_lims;
       end
   
       function map_lims_changed(app)
           % rescale elements in map_ax so they will be visable in current XLims
           % If ROI_time is too small, prioritize movment over shape change
            
           x_range = diff(app.global_map_ax.XLim);
           
           %%%%%%%%% update ROI_time size
           % collect relevant sizes
           actual_win_size = app.model.events_time2disp;
           
           ROI_base_dur = diff(actual_win_size);

           % compute new ROI display size
           % scale it so it is 0.5% of axe limits, but not less then its actual size
           scale_factor = max(1,(x_range.*0.005)./ROI_base_dur);
           disp_dur = ROI_base_dur.*scale_factor;
           time2add = disp_dur - ROI_base_dur;

           % set new time
           app.ROI_time.Position(:,1) = actual_win_size + [-time2add time2add]/2;

           %%%%%%%%% update ROI_time interactions
           % this was only tested for figure defined in pixcels
           units = app.app_base.UIFigure.Units;


           % find ROI edges location on screen
           pos = app.ROI_time.Position;
           edge_real_pos = axescoord2figurecoord(pos(:,1), pos(:,2), app.global_map_ax);
           line_size = diff(edge_real_pos);

           % by figure units, choose what limit to use
           switch units
               case {'pixels','points'}
                   limit2use = 50;
               case 'normalized'
                   limit2use = 0.001;
               case 'inches'
                   limit2use = 10/72;
               case 'centimeters'
                   limit2use = 0.035306*10;
           end

           % modify line interactions
           if  line_size < limit2use
               % small ROI - only let movment
               app.ROI_time.InteractionsAllowed = "translate";
           else
               % ROI may be big enough to let user select interaction
               app.ROI_time.InteractionsAllowed = "all";
           end

           %%%%%%%%% Update event size
           app.update_global_event_map();
       end
   
       function move_ROI(app)
           new_time = app.model.events_time2disp;
           if ~isempty(app.ROI_time.Position)
               y_pos = app.ROI_time.Position(1,2);
           else
               y_pos = 0;
           end
           app.ROI_time.Position = [new_time(1) y_pos;new_time(2) y_pos];

           app.map_lims_changed()
       end
   end
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%% Respond to event changes %%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   methods
       function update_global_event_map(app)
            % Update event map on map_axe to display existing events
            % Scale each event so it will show on the "minimap"

            % scale events durations, so longest event is 1% of view
            % (but no smaller than actual)
            all_events = app.model.events_time.values;
            biggest_evt_size = max( diff(vertcat(all_events{:}),1,2) );
            x_range = diff(app.global_map_ax.XLim);
            scale_factor = max(1,(x_range.*0.005)./biggest_evt_size); % do not let scaling be less than 1 - event can't be smaller than truth
            
            % collect what events to modify
            event_types = app.model.events_types;
            nTypes = numel(event_types);
            
            % update lines
            Y_vals = rescale(1:nTypes, -0.5, 0.5);
            for iType = 1:nTypes
                now_type = event_types(iType);

                % collect all events of type
                events_t = app.model.events_time{now_type};

                % convert to nan seperated vec, create matching Y values
                events_t_nansep = app.model.events_nan_vec(events_t);
                Y_data = Y_vals(iType).*ones(size(events_t_nansep)); % spread events between -0.5 & 0.5

                % rescale events
                event_dur = diff(events_t,1,2);
                disp_event_dur = event_dur.*scale_factor;
                time2add = disp_event_dur - event_dur;
                X_data = [events_t(:,1)-(time2add/2) events_t(:,2)+(time2add/2)];
                X_data = app.model.events_nan_vec(X_data);

                % plot the events
                [app.event_lines.global(now_type).XData, app.event_lines.global(now_type).YData] = ...
                    deal(X_data, Y_data);

                % keep datatips refer to real event locations
                app.event_lines.global(now_type).DataTipTemplate.DataTipRows(1).Value = events_t_nansep;

            end
       end
        
       function update_local_map_lines(app)
           % update local map when time window is changed

           % collect what events to modify
           events_types = app.model.events_types;
           nTypes = numel(events_types);

           % update lines
           Y_vals = rescale(1:nTypes, -0.5, 0.5);
           for iType = 1:nTypes
               now_type = events_types(iType);

               % collect data, convert to vec, create matching Y value
               event2disp = app.model.events2disp{now_type};
               X_data = app.model.events_nan_vec(event2disp);
               Y_data = Y_vals(iType)*ones(size(X_data)); % spread events between -0.5 & 0.5

               [app.event_lines.local(now_type).XData, app.event_lines.local(now_type).YData] = ...
                   deal(X_data,Y_data);
           end

           % make sure event_ax do not grow beyond its requested limits
           app.local_map_ax.XLim = app.model.events_time2disp;
       end
       
       function remake_evt_lines(app)
           % add & remove lines as changed in model

           % collect data
           model_types = app.model.events_types;
           current_types = app.lines_types;
           
           % find which types were removed
           types_removed = ~ismember(current_types, model_types);
           if any(types_removed)
               app.rmv_map_lines(current_types(types_removed))
           end

           % find which types were added
           types_added = ~ismember(model_types, current_types);
           if any(types_added)
               app.make_map_lines(model_types(types_added))
           end
       end

       function make_map_lines(app, event_types, event_colors)
           % create new lines for a list of event types.

           %%%%%% Input management

           % validate no type already exist
           dup_lines = ismember(event_types, app.lines_types);
           if any(dup_lines)
               error("Line for event types {%s} already exist - delete them before adding", ...
                   join(event_types(dup_lines),', '))
           end

           % use default colors if no colors given
           if ~exist("event_colors","var")
               % collect colors of existing lines
               if app.local_exists
                   all_line = app.event_lines.local;
               elseif app.global_exists
                   all_line = app.event_lines.global;
               end
               
               % include white background
               if isempty(all_line)
                   all_colors = [1 1 1];
               else
                   all_colors = vertcat(all_line.Color);
                   all_colors = [1 1 1; all_colors];
               end

               % create new colors, distingutiable from existing
               nTypes = numel(event_types);
               event_colors = distinguishable_colors(nTypes, all_colors);
           elseif numel(event_types) ~= size(event_colors,1)
               error('event_colors must have a row for every event_type')
           end

           %%%%%% Add new lines
           warning('off','MATLAB:table:RowsAddedExistingVars') % do not warn on new row names
           for iType = 1:numel(event_types)
               current_type = event_types(iType);

               if app.global_exists
                   app.event_lines.global(current_type) = matlab.graphics.chart.primitive.Line(...
                       'Parent',app.global_map_ax,...
                       'Color',event_colors(iType,:),'LineWidth',3,...
                       'DisplayName',current_type);
               end

               if app.local_exists
                   app.event_lines.local(current_type) = matlab.graphics.chart.primitive.Line(...
                       'Parent',app.local_map_ax,...
                       'Color',event_colors(iType,:),'LineWidth',3,...
                       'DisplayName',current_type,...
                       'ButtonDownFcn',@(~,evt) app.notify_line_pressed(evt));
               end
                
               % bind the lines colors together
               type_lines = app.event_lines{current_type,:};
               type_lines = type_lines(isgraphics(type_lines));
               addlistener(type_lines, ...
                       'Color','PostSet',...
                       @(~, evt) app.match_line_colors(evt, type_lines));
           end
           warning('on','MATLAB:table:RowsAddedExistingVars')
            
           % increase / decrease number of legend if needed (up to 3)
           if ~isempty(app.event_lg)
               nEntries = numel(app.event_lg.String);
               if nEntries > 6 && nEntries < 12
                   app.event_lg.NumColumns = 2;
               elseif nEntries >= 12
                   app.event_lg.NumColumns = 3;
               else
                   app.event_lg.NumColumns = 1;
               end
           end
       end

       function rmv_map_lines(app, events_type)
            % remove the lines linked with certain event_type.
            % if no event type is given, remove all that exist.
           
            % set default
            if ~exist("events_type","var")
                events_type = app.lines_types;
            end

            % error when given type that does not exist
            types_not_existing = ~ismember(events_type, app.lines_types);
            if any(types_not_existing)
                error("Event type {%s} does not exist in event lines",...
                    join(events_type(types_not_existing), ', '))
            end

            % delete graphics & table rows
            delete(app.event_lines{events_type,:})
            app.event_lines(events_type,:) = [];
       end
       
       function notify_line_pressed(app,evt)
           evt_data = general_evt_data();
           evt_data.data.evt = evt;
           notify(app, "line_pressed", evt_data)
       end
   end


   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%% Legend Callbacks %%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods
       % ItemHitFcn: event_ax_legend
        function event_ax_legend_on_click(app,evt)
            % Managed item click interactions -
            % change color on left click
           

            switch evt.SelectionType
                case 'normal'
                    % let user pick a new color for clicked element
                    current_color = evt.Peer.Color;
                    new_color = uisetcolor(current_color, evt.Peer.DisplayName + " color");
                    evt.Peer.Color = new_color; % note that due to listners, this should update matching lines

                    % for some reason, uisetcolor return focus to command
                    % window. return it to uifigure
                    focus(app.app_base.UIFigure)
            end


        end
   end


   methods (Hidden, Static)
       % Listener to Lines Color PostSet: global_map_lines, local_map_lines
       function match_line_colors(evt, lines2change)
           % convert colors of line2change to colors in the event
           [lines2change.Color] = deal(evt.AffectedObject.Color);
       end
   
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%% Utilites %%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   methods
       function reset_global_map_end(app)
           app.global_map_ax.InteractionOptions.RestoredXLimits = ...
               [0 app.model.longest_data];
           % zoom(app.global_map_ax, "out")
       end

   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%% Get Set Methods %%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods
       function res = get.local_exists(app)
           res = ~isempty(app.local_map_ax) && isvalid(app.local_map_ax);
       end

       function res = get.global_exists(app)
           res = ~isempty(app.global_map_ax) && isvalid(app.global_map_ax);
       end

       function res = get.all_lines(app)
           res = app.event_lines.Variables;
       end

       function res = get.lines_types(app)
           res = app.event_lines.Properties.RowNames;
           res = string(res(:)');
       end
   end
end