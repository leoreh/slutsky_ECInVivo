classdef event_marker_4_reduct_displayer < event_marker

    properties
        lines_context_menu              matlab.ui.container.ContextMenu
        rmv_event_menu                  matlab.ui.container.Menu
        show_on_reduce_axe_menu         matlab.ui.container.Menu

        reduce_displayer                reduct_displayer
        event_center_sample (1,1)       struct
    end
     
    % constructor
    methods
        function app = event_marker_4_reduct_displayer(reduce_displayer_obj, event_center_sample, varargin)
            
            % build the app, add properties
            app@event_marker(varargin{:})
            app.reduce_displayer = reduce_displayer_obj;
            app.event_center_sample = event_center_sample;
            
            % change to unsaveable
            app.saveable = false;
            
            % add context menu to event lines
            app.lines_context_menu = uicontextmenu(app.UIFigure);
            % app.rmv_event_menu = uimenu(app.lines_context_menu,"Text","Remove event");
            app.show_on_reduce_axe_menu = uimenu(app.lines_context_menu,"Text","Show event on reduction mqp");
            [app.event_line.ContextMenu] = deal(app.lines_context_menu);
        end
    end
    
    % superclass override
    methods
        function event_clicked(app,evt)
            % overwriting event marker for behaviors that match
            % reduct_displayer - open context menu and pass the Hit eveent 

            switch app.UIFigure.SelectionType
                case 'alt'
                    % collect event type
                    event_type = evt.Source.DisplayName;

                    % collect time clicked and find which event it is refering
                    t_clicked = evt.IntersectionPoint(1);
                    events_starts = app.events_time.(event_type)(:,1);
                    events_ends   = app.events_time.(event_type)(:,2);
                    relv_event = t_clicked >= events_starts & t_clicked <= events_ends;

                    % pass event type & which event to context menus
                    % app.rmv_event_menu.MenuSelectedFcn = @(~,~) app.rmv_event(event_type,relv_event);
                    app.show_on_reduce_axe_menu.MenuSelectedFcn = @(~,~) app.show_on_reduce_axe(event_type,relv_event);
                case 'normal' % left click
                    datatip(evt.Source,evt.IntersectionPoint(1),evt.IntersectionPoint(2));
            end

        end
    end
    
    % menu functions
    methods
        function rmv_event(app,event_type,relv_event)
            % remove cliked on event
            
            % insert events without deleted
            events_starts(relv_event) = [];
            events_ends(relv_event) = [];
            app.events_time.(event_type) = sortrows([events_starts, events_ends],[1 2],"ascend");
            
            % remove event sample
            app.event_center_sample.(event_type)(relv_event) = [];
            
            % refresh plot
            app.plot_time(app.event_ax.XLim)
            app.update_event_map()
        end

        function show_on_reduce_axe(app,event_type,relv_event)
            % show clicked event on reduct_displayer reduction_axe

            % collect event sample
            relevant_samp = app.event_center_sample.(event_type)(relv_event);

            % if more than 1 event was selected, take only the first one
            relevant_samp = relevant_samp(1);

            % pass to relevant function in reduct_displayer
            app.reduce_displayer.wv_on_other_disp_MenuSelected(relevant_samp,matlab.graphics.chart.primitive.Line)

        end
    end
    
end