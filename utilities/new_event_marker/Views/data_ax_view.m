classdef data_ax_view < handle
    % create axes & lines for data view. This includes:
    % 1) data_axes, axes with 1 line. Used to display the actual signal.
    % 2) event_ax, axes with 1 line per event type.

    properties
        app_base
        model data2disp_model = data2disp_general_model()
    end
    
    properties
        data_axes   cell_dict = cell_dict(string.empty, {}) % str > cell containing axe
        data_lines  cell_dict = cell_dict(string.empty, {}) % str > cell containing line

        in_win_change_listener event.listener
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% Setup, Adding & Removing Componants %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function app = data_ax_view(app_base, model)
            % connect view to model & base
            
            arguments
                app_base
                model data2disp_model
            end

            % save handles to main objects
            app.app_base = app_base;
            app.model = model;
        end

        function set_view(app)
            % prepare view to work -
            % add event_ax, listen to "in_win_change" in model.

            % change view when model changed
            app.in_win_change_listener = listener(app.model, "in_win_change", @(~,~) app.replot_data);
        end

        function make_data_axes(app, data_ax_p, row_loc, col_loc, data_labels, disp_XAxis)
            % make new axes to match data. Replace old axes. used by
            % app_base to dynamicly add data upon "data_loaded" event.
            %
            % INPUTS:
            %   data_ax_p - Array of handles to grid layout containers (uigridlayout).
            %               Parent for each new axe.
            %   row_loc, col_loc -
            %               Cell arrays, for each axe its position in its
            %               parent layout.
            %   data_labels -
            %               String row vec, for each axe how to title it.
            %               Note that currently "replot_data" assume same
            %               label for axes & model, if you want diffrent
            %               names you need to add a converter.
            %   disp_XAxis - 
            %               Logical flags, for each data_label if to present the XAxis or not.
            %               Default: true for all (show all)
            % Note that all inputs must have the same number of elements.

            % make sure all inputs have the same number of elements
            nReq = numel(data_ax_p);
            if any([numel(row_loc), numel(col_loc), numel(data_labels)] ~= nReq)
                error('all inputs must have the same number of elements')
            end
            
            if any(~ismember(data_labels, app.model.data_labels))
                error('Using axes labels that do not exist in model is currently unsupported. To support that, create a converter.')
            end

            % default disp_XAxis to true for all
            if ~exist("disp_XAxis","var") || isempty(disp_XAxis)
                disp_XAxis = true(nReq,1);
            end

            % create objects
            for iAx = numel(data_labels):-1:1
                current_label = data_labels(iAx);
                % create & place axe
                 ax = uiaxes(...
                    'Parent', data_ax_p(iAx),...
                    'YLimMode','manual', 'Clipping','off',...
                    "PositionConstraint","innerposition",...
                    'Tag',data_labels{iAx});
                ax.Layout.Row = row_loc{iAx};
                ax.Layout.Column = col_loc{iAx};
                
                % manage how axe look & interact
                label_obj = xlabel(ax, current_label, ...
                    "FontWeight","bold","FontSize",14);
                label_obj.Visible = 'on';
                ax.XLimitMethod = 'tight';
                ax.XAxis.Exponent = 0;
                % ax.XAxis.Visible = 'off';
                ax.Interactions = dataTipInteraction;
                ax.YTick = 0; % This is important to solve a resizing bug - all axes in the app need to have the same Y ticks

                % set axis ylims to match +-2 std in a 30 min / 10% period in mid recording
                rec_len = app.model.longest_data;
                mid_rec = rec_len./2;
                sample_size = min(60*30, rec_len.*0.1); % do not sample more than 30 minutes
                sample_size = max(sample_size,5); % at least sample 10 secound
                sample_period = mid_rec + [-sample_size sample_size]/2;
                mid_sig = app.model.collect_data(current_label, sample_period);
                avg_sig = mean(mid_sig);
                std_sig = std(mid_sig);
                if std_sig == 0
                    % deal with constant signal
                    std_sig = 0.1;
                end
                max_range = [avg_sig-10*std_sig, avg_sig+10*std_sig];
                ax.YLim = max_range;

                % create line in axe
                app.data_lines{current_label} = matlab.graphics.chart.primitive.Line(...
                    'Parent', ax, ...
                    'Color','b','LineWidth',1, 'Tag',current_label);
                
                % hide XAxis if requested
                if ~disp_XAxis(iAx)
                    ax.XAxis.Visible = "off";
                end

                % save ax in app
                app.data_axes{current_label} = ax;
            end
        end

        function rmv_data_axes(app, data_labels)
            % remove the axes linked with certain dataset.
            % if no data label is given, remove all that exist in dataset
            if ~exist("data_labels","var")
                data_labels = string(keys(app.data_axes));
            end
            
            % delete existing axes & clear cell_dicts
            delete([app.data_axes{data_labels}])
            app.data_axes(data_labels) = [];
            app.data_lines(data_labels) = [];
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% Model Window Change Response %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function replot_data(app)
            % push data in window to displays.
            % Assume that labels are the same between model and axes.
            
            axes_labels = string(app.data_axes.keys);

            % modify displayed data in data_axes
            for iLabel = app.model.data_labels
                % skip if label do not exist in axes
                if ~ismember(iLabel, axes_labels)
                    continue
                end
                
                % collect data, find matching time points
                Y_data = app.model.data2disp.data{iLabel};
                X_data = app.model.data2disp.actual_time_window{iLabel};
                X_data_full = linspace(X_data(1), X_data(2), numel(Y_data));

                % push to plot
                [app.data_lines{iLabel}.XData, ...
                    app.data_lines{iLabel}.YData] = ...
                    deal(X_data_full,Y_data);

                % update X limits to match
                app.data_axes{iLabel}.XLim = app.model.data2disp.req_time_window{iLabel};
            end

        end
    
    end
end