classdef IED_data_ax < data_ax_view

    properties
        thr_line cell_dict = cell_dict(string.empty, {}) % str > cell containing 2 lines (for high & low threshold)
        detection_marker cell_dict = cell_dict(string.empty, {}) % str > cell containing scatter

        detection_changed_listener event.listener
    end

    methods
        function set_view(app)
            % add listener to "detection_changed" event
            set_view@data_ax_view(app)

            % change marker when model changed it
            app.detection_changed_listener  = listener(app.model, "detection_changed", @(~,~) app.replot_marker);
        end

        function make_data_axes(app, data_ax_p, row_loc, col_loc, data_labels, disp_XAxis)
            % add threshold line to axes based on model.ied.signal2use
            % add marker of current detection point in all signals

            % create all axes as super
            make_data_axes@data_ax_view(app, data_ax_p, row_loc, col_loc, data_labels, disp_XAxis)

            % add markers, and lines as needed
            sig2use = app.model.ied.sig2use;
            for iAx = 1:numel(data_labels)
                iLabel = data_labels(iAx);
                if contains(iLabel, sig2use)
                    % make threshold lines for signals that threshold was
                    % analysed on. Note that
                    app.thr_line{iLabel} = [...
                        matlab.graphics.chart.primitive.Line(...
                        'Parent', app.data_axes{iLabel}, ...
                        'Color','r','LineWidth',1, 'LineStyle', '--', ...
                        'Tag',"Thr "+ iLabel)
                        
                        matlab.graphics.chart.primitive.Line(...
                        'Parent', app.data_axes{iLabel}, ...
                        'Color','r','LineWidth',1, 'LineStyle', '--', ...
                        'Tag',"-Thr "+ iLabel)...
                        ];
                end

                app.detection_marker{iLabel} = matlab.graphics.chart.primitive.Scatter(...
                    'Parent', app.data_axes{iLabel}, ...
                    'Marker', '*', 'CData',[1 0 0], 'SizeData', 50, ...
                    'Tag',iLabel);
                uistack(app.detection_marker{iLabel},"top")
            end
        end

        function replot_data(app)
            replot_data@data_ax_view(app)

            % modify data displayed in thr_lines
            thr_labels = string(app.thr_line.keys);
            thr_labels = thr_labels(:)'; % force row vec
            for iLabel = thr_labels
                % collect threshold values for this axe
                % assume label is the same for axe & data
                req_win = app.model.data2disp.req_time_window{iLabel};
                thr_val = app.model.collect_threshold(req_win);

                X_data = app.model.data2disp.actual_time_window{iLabel};
                X_data_full = linspace(X_data(1), X_data(2), numel(thr_val));
                [app.thr_line{iLabel}(1).XData, app.thr_line{iLabel}(1).YData] = ...
                    deal(X_data_full, thr_val);

                [app.thr_line{iLabel}(2).XData, app.thr_line{iLabel}(2).YData] = ...
                    deal(X_data_full, -thr_val);
            end
        end

        function replot_marker(app)
            % change marker location when model moved marker

            axes_labels = string(app.data_axes.keys);
            detection_time = app.model.collect_detection();
            for iLabel = app.model.data_labels
                % skip if label do not exist in axes
                if ~ismember(iLabel, axes_labels)
                    continue
                end
                
                % find the closest sample in this dataset
                fs = app.model.data_sources.fs(iLabel);

                % find the sample matching the detection
                matching_samp = round(detection_time.*fs);
                % make sure this samp is not out of data
                nData_points = length(app.model.data_sources.data{iLabel});
                if (matching_samp > nData_points) || (matching_samp < 0)
                    matching_samp = [];
                end
                XData = matching_samp./fs;
                YData = app.model.data_sources.data{iLabel}(matching_samp);

                % push to marker
                [app.detection_marker{iLabel}.XData, app.detection_marker{iLabel}.YData] = ...
                    deal(XData, YData);
            end
        end
    end

end