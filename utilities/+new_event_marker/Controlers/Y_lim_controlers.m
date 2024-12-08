classdef Y_lim_controlers < handle

    properties
        app_base
    end

    properties
        data_axes      cell_dict = cell_dict(string.empty, {}) % str > cell containing axe
        axes_editboxes cell_dict = cell_dict(string.empty, {}) % str > cell containing editbox

        zoom_together_tree  matlab.ui.container.CheckBoxTree
        axes_checkboxes cell_dict = cell_dict(string.empty, {}) % str > cell containing TreeNode
        check_all               matlab.ui.container.TreeNode
    end

    methods
        function app = Y_lim_controlers(app_base)
            % connect controler to model & base

            arguments
                app_base
            end

            % save handles to main objects
            app.app_base = app_base;
        end

        function set_controlers(app, options)
            arguments
                app

                options.zoom_together_tree
            end

            %%%%%%% Connect controlers with defaults
            
            % convert options to cell_dict
            components_inputs = cell_dict( string(fieldnames(options)), struct2cell(options));
            
            % create components_info_tab
            components_names = "zoom_together_tree";
            def_components = app.def_control("zoom_together_tree");
            all_callbacks = {@(~,~) app.Ylink_axes};
            callbacks_names = "CheckedNodesChangedFcn";
            components_info_tab = table(def_components, all_callbacks, callbacks_names, ...
                'RowNames',components_names);

            % push to app
            add_create_components(components_inputs, app, components_info_tab);

            % if default controler was used, insert everything under
            % "check_all" checkbox
            if isfield(options, "zoom_together_tree") && iscell(options.zoom_together_tree)
                app.check_all = app.zoom_together_tree.Children(1);
            end
            

        end
    end

    methods
        function add_ax(app, axes_cell_dict, editboxes_in, checkbox_flag)
            % add axes to Y controler, connect them to editboxes & checkbox
            % axes_cell_dict - cell_dict, cell_dict, with keys as ax_labels.
            %                  values are axes connected to those labels.
            % editboxes_in   - cell_dict, with keys as ax_labels. Note that
            %                  while editboxes_in keys must exist in
            %                  axes_cell_dict, the other way is not required.
            %                  values either editbox handle or cell, see add_create_components.
            % checkbox_flag  - logical vector, if to create checkbox matching axe or not
            
            %%%%% validate inputs
            axes_labels = string(axes_cell_dict.keys);
            axes_labels = axes_labels(:);

            % make sure all labels are in editboxes_in
            if any( ~ismember(string(editboxes_in.keys), axes_labels))
                error('fieldnames in editboxes_in must be axes labels (axes_cell_dict keys)')
            end

            % make sure there is a flag for each entry
            if numel(axes_labels) ~= numel(checkbox_flag)
                error('checkbox_flag must have 1 flag for each axe in axes_cell_dict')
            end
            
            %%%%% create & connect componants
            % enter axes to axes dictionary
            app.data_axes(axes_labels) = axes_cell_dict(axes_labels);

            % add checkboxes
            app.add_checkboxes(axes_labels(checkbox_flag))

            % add editboxes
            app.add_editboxes(editboxes_in, axes_cell_dict)
        end

        function add_checkboxes(app, ax_labels)
            % data_labels - string vector, which labels 2 add

            nLabels = numel(ax_labels);
            for iLabel = 1:nLabels
                now_data = ax_labels(iLabel);
                
                % find where to app - under check all, or in tree
                if ~isempty(app.check_all) && isvalid(app.check_all)
                    node_parent = app.check_all;
                else
                    node_parent = app.zoom_together_tree;
                end
                
                % add treenode
                app.axes_checkboxes{now_data} = ...
                    matlab.ui.container.TreeNode('Parent', node_parent, ...
                    'Text', now_data);

                % remove checkbox if axe is deleted
                addlistener(app.data_axes{now_data}, "ObjectBeingDestroyed",...
                    @(~,~) delete(app.axes_checkboxes{now_data}));
            end

            % expend tree
            expand(app.zoom_together_tree)
        end

        function add_editboxes(app, editboxes_in, axes_cell_dict)
            % editboxes_in - struct, with fieldnames as ax_labels.
            %                values either editbox handle or cell, see add_create_components
            
            ax_labels = string(keys(editboxes_in));
            ax_labels = ax_labels(:)';

            componants_info = table();
            warning('off','MATLAB:table:RowsAddedExistingVars')
            for iAx = ax_labels
                % create default componant
                new_ctrl = app.def_control("axes_editboxes");
                box_label = sprintf('%s ylims',iAx);
                new_ctrl.Tooltip = {box_label};
                new_ctrl.Placeholder = box_label;
                
                % place in table
                componants_info.def_components(iAx) = new_ctrl;
                componants_info.all_callbacks(iAx) = {@(~,evt) app.lim_box_change(evt, iAx)};
                componants_info.callbacks_names(iAx) = "ValueChangedFcn";
            end
            warning('on','MATLAB:table:RowsAddedExistingVars')

            % place in a struct
            comp_dict = add_create_components(editboxes_in, [], componants_info);

            % transform to cell_dict
            app.axes_editboxes(ax_labels) = comp_dict(ax_labels);

            % reverse databind axes & editboxes
            for iAx = ax_labels
                axes_cell_dict{iAx}.YAxis.LimitsChangedFcn = @(~,~) app.update_lim_editbox(iAx);
                app.update_lim_editbox(iAx);

                % remove checkbox if connected axe is deleted
                addlistener(axes_cell_dict{iAx}, "ObjectBeingDestroyed",...
                    @(~,~) delete(app.axes_editboxes{iAx}));
            end
        end
    end
    
    methods
        function Ylink_axes(app)
            % get from checkboxes what data should be y linked, and link
            % them

            % unlink axes
            linkaxes([app.data_axes{string(app.axes_editboxes.keys)}],'off');

            if isempty(app.zoom_together_tree.CheckedNodes)
                return
            else
                % collect from checkbox what 2 'y' sync
                data_with_checkboxes = string(app.axes_checkboxes.keys);
                checked_nodes = ismember(...
                    data_with_checkboxes,...
                    {app.zoom_together_tree.CheckedNodes.Text});

                % sync if needed
                axes2link = data_with_checkboxes(checked_nodes);
                linkaxes([app.data_axes{axes2link}],'y')
            end
            

        end

        function lim_box_change(app, evt, ax_label)

            % parse text 2 numbers
            in_parsed = parse_2num_editbox(evt.Value);

            % validate input
            if numel(in_parsed) ~= 2 || any(isnan(in_parsed)) || any(isinf(in_parsed)) || ...
                    ~issorted(in_parsed,"strictmonotonic")
                % wrong num of elements or bad transfer or bad range- flash red, restore
                app.axes_editboxes{ax_label}.BackgroundColor = 'r';
                pause(0.1)
                app.axes_editboxes{ax_label}.BackgroundColor = 'w';
                app.axes_editboxes{ax_label}.Value = evt.PreviousValue;
            else
                % flash green, sort values, change YLim, make text look nice
                app.axes_editboxes{ax_label}.BackgroundColor = 'g';
                app.data_axes{ax_label}.YLim = in_parsed;
                pause(0.1)
                app.axes_editboxes{ax_label}.BackgroundColor = 'w';
                app.axes_editboxes{ax_label}.Value = num2str(in_parsed,'[%.3f %.3f]');

                % move focus to figure
                focus(app.app_base.UIFigure)
            end

        end
    
        % Ylim listner: data_axes{iData}
        function update_lim_editbox(app,iAx)
            % simply fix the editbox to match the new ylim value
            app.axes_editboxes{iAx}.Value = num2str(app.data_axes{iAx}.YLim,'[%.3f %.3f]');
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Keypress Responce Callbacks %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function [help_text] = key_responce(app, key_pressed)
            % create keyboard shortcuts

            % create help text
            if nargout > 0
                help_text = ...
                    "control+d -> zoom out by 10% on all axes" + newline + ...
                    "control+i -> zoom in by 10% on all axes";
                return
            end

            % create key responce
            all_keys = string(app.data_axes.keys());
            all_keys = all_keys(:)';
            switch key_pressed
                case "control+d"
                    % zoom out by 10% on all axes
                    for iAx = all_keys
                        ylims = app.data_axes{iAx}.YLim;
                        y_win = diff(ylims);
                        y_cent = mean(ylims);
                        new_ylims = y_cent + ([-y_win y_win]/2)*1.1;
                        app.data_axes{iAx}.YLim = new_ylims;
                    end
                case "control+i"
                    % zoom in by 10% on all axes
                    for iAx = all_keys
                        ylims = app.data_axes{iAx}.YLim;
                        y_win = diff(ylims);
                        y_cent = mean(ylims);
                        new_ylims = y_cent + ([-y_win y_win]/2)*0.9;
                        app.data_axes{iAx}.YLim = new_ylims;
                    end
            end
        end
    
    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% Helper Functions %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Hidden)
        function control_obj = def_control(~, obj_name)
            switch obj_name
                case "zoom_together_tree"
                    control_obj = matlab.ui.container.CheckBoxTree();
                    control_obj.Tooltip = 'Choose which dataset have linked Y axis';
                    % app.zoom_together_checkbox.CheckedNodesChangedFcn = @(~,~) zoom_together_checkbox_checked_nodes_changed(app);

                    check_all = matlab.ui.container.TreeNode('Parent',control_obj);
                    check_all.Text = 'All';
                case "axes_editboxes"
                    control_obj = matlab.ui.control.EditField();
                    control_obj.HorizontalAlignment = 'center';
            end
        end
    
    end
end