function [obj] = add_create_components(components_inputs, obj, components_info)
    % ADD_CREATE_COMPONENTS Add and create GUI components to object
    %
    % add_create_components(components_inputs, obj, default_components, all_callbacks)
    %
    % This function adds graphical user interface (GUI) components specified in
    % the components_inputs struct to the obj handle. It creates default component
    % objects if needed or uses user-provided handles. The components are added
    % to the layout locations specified and callbacks are connected.
    %
    % INPUTS:
    %   components_inputs - 
    %           cell_dict (string => cell) containing keys that are obj properties,
    %           and values that are either:
    %               1) Cell array {Parent, Row, Column} to place default component
    %                  Parent must be a handle for a grid layout container (uigridlayout)
    %               2) Handle to user created object
    %
    %   obj -   Handle to object to add components too, or empty
    %
    %   components_info_tab - 
    %           Table, containing components info. Row names are obj
    %           properties that the component should be saved as.
    %           Table have the variables:
    %
    %       def_components - 
    %               handles to default components
    % 
    %       all_callbacks - 
    %               Cell array of function handles, callbacks to
    %               connect to components
    % 
    %       callbacks_names - 
    %               Text vector, which callback to add the callback to.
    %
    %  OUTPUTS:
    %   obj - handle to input object , or cell_dict with same keys as
    %         components_inputs 
    %
    % Example:
    %
    % % create object to insert componants to
    % f = uifigure();
    % addprop(f,"editbox");
    % addprop(f,"uidatepicker");
    % base_grid = uigridlayout(f,[1,2]);
    %
    % % create struct holding componants data
    % components = cell_dict(string.empty,{});
    % components{"editbox"} = uieditfield(f,"Position",[0 0 100 400],"Value","5");
    % components{"uidatepicker"} = {base_grid, 1, 2};
    %
    % % create a table holding components general info
    % components_names = ["editbox" ; "uidatepicker"];
    % def_components = [matlab.ui.control.EditField ; matlab.ui.control.DatePicker];
    % all_callbacks = {@(~,~)disp("Editbox Changed") ; @(~,~)disp("Datepicker Changed") };
    % callbacks_names = ["ValueChangedFcn" ; "ValueChangedFcn"];
    % components_info = table(def_components, all_callbacks, callbacks_names, 'RowNames',components_names);
    %
    % add_create_components(components, f, components_info);
    %
    % % note that you can remove any of the fields, and the function will
    % % still run OK.
    
    arguments
        components_inputs cell_dict
        obj
        components_info table
        % default_components matlab.graphics.Graphics
        % all_callbacks cell
        % callbacks_names string
    end

    % validate table
    var_names = components_info.Properties.VariableNames;
    try
        good_names = ["def_components", "all_callbacks", "callbacks_names"];
        for iVar = 3:-1:1
            new_vars(iVar) = validatestring(var_names{iVar}, good_names);
        end
        components_info.Properties.VariableNames = new_vars;
    catch ME
        if ME.identifier == "MATLAB:unrecognizedStringChoice"
            error(replace(ME.message,"input","variable name"))
        else
            rethrow(ME)
        end
    end
    
    % if no object was given, return cell_dict
    if isempty(obj)
        obj = cell_dict(string.empty,{});
    end
    
    objcts2add = keys(components_inputs);
    for iComponent = string(objcts2add)'
        %%% insert component to object
        if iscell(components_inputs{iComponent})
            % take the created default object, place it in the location specified
            % in the cell array
            current_componant = components_info.def_components(iComponent);
            current_componant.Parent = components_inputs{iComponent}{1};
            current_componant.Layout.Row = components_inputs{iComponent}{2};
            current_componant.Layout.Column  = components_inputs{iComponent}{3};

        elseif isa(components_inputs{iComponent}, class(components_info.def_components(iComponent)))
            % use the user given handle
            current_componant = components_inputs{iComponent};

        else
            error("%s is expected to be of class %s, or cell array {GridLayoutParent, Row, Col}")
        end

        %%% add callback
        callback_name = components_info.callbacks_names(iComponent);
        current_componant.(callback_name) = components_info.all_callbacks{iComponent};

        %%% push to object
        if ~isa(obj,'handle')
            % cell_dict case
            obj{iComponent} = current_componant;
        else
            % object handle case
            obj.(iComponent) = current_componant;
        end
    end
end