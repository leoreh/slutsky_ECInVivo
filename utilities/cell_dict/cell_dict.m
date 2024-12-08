classdef cell_dict < matlab.mixin.indexing.RedefinesBrace & ...
        matlab.mixin.indexing.RedefinesParen
    
    properties
        data
    end

    methods
        function obj = cell_dict(varargin)
            %CELL_DICT Wrap container.Map/dictionary to provide consistent
            %indexing when value type is any/cell.
            % NOTES:
            %   * if source data is container.map, this is basicly a handle
            %     class. If it is dictionary, it is a value class. Use with
            %     care.
            %   * If no existing object is given, will prefere using
            %     dictionary over container.map, based on MATLAB version.
            %   * Generally, if you don't care about older MATLAB versions,
            %     use dictionary directly.
            %   
            % SYNTAX:
            %   obj = cell_dict()                   % no entries
            %   obj = cell_dict(data)               % convert existing object
            %   obj = cell_dict(keys, values)       % create from keys & values
            %   obj = cell_dict(k1,v1,...,kN,vN)    % create from keys & values
            %
            % DESCRIPTION:
            %   Constructor to initialize the cell_dict with data. Supports
            %   creating a cell_dict with no entries, initializing from an existing
            %   containers.Map/dictionary, or creating from key/value pairs.
            %
            % Examples:
            % Create an empty cell_dict 
            % d = cell_dict(); 
            % 
            % % Create from key/value pairs
            % d = cell_dict("key1", {[1, 2]}, "key2", {{'a','b'}});
            % 
            % % Access values
            % val = d{"key1"}; 
            % 
            % % Access elements of value cell arrays
            % val1 = d{"key1"}(1);
            % 
            % % Assign values
            % d{"key3"} = {3, 4};
            % 
            % 
            % % Display contents
            % disp(d);
            % 
            % % Initialize from existing containers.Map/dictionary
            % data = containers.Map({'k1','k2'}, {[1 2], 'val'});
            % d = cell_dict(data);
            %
            % See Also dictionary, containers.Map 
            %
            % Author: Lior de Marcas (LdM: https://github.com/Liordemarcas)
            % Date: 2024-02-15

            switch nargin()
                case 0
                    % create an empty container
                    if isMATLABReleaseOlderThan("R2022b")
                        obj.data = containers.Map;
                    else
                        obj.data = dictionary;
                    end
                case 1
                    % input is the already created object
                    obj.data = varargin{1};
                otherwise
                    % input is (keys, values) or (key,  value, key, value ...)
                    % responce is simillar in both cases.
                    if isMATLABReleaseOlderThan("R2022b")
                        keys = [varargin{1:2:end}];
                        if ischar(keys) || iscellstr(keys) %#ok<ISCLSTR> 
                            keys = string(keys);
                        end
                        values = [varargin{2:2:end}];
                        if any(size(keys) ~= size(values))
                            % use dictionary built in error for this case
                            error('MATLAB:dictionary:KeyValueDimsMustMatch',...
                                'Dimensions of the key and value must be compatible, or the value must be scalar.')
                        end
                        if isscalar(values)
                            values = [values{:}];
                        end
                        if isempty(keys) && isempty(values)
                            obj.data = containers.Map;
                        else
                            obj.data = containers.Map(keys,values,'UniformValues',false);
                        end
                    else
                        obj.data = dictionary(varargin{:});
                    end
            end

            % make sure value type is cell
            [~,valueType] = types(obj);
            if valueType ~= "cell" && (numEntries(obj) ~= 0)
                error('Value type must be cell!')
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Implement  parentheses methods %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function out = cat(varargin) %#ok<STOUT> % should just error
            error('MATLAB:class:concatenationScalar', ...
                "Concatenating 'cell_dict2' objects is not supported because 'cell_dict' objects can only be scalar.")
        end
    
        function varargout = size(obj,varargin)
             [varargout{1:nargout}] = builtin('size',obj,varargin{:});
        end
    end

    methods (Static)
        function obj = empty()
            % no real empty - just init an object with an empty dict
            obj = cell_dict();
        end
    end

    methods (Access = protected)
        function varargout = parenReference(obj, indexOp)
            if isa(obj.data, 'containers.Map')
                % use brace indexing
                keys = [indexOp(1).Indices{:}];
                [varargout{1}{1:numel(keys)}] = deal(obj{keys});
            else
                % faster dictionary responce, let multiple outputs
                [varargout{1:nargout}] = obj.data.(indexOp);
            end
        end
        
        function obj = parenAssign(obj, indexOp, in_cell)
            if ~iscell(in_cell)
                error('MATLAB:dictionary:CannotSliceValues', ...
                    "Unable to use '%s' as value for dictionary with 'cell' value type.", ...
                    class(in_cell))
            end
            
            if isa(obj.data, 'containers.Map')
                % use brace indexing
                keys = [indexOp(1).Indices{:}];
                [obj{keys}] = deal(in_cell{:});
            else
                % faster dictionary responce
                obj.data.(indexOp) = in_cell;
            end
        end

        function obj = parenDelete(obj, indexOp)
            if isa(obj.data, 'containers.Map')
                keys = [indexOp(1).Indices{:}];
                if isstring(keys)
                    keys = cellstr(keys);
                end
                warning("off",'MATLAB:Containers:Map:NoKeyToRemove')
                remove(obj.data, keys);
                warning("on",'MATLAB:Containers:Map:NoKeyToRemove')
            else
                % faster dictionary responce
                obj.data.(indexOp) = [];
            end
        end

        function n = parenListLength(obj, indexOp, indexContext)
            if isa(obj.data, 'containers.Map')
                % containers.Map will always return 1 output
                n = 1;
            else
                % see how many dictionary will return
                n = listLength(obj.data,indexOp,indexContext);
            end
            % n = 1;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%% Implement Brace methods %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = protected)
        function varargout = braceReference(obj, indexOp)
            % Custom implementation of brace indexing
            key = [indexOp(1).Indices{:}];
            if ~isscalar(key) && ~isscalar(indexOp)
                % use matlab error for comma seperated
                error("MATLAB:index:expected_one_output_from_intermediate_indexing", ...
                    "Intermediate brace '{}' indexing produced a comma-separated list with 2 values, but it must produce a single value when followed by subsequent indexing operations.")
            end
            if isa(obj.data, 'containers.Map')
                
                % containers.Map don't play well with strings, transfer to char
                if isstring(key)
                    key = cellstr(key);
                end

                % If the data is a containers.Map, access it directly
                for iKey = nargout:-1:1
                    varargout{iKey} = obj.data(key{iKey});
                end
            elseif isMATLABReleaseOlderThan("R2023b")
                % Older dictionary, use 2 step referencing to avoid error
                value = obj.data(key);
                [varargout{1:nargout}] = value{:};
            else
                % New dict type
                [varargout{1:nargout}] = obj.data{key};
            end

            if ~isscalar(indexOp)
                % because of error before - only 1 output
                varargout{1} = varargout{1}.(indexOp(2:end));
            end
        end
        
        function obj = braceAssign(obj, indexOp, varargin)
            % Custom implementation of brace assignment
            key = [indexOp(1).Indices{:}];
            if ~isscalar(key) && ~isscalar(indexOp)
                % use matlab error for comma seperated
                error("MATLAB:index:expected_one_output_from_intermediate_indexing", ...
                    "Intermediate brace '{}' indexing produced a comma-separated list with 2 values, but it must produce a single value when followed by subsequent indexing operations.")
            end
            
            if isa(obj.data, 'containers.Map')
                if ~isscalar(indexOp)
                    % multiple index operation - error block cases with
                    % more than 1 input
                    val = obj.data(key);
                    val.(indexOp(2:end)) = varargin{1};
                    obj.data(key) = val;
                else
                    % loop through
                    for iKey = numel(key):-1:1
                        obj.data(key(iKey)) = varargin{iKey};
                    end
                end
                
               
            elseif isMATLABReleaseOlderThan("R2023b")
                if ~isscalar(indexOp)
                    % multiple index operation - error block cases with
                    % more than 1 input.
                    % load subvar, convert it in place, than push back.
                    value = obj{key};
                    value.(indexOp(2:end)) = varargin{1};
                    obj.data(key) = {value};
                else
                    obj.data(key) = varargin;
                end
                
            else
                if ~isscalar(indexOp)
                    % multiple index operation - error block cases with
                    % more than 1 input
                    obj.data{key}.(indexOp(2:end)) = varargin{1};
                else
                    [obj.data{key}] = deal(varargin{:});
                end
            end

            
        end
        
        function len = braceListLength(~,indexOp,~)
            % Assume that wanted ouput is simply number of requested indecies
            len = numel(indexOp(1).Indices{:});
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% Implement dictionary methods %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods
        function obj = insert(obj, key, value, options)
            % See https://www.mathworks.com/help/matlab/ref/dictionary.insert.html
            
            arguments
                obj
                key
                value
                options.Overwrite (1,1) logical {mustBeNumericOrLogical} = true
            end
            
            % if user want to prevent overwrite, remove existing keys
            if ~options.Overwrite
                keys2rmv = isKey(obj, key);
                key(keys2rmv) = [];
                value(keys2rmv) = [];
            end

            % write new values
            obj(key) = value;
        end
    
        function obj = remove(obj, key)
            % See https://www.mathworks.com/help/matlab/ref/dictionary.remove.html
            obj(key) = [];
        end
        
        function k = keys(obj, out_type)
            % See https://www.mathworks.com/help/matlab/ref/dictionary.keys.html
            k = keys(obj.data);
            if exist("out_type","var") && strcmpi(out_type, "cell") && ~iscell(k)
                k = num2cell(k);
            end
        end
    
        function E = entries(obj, out_format)
            % See https://www.mathworks.com/help/matlab/ref/dictionary.keys.html
            
            % default to table format (how we collect)
            if nargin == 1
                out_format = "table";
            end

            % collect entries as table
            if isa(obj.data,"containers.Map")
                all_keys = keys(obj.data)';
                all_vals = values(obj.data)';
                E = table(all_keys,all_vals,'VariableNames',["Key","Value"]);
            else
                E = entries(obj.data);
            end
            
            % convert as needed
            if strcmpi(out_format,"struct")
                E = table2struct(E);
            elseif strcmpi(out_format,"cell")
                E = table2cell(E);
            end
        end
    
        function v = values(obj, out_type)
            % See https://www.mathworks.com/help/matlab/ref/dictionary.values.html

            v = values(obj.data);
            if exist("out_type","var") && strcmpi(out_type, "cell") && ~iscell(v)
                % should actully never be triggered, as valueType should always be
                % cell in this object. However keeping for preserving logic.
                v = num2cell(v);
            end
        end
    
        function tf = isKey(obj, key)
            % See https://www.mathworks.com/help/matlab/ref/dictionary.isKey.html
            if isstring(key) && isa(obj.data, 'containers.Map')
                key = cellstr(key);
            end
            tf = isKey(obj.data, key);
        end
    
        function n = numEntries(obj)
            % See https://www.mathworks.com/help/matlab/ref/dictionary.numEntries.html
            if isa(obj.data,"containers.Map")
                n = length(obj.data);
            else
                n = numEntries(obj.data);
            end
        end
    
        function [keyType,valueType] = types(obj)
            % See https://www.mathworks.com/help/matlab/ref/dictionary.types.html
           
            if isa(obj.data,"containers.Map")
                keyType = obj.data.KeyType;
                if keyType == "char"
                    % dictionary keytype is always converted to "string",
                    % containers.Map is the other way around.
                    keyType = "string";
                end

                valueType = obj.data.ValueType;
                if valueType == "any"
                    % "any" does not exist in dictionary keytype, but cell
                    % is the same
                    valueType = "cell";
                end

            else
                [keyType,valueType] = types(obj.data);
            end
        end
    
        function value = lookup(obj,key,options)
            % See https://www.mathworks.com/help/matlab/ref/dictionary.types.html
            
            arguments
                obj
                key
                options.FallbackValue
            end
            
            % since lookup does not exist for dictionary pre-2023b, we will
            % teat old dictionary the same as containers.Map
            if isMATLABReleaseOlderThan("R2023b")
                % loop & collect each value
                for iVal = numel(key):-1:1
                    try
                        value(iVal) = obj(key(iVal));
                    catch ME
                        if isfield(options, "FallbackValue") && ...
                                ismember(ME.identifier, ["MATLAB:Containers:Map:NoKey" "MATLAB:dictionary:ScalarKeyNotFound"])
                            % failed due to key not found, use fallback
                            value(iVal) = options.FallbackValue;
                        else
                            % no fallback OR other failer cause
                            rethrow(ME)
                        end
                    end
                end

            else
                if isfield(options, "FallbackValue")
                    value = lookup(obj.data, key, "FallbackValue", options.FallbackValue);
                else
                    value = lookup(obj.data, key);
                end
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%% Generate nice display %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function disp(obj)
            %DISP Custom command window display for cell_dict object
            
            fprintf('cell_dict that contains:\n')

            if isa(obj.data, 'containers.Map')
                % Collect keys and values
                all_entries = entries(obj,"table");
                [keyType,valueType] = types(obj);

                % Display in dictionary format
                fprintf('  %s (%s ⟼ %s) with %d entries:\n', ...
                    'containers.Map',keyType, valueType, numEntries(obj));

                % display each entry
                for iEntry = 1:numEntries(obj)
                        % collect values
                        key = all_entries.Key{iEntry};
                        value = all_entries.Value(iEntry);

                        % format values
                        key_str = strtrim(formattedDisplayText(key));
                        if keyType == "string"
                            % fix missing "" in char/string key type
                            key_str = '"' + key_str + '"';
                        end
                        value_str = strtrim(formattedDisplayText(value));

                        % print line
                        fprintf('    %s ⟼ %s\n', key_str, value_str);
                end

            else
                % just show dictionary format
                disp(obj.data)
            end
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% Register Obj as Dictionary %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % methods
    %     function res = isa(~, data_type)
    %         % return true when asked if this is a dictionary
    %         if any(ismember(data_type, ["dictionary", "cell_dict"]))
    %             res = true;
    %         else
    %             res = false;
    %         end
    %     end
    % 
    %     function mustBeA(obj, data_type)
    %         % prevent an error when asked if this is dictionary
    %         if any(ismember(data_type, "dictionary"))
    %             % do not raise an error
    %         else
    %             try
    %                 builtin("mustBeA",obj,data_type)
    %             catch ME
    %                 throwAsCaller(ME)
    %             end
    %         end
    %     end
    % end
end