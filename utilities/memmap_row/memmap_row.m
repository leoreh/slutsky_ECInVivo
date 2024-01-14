classdef memmap_row < handle
    % wrapper to be able to split a memapfile to rows, and acess them
    properties
        mapped_file %(1,1) memmapfile
        row2indx    %(1,1) double
        data_name   %(1,:) string
        treat_as_double = true %(1,1) logical % no matter the underlining data type, output double
    end
    methods
        function obj = memmap_row(mapped_file, row2indx)
            % take a mapped file (using "memmapfile"), and what row to index to.
            % give way to index only this row using simple parantesis index.
            % 

            % simply validate that inputs are what they are expected to be, and
            % assign them into properties

            arguments
                mapped_file (1,1) memmapfile
                row2indx    (1,1) double
            end

            % validate that mapped_file has only 1 data in it
            data_names = fieldnames(mapped_file.Data);
            if numel(data_names) ~= 1
                error('Can seperate only files with 1 data name under them')
            end
            if numel(mapped_file.Data)~=1
                error('Can seperate only files with no repeates')
            end
            if mapped_file.Writable
                error('Writing to mapped file is not supported')
            end

            % save props
            obj.mapped_file = mapped_file;
            obj.row2indx    = row2indx;
            obj.data_name   = data_names{:};
        end
        
    end
    
    % object management methods
    methods(Hidden)
        function varargout = subsref(obj,s)
            switch s(1).type
                case '.'
                    [varargout{1:nargout}] = builtin('subsref',obj,s);
                case '()'
                    if length(s) == 1
                        % Implement obj(indices)
                        varargout{1} = obj.mapped_file.Data.(obj.data_name)(obj.row2indx,s.subs{:});
                        if iscellstr(s.subs) && s.subs == ":"
                            % due to the (obj.row2indx,s.subs{:}), output
                            % will be row even in the (:) case. force col
                            % vec in this case.
                            varargout{1} = varargout{1}';
                        end
                    else
                        % Use built-in for any other expression
                        [varargout{1:nargout}] = builtin('subsref',obj,s);
                    end
                case '{}'
                    [varargout{1:nargout}] = builtin('subsref',obj,s);
                otherwise
                    error('Not a valid indexing expression')
            end

            % if needed, convert numerics to double
            if obj.treat_as_double
                for iOutput = 1:numel(varargout)
                    if isnumeric(varargout{iOutput})
                        varargout{iOutput} = double(varargout{iOutput});
                    end
                end
            end
        end

        function ind = end(obj,k,n)
            % deal with passing "end" into the index
            sz = size(obj);
            if k < n
                ind = sz(k);
            else
                ind = prod(sz(k:end));
            end
        end

        function varargout = size(obj,varargin)
            % simply overloading size so it would treat the data instead of the obj

            % use builtin on data
            [varargout{1:nargout}] = builtin('size', obj.mapped_file.Data.(obj.data_name), varargin{:});

            % change 1st dim to 1
            if nargout <= 1 && isempty(varargin)
                % call type ans = size(obj)
                varargout{1}(1) = 1;
            elseif isempty(varargin)
                % call type [varargout] = size(obj)
                varargout{1} = 1;
            elseif nargout <= 1
                % call type ans = size(obj,varargin)
                varargout{1}([varargin{:}]==1) = 1;
            else
                % call type [varargout] = size(obj,varargin)
                [varargout{[varargin{:}]==1}] = deal(1);
            end
        end

        function big_dim = length(obj)
            if prod(size(obj))~=0 %#ok<PSIZE> % numel and size behave diffrently here
                big_dim = max(size(obj));
            else
                big_dim = 0;
            end
        end
    
    end
    
    % Built-in that need help dumping all of the memory data into them.
    % Always same pattern: run through inputs, collect memmap_row data,
    % place in built-in
    methods(Hidden)
        % basic math
        function res = plus(varargin)
            varargin = collect_all_rows(varargin);
            res = builtin("plus",varargin{:});
        end
        
        function res = minus(varargin)
            varargin = collect_all_rows(varargin);
            res = builtin("minus",varargin{:});
        end
        
        function res = uminus(varargin)
            varargin = collect_all_rows(varargin);
            res = builtin("uminus",varargin{:});
        end
        
        function res = mtimes(varargin)
            varargin = collect_all_rows(varargin);
            res = builtin("mtimes",varargin{:});
        end
        
        function res = times(varargin)
            varargin = collect_all_rows(varargin);
            res = builtin("times",varargin{:});
        end
        
        function res = mldivide(varargin)
            varargin = collect_all_rows(varargin);
            res = builtin("mldivide",varargin{:});
        end
        
        function res = mrdivide(varargin)
            varargin = collect_all_rows(varargin);
            res = builtin("mrdivide",varargin{:});
        end
        
        function res = ldivide(varargin)
            varargin = collect_all_rows(varargin);
            res = builtin("ldivide",varargin{:});
        end
        
        function res = rdivide(varargin)
            varargin = collect_all_rows(varargin);
            res = builtin("rdivide",varargin{:});
        end
        
        function res = power(varargin)
            varargin = collect_all_rows(varargin);
            res = builtin("power",varargin{:});
        end
        
        function res = sum(varargin)
            varargin = collect_all_rows(varargin);
            res = builtin("sum",varargin{:});
        end
   
        
        % derivative & integral
        function res = diff(varargin)
            varargin = collect_all_rows(varargin);
            res = builtin("diff",varargin{:});
        end
        
        function res = cumsum(varargin)
            varargin = collect_all_rows(varargin);
            res = builtin("cumsum",varargin{:});
        end

        
        % basic logic
        function res = eq(varargin)
            varargin = collect_all_rows(varargin);
            res = builtin("eq",varargin{:});
        end
        
        function res = ne(varargin)
            varargin = collect_all_rows(varargin);
            res = builtin("ne",varargin{:});
        end
        
        function res = lt(varargin)
            varargin = collect_all_rows(varargin);
            res = builtin("lt",varargin{:});
        end
        
        function res = gt(varargin)
            varargin = collect_all_rows(varargin);
            res = builtin("gt",varargin{:});
        end
        
        function res = le(varargin)
            varargin = collect_all_rows(varargin);
            res = builtin("le",varargin{:});
        end
        
        function res = ge(varargin)
            varargin = collect_all_rows(varargin);
            res = builtin("ge",varargin{:});
        end
        

        % sort-based function
        function varargout = min(varargin)
            varargin = collect_all_rows(varargin);
            [varargout{1:max(1,nargout)}] = builtin("min",varargin{:});
        end
        
        function varargout = max(varargin)
            varargin = collect_all_rows(varargin);
            [varargout{1:max(1,nargout)}] = builtin("max",varargin{:});
        end

        function varargout = sort(varargin)
            varargin = collect_all_rows(varargin);
            [varargout{1:max(1,nargout)}] = builtin("sort",varargin{:});
        end
        

        % shape manipulations
        function varargout = transpose(varargin)
            varargin = collect_all_rows(varargin);
            [varargout{1:max(1,nargout)}] = builtin("transpose",varargin{:});
        end
        
        function varargout = ctranspose(varargin)
            varargin = collect_all_rows(varargin);
            [varargout{1:max(1,nargout)}] = builtin("ctranspose",varargin{:});
        end
        

        % graphics
        function res = plot(varargin)
            varargin = collect_all_rows(varargin);
            
            if nargout == 0
                % do not define res if no outputs were asked
                builtin("plot",varargin{:});
            else
                res = builtin("plot",varargin{:});
            end
            
        end
        

        % convert to other numeric types
        function res = double(obj)
            res = builtin("double",obj.mapped_file.Data.(obj.data_name));
        end
    end
    
    % type tests answers for object
    methods (Hidden)
        function answer = isnumeric(~)
            % force everything to treat this object as numeric
            answer = true;
        end
    end


    
end