classdef reduct_displayer_sig_display < multi_window_lock

    properties
        reduce_displayer                reduct_displayer
    end

    methods
        function app = reduct_displayer_sig_display(varargin)
            % collect & save connection to reduce_displayer
            
            % find reduce_displayer
            init_time_pos = cellfun(@(x) (isa(x,'reduce_displayer'))  && x == "init_time",...
                varargin,'UniformOutput',true);
            if ~any(init_time_pos)
                error('Must have reduce_displayer object')
            end
            reduce_displayer = varargin{init_time_pos};
            varargin{init_time_pos} = [];

            % call super constroctor
            app@multi_window_lock(varargin{:})

            app.reduce_displayer = reduce_displayer;
        end

    end
end