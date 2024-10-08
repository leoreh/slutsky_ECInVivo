classdef data < handle & matlab.mixin.CustomDisplay & matlab.mixin.Copyable
    % Hold IED data (i.e. source signal, detection location & rate) and
    % manage the detection steps.
    % Note that this is a handle class - it is updated directly from
    % curation_window during curation.
    %
    % Based on getIIS by LH (see +IED/legacy folder)
    % By: LdM
    % Published: 230827


    properties

        %%%% props for obj managment %%%%

        status          % text scalar, stage of analysis: "init", "pre_curation", "during_curation", "curated", "analyzed"
        last_mark       % numeric pos scalar, last discharge viewed during curation

        %%%% props that hold data %%%%
        % (extracted & source data & data info)
        
        % input data
        data_sources    % table, with variables "data" (cell) & "fs" (double),
                        % row names are data labels. 
        sig2use         % label of the channel to perform all the calculatons on
        
        % analysis results
        pos             % discharges detected idx (in samples)
        accepted        % logical vec, true when discharge is accepted (default: true)
        rate            % calculated discharge rate, see times2rate
        edges           % edges of bins used for rate, see times2rate
        cents           % centers of bins used for rate, see times2rate
        manual_marked_events % cell_dict, events that were marked manually during curation

        %%%% props that hold analysis parameters %%%
        
        % thr - thr used for detection, either 2 elem vector [z-score mV] or 1 elem [z-score for moving z socre threshold]
        thr double {mustBeNonnegative}           

        % thrDir - text scalar, direction of thr, "positive", "negative" or "both"
        thrDir (1,1) string     {mustBeMember(thrDir, {'positive','negative','both'})} = "both"

        % marg - instruction of how much to clip around each discharge, in [s]
        marg double {mustBePositive} = []

        % binsize - numeric scalar, size of bin for rate calc in [s]
        binsize double {mustBePositive} = []

        % smf - smoothing window for rate, in [bins]
        smf double {mustBePositive, mustBeInteger} = []

        %%% props that give metadata about results %%%%

        file_loc string % full path to file, if exist (empty otherwise)
        git_created     % slutsky_ECInVivo & cellexplorer versions upon creation
        git_last_step   % slutsky_ECInVivo & cellexplorer versions upon last analysis step
    end

    properties (Dependent)
        sig % quick handle to the signal to perform calculations
        fs  % quick handle to the sampling frequency of the same 
    end
    

    %%%%%%%%%%%%%%%% Constructor  %%%%%%%%%%%%%%%%%%%%
    methods (Access = public)
        function obj = data(data, fs, label, options)
            % Constructor, simply fill in the properties with user data.
            % assuming user called it before detection, or called detection
            % directly.
            %   INPUTS:
            %       see IED_data properties.
            %   OUTPUTS:
            %       IED_data obj.

            arguments (Repeating)
                data    {mustBeA(data, ["int8", "int16", "int32", "int64", "uint8", "uint16", "uint32", "uint64",...
                    "single", "double", "memmap_row"])}
                fs      (1,1) double {mustBePositive}
                label   (1,1) string {mustBeTextScalar}
            end

            arguments
                options.sig2use

                options.thr     double  {mustBeNonnegative}
                options.thrDir  string  {mustBeMember(options.thrDir, {'positive','negative','both'})}
                options.binsize double  {mustBePositive}
                options.marg    double  {mustBePositive} 
                options.smf     double  {mustBePositive} 

                options.file_loc string
            end

            % fill initial values
            obj.status = "init";

            % fill user values
            % for iField = fieldnames(p.Results)'
            for iField = fieldnames(options)'
                obj.(iField{:}) = options.(iField{:});
            end
            
            % fill data values
            data_sources = table('Size',[0 2],'VariableTypes',["cell", "double"], 'VariableNames', ["data", "fs"]);
            warning('off','MATLAB:table:RowsAddedExistingVars')
            for iData = 1:numel(data)
                data_sources.data{label{iData}} = data{iData};
                data_sources.fs(label{iData}) = fs{iData};
            end
            obj.data_sources = data_sources;

            % add git info
            try
                obj.git_created = get_gits_status(["slutsky_ECInVivo", "CellExplorer"]);
            catch err
                obj.git_created = join(["Error: " err.message]);
            end
        end

    end

    %%%%%%%%%%%%%%%% Class Manipulations %%%%%%%%%%%%%
    methods (Access = public)

        %%%%%%% Pointers to Package Functions
        function obj = detect(obj, varargin)
            % just a cover to point for package function IED.detect

            obj = IED.detect(obj,varargin{:});
        end

        function obj = curate(obj, varargin)
            % just a cover to point for package function IED.detect

            obj = IED.curate(obj,varargin{:});
        end

        function obj = analyze(obj, varargin)
            % just a cover to point for package function IED.analyse

            obj = IED.analyze(obj,varargin{:});
        end
    
    end
    
    %%%%%%%%%%% Collecting Data From Class %%%%%%%%%%%
    methods (Access = public)

        [accepeted_vec] = filter_accpeted(obj, accepted_filt, accepeted_vec)

        [pos_labels] = match_event_state(obj, labels, accepted)
   
        [clipped_discharges, tstamps] = extract_discharges(obj,time_marg,accepted)

        [events_width_half_amp, events_amp, event_width_10_amp] = calculate_props(obj, margs)
    end
    
    %%%%%%%%%%%%%% Modify Class Display %%%%%%%%%%%%%%
    methods (Access=protected)
        function out = getHeader(obj)
            % just create header to the obj that match the package
            out = matlab.mixin.CustomDisplay.getClassNameForHeader(obj);
            out = replace(out,">data<",">IED.data<");
            out = [out ' with properties:'];
        end

        function out = getFooter(obj)
            if numel(obj) > 1
                out = '';
                return
            end
            detected_events = numel(obj.pos);
            accepted_events = sum(obj.accepted);
            out = sprintf('Status: %s, Detected: %d, Accepted: %d\n',...
                obj.status,detected_events,accepted_events);
        end
    end

    %%%%%%%%%%%%%% Save & Load %%%%%%%%%%%%%%
    methods (Static)
        function obj = loadobj(s)
            % deal with errors during loading process
            if isstruct(s)
                if ~isfield(s,"sig2use")
                    s.sig2use = "LFP";
                end

                if ~isfield(s, "data_sources")
                    % match legacy version to new one
                    
                    % collect data-repeating values
                    data{1} = s.sig;
                    data{2} = s.fs;
                    data{3} = "LFP";
                    if isfield(s,"emg") && ~isempty(s.emg)
                        data{4} = s.emg;
                        data{5} = s.fs;
                        data{6} = "EMG";
                    end
                else
                    data = {};
                end
                
                if isempty(s.thrDir)
                    s.thrDir = "both";
                end
                % collect name-value arguments
                all_fields = fieldnames(s);
                fields2use = ["sig2use", "thr", "thrDir", "binsize", "marg", "smf", "file_loc"];
                fields2remove = all_fields(~ismember(all_fields, fields2use));
                in_struct = rmfield(s, fields2remove);
                name_val_arg = namedargs2cell(in_struct);
                
                % create object
                obj = IED.data(data{:}, name_val_arg{:});

                %%%%% correct parameters that are outside of obj creation
                
                % git data
                if isfield(s, "git_created")
                    obj.git_created = s.git_created;
                end
                try
                    obj.git_last_step = get_gits_status(["slutsky_ECInVivo", "CellExplorer"]);
                catch err
                    obj.git_last_step = join(["Error: " err.message]);
                end
                
                % object work status
                obj.status = s.status;
                obj.last_mark = s.last_mark;
                
                % obj detection & curation
                obj.pos = s.pos;
                obj.accepted = s.accepted;
                
                % obj analysis results
                obj.rate = s.rate;
                obj.edges = s.edges;
                obj.cents = s.cents;
            else
                % no error - just move on
                obj = s;
            end
        end

    end

    %%%%%%%%%%%%%% Set & Get props %%%%%%%%%%%%%%
    methods
        function res = get.sig(obj)
            % save legacy option of sig reference - jump to actual signal 2 use
            res = obj.data_sources.data{obj.sig2use};

        end
        
        function set.sig(obj, value)
            % pass the value to signal 2 use
            obj.data_sources.data{obj.sig2use} = value;
        end

        function res = get.fs(obj)
            % save legacy option of sig reference - jump to signal 2 use fs
            res = obj.data_sources.fs(obj.sig2use);

        end

        function set.fs(obj, value)
            % pass the value to signal 2 use
            obj.data_sources.fs(obj.sig2use) = value;
        end
    end

end

% EOF