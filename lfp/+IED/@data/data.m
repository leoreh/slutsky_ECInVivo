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
        sig             % original signal, numeric vector or memmapfile_row
        emg             % emg signal, numeric vector or memmapfile_row.
                        %   MUST be with the same fs as sig, and of the
                        %   same length.
        fs              % numeric pos scalar, signal sampling freq
        
        % analysis results
        pos             % discharges detected idx (in samples)
        accepted        % logical vec, true when discharge is accepted (default: true)
        rate            % calculated discharge rate, see times2rate
        edges           % edges of bins used for rate, see times2rate
        cents           % centers of bins used for rate, see times2rate

        %%%% props that hold analysis parameters %%%

        thr             % 2 elem vector, thr used for detection, [z-score mV]
        thrDir          % text scalar, direction of thr, "positive", "negative" or "both"
        marg            % instruction of how much to clip around each discharge, in [s]
        binsize         % numeric scalar, size of bin for rate calc in [s]
        smf             % smoothing window for rate, in [bins]

        %%% props that give metadata about results %%%%

        file_loc        % full path to file, if exist (empty otherwise)
        git_created     % slutsky_ECInVivo & cellexplorer versions upon creation
        git_last_step   % slutsky_ECInVivo & cellexplorer versions upon last analysis step
    end
    

    %%%%%%%%%%%%%%%% Constructor  %%%%%%%%%%%%%%%%%%%%
    methods (Access = public)
        function obj = data(varargin)
            % Constructor, simply fill in the properties with user data.
            % assuming user called it before detection, or called detection
            % directly.
            %   INPUTS:
            %       see IED_data properties.
            %   OUTPUTS:
            %       IED_data obj.

            p = inputParser;
            % input data
            addParameter(p, 'sig', [], @(x) isnumeric(x) || isa(x,'memmap_row'))
            addParameter(p, 'emg', [], @(x) isnumeric(x) || isa(x,'memmap_row'))
            addParameter(p, 'fs',  [], @isnumeric)

            % analysis parms
            addParameter(p, 'thr', [], @isnumeric)
            addParameter(p, 'thrDir', 'both', @(x) ismember(x,{'positive','negative','both'}))
            addParameter(p, 'binsize', [], @isnumeric)
            addParameter(p, 'marg', [], @isnumeric)
            addParameter(p, 'smf', [], @isnumeric)

            % metadata
            addParameter(p, 'file_loc', "", @isstr)
            
            p.parse(varargin{:})

            % fill initial values
            obj.status = "init";

            % fill user values
            for iField = fieldnames(p.Results)'
                obj.(iField{:}) = p.Results.(iField{:});
            end

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
end

% EOF