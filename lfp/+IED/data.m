classdef data < handle & matlab.mixin.CustomDisplay
    % Hold IED data (i.e. source signal, detection location & rate) and
    % manage the detection steps.
    % Note that this is a handle class - it is updated diractly from
    % curation_window during curation.
    %
    % Based on getIIS by LH (see +IED/legacy folder)
    % By: LdM
    % Published: 230827
    
    
    properties
        
        %%%% props for obj managment %%%%
        status          % text scalar, stage of analysis: "init", "pre_curation", "during_curation", "curated", "analysed"
        last_mark       % numeric pos scalar, last discharge viewed during curation

        %%%% props that hold data %%%%
        % (extracted & source data & data info)

        pos             % discharges detected idx (in samples)
        sig             % original signal, file2read, numeric vector or memmapfile_row
        fs              % numeric pos scalar, signal sampling freq
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

        %%% props that saved res file location %%%%
        file_loc        % full path to file, if exist (empty otherwise)
    end
    

    methods (Access = public)
        function obj = data(varargin)
            % Constractor, simply fill in the properties with user data.
            % assuming user called it before detection, or called detection
            % directly.
            %   INPUTS:
            %       see IED_data properties.
            %   OUTPUTS:
            %       IED_data obj.
            
            p = inputParser;
            addParameter(p, 'sig', [], @(x) isnumeric(x) || isa(x,'memmap_row') || isfile(x))
            addParameter(p, 'fs', 1250, @isnumeric)
            addParameter(p, 'thr', [10 0], @isnumeric)
            addParameter(p, 'thrDir', 'positive', @(x) ismember(x,{'positive','negative','both'}))
            addParameter(p, 'binsize', [], @isnumeric)
            addParameter(p, 'marg', 0.05, @isnumeric)
            addParameter(p, 'file_loc', "", @isstr)
            p.parse(varargin{:})
            
            % fill initial values
            obj.status = "init";

            % fill user values
            obj.sig     = p.Results.sig;
            obj.fs      = p.Results.fs;
            obj.thr     = p.Results.thr;
            obj.thrDir  = p.Results.thrDir;
            if isempty(p.Results.binsize)
                % default bin size for 1 min
                obj.binsize = obj.fs*60;
            else
                obj.binsize = p.Results.binsize;
            end
            obj.marg    = p.Results.marg;

            
        end
        
        function [clipped_discharges, tstamps] = extract_discharges(obj)
            % extract IED discharges according to requested margins around.
            % Only extract accepted discharges!
            % them.
            %   OUTPUT:
            %       tstamps            - time in [ms] for each sample in clipped_discharges
            %       clipped_discharges - discharge X sample, voltage response

            margs = floor(obj.marg * obj.fs); % margs [samples]; obj.marg [ms]
            tstamps = linspace(-obj.marg, obj.marg, margs * 2 + 1); % in [ms]
            true_pos = obj.pos(obj.accepted);
            for iDischarges = numel(true_pos):-1:1
                area2clip = (true_pos(iDischarges)-margs) : (true_pos(iDischarges) + margs);
                clipped_discharges(iDischarges,:) = obj.sig(area2clip);
            end
        end

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

%         function obj = saveobj(obj)
%             % make sure you close manual curation window before saving the
%             % obj, so it won't be problematic
%             if obj.curation_open
%                 f = warndlg('Close curation window before saving the IED data!','close curation before saving');
%                 waitfor(f)
%             end
%         end
    end

    methods (Access=protected)
        function out = getHeader(obj)
            % just create header to the obj that match the package
            out = matlab.mixin.CustomDisplay.getClassNameForHeader(obj);
            out = replace(out,">data<",">IED.data<");
            out = [out ' with properties:'];
        end
    end
end