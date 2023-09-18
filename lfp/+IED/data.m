classdef data < handle & matlab.mixin.CustomDisplay
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
            % Constructor, simply fill in the properties with user data.
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
            addParameter(p, 'marg', [], @isnumeric)
            addParameter(p, 'smf', [], @isnumeric)
            addParameter(p, 'file_loc', "", @isstr)
            p.parse(varargin{:})
            
            % fill initial values
            obj.status = "init";

            % fill user values
            obj.sig     = p.Results.sig;
            obj.fs      = p.Results.fs;
            obj.thr     = p.Results.thr;
            obj.thrDir  = p.Results.thrDir;
            obj.binsize = p.Results.binsize;
            obj.marg    = p.Results.marg;
            obj.smf     = p.Results.smf;
            
        end
        
        function [clipped_discharges, tstamps] = extract_discharges(obj,time_marg,accepted)
            % extract IED discharges according to requested margins around.
            % Only extract accepted discharges!
            % them.
            %
            %   INPUTS:
            %       obj                - ied object calling this method.
            %       time_marg          - positive scalar, requested time window around each IED in [S].
            %                            If empty or not given, use marg in ied object (default).
            %       accepted           - logical vector, same size of ied.pos.
            %                            If not given, use accpepted in ied object (default).
            %   OUTPUT:
            %       tstamps            - time in [ms] for each sample in clipped_discharges
            %       clipped_discharges - discharge X sample, voltage response
            
            % use object margs if nothing else specefied
            if ~exist("time_marg","var") || isempty(time_marg)
                time_marg = obj.marg;
            end
            if ~exist("accepted","var")
                true_pos = obj.pos(obj.accepted);
            else
                true_pos = obj.pos(accepted);
            end
            % convert from time margings to samples
            margs = floor(time_marg * obj.fs); % margs [samples]; obj.marg [ms]
            
            % create time stamps to match margings
            tstamps = linspace(-time_marg, time_marg, margs * 2 + 1); % in [ms]
            
            % extract waveforms
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
        
        function [accepeted_vec] = filter_accpeted(obj, accepted_filt, accepeted_vec)
            % Match a logical vector that was calc only on accepted IEDs,
            % to a full accpeted vector. Used to fix accpeted after you
            % Used for secound - stage, to remove IEDs that were accepted
            % earlier.
            % Note that order is important - values in accepted_filt should
            % correspond with find(accpeted_vec).
            %
            % INPUT:
            %   obj - ied object calling this function.
            %   accepted_filt - logical vector, numel(accepted_filt) == sum(accpeted_vec).
            %                   true for values that should remine accpeted.
            %   accpeted_vec  - logical vector, numel(accpeted_vec) == numel(ied.pos).
            %                   what positions of IED are accepted originaly.
            %                   if not given, use ied.accepted.
            %
            % OUTPUT:
            %   accpeted_vec  - same as accpeted_vec input, but only accpeted
            %                   values that were true in accepted_filt are
            %                   still true.
            %
            % EXAMPLE:
            %   %%% General: %%%
            %   accpeted_vec = [0 1 1 0];
            %   accepted_filt = [0 1];
            %   new_accepted = filter_accpeted(obj, accepted_filt, accpeted_vec)
            %
            %   %new_accepted = [0 0 1 0]
            %   
            %   %%% Assign new accpeted %%%
            %   ied.accpeted = [0 1 1 0];
            %   ied.accpeted = filter_accpeted(obj, accepted_filt)
            
            if ~exist("accepeted_vec","var")
                accepeted_vec = obj.accepted;
            elseif numel(accepeted_vec) ~= numel(obj.pos) || ~islogical(accepeted_vec) || ~isvector(accepeted_vec)
                error('accpeted_vec need to be a logical vector, & match size with ied.pos')
            end
            
            if numel(accepted_filt) ~= sum(accepeted_vec)
                error("accepted_filt must have a value for each true val in accpeted_vec, total %d",sum(accepeted_vec))
            end
            
            true_pos = find(accepeted_vec);
            pos2rmv = true_pos(~accepted_filt);
            accepeted_vec(pos2rmv) = false; 

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

        function out = getFooter(obj)
            detected_events = numel(obj.pos);
            accepted_events = sum(obj.accepted);
            out = sprintf('Status: %s, Detected: %d, Accepted: %d',...
                obj.status,detected_events,accepted_events);
        end
    end
end

% EOF