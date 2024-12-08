classdef IED_data2disp < data2disp_model

    properties (Dependent)
        data_sources
    end

    properties
        ied IED.data
    end
    
    properties (Hidden, Transient)
        thr_vals = []
        temp_event_types (:,1) string = string.empty()
    end

    events
        detection_changed
    end

    properties(Hidden)
        init_accepted % accepted status for detection before any curation
        use_only_init_accepted (1,1) logical = false
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% Deal With Loading New Data %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function add_data(obj, options)
            arguments
                obj

                options.ied
                options.init_time
                
                options.only_init_accepted (1,1) logical = false

                options.temp_events
            end

            % push in ied object
            obj.ied = options.ied;

            if options.only_init_accepted
                obj.init_accepted = obj.ied.accepted;
                obj.use_only_init_accepted = true;
            end
            
            % use default marg if needed
            if isempty(obj.ied.marg)
                obj.ied.marg = 0.15;
                fprintf("setting ied.marg to be %g\n", obj.ied.marg)
            end

            % write down the detected events
            if ~isempty(obj.ied.pos)
                obj.events_time{"Detected"} = ([-obj.ied.marg obj.ied.marg]./2) + ([obj.ied.pos obj.ied.pos]./obj.ied.fs);
            else
                obj.events_time{"Detected"} = double.empty(0,2);
            end
            % insert the manual marked events in the ied
            if ~isempty(obj.ied.manual_marked_events)
                mm_events = obj.ied.manual_marked_events;
                obj.events_time = obj.events_time.insert(string(mm_events.keys), mm_events.values);
            end
            % insert events2display which should not be kept in the ied object
            if isfield(options, "temp_events")
                temp_events = obj.validate_events(options.temp_events);
                obj.events_time = obj.events_time.insert(string(temp_events.keys), temp_events.values);
                obj.temp_event_types = string(temp_events.keys);
            end

            % remove default event type
            if isempty(obj.ied.manual_marked_events) || ~isKey(obj.ied.manual_marked_events,"Event")
                obj.rmv_event_type("Event");
            end

            % init time 2 display
            obj.move_time(options.init_time, [], false)

            % finish loading, notify about it by order:
            notify(obj, "data_loaded")
            notify(obj, "events_types_changed")
            % notify(obj, "events_changed")
            notify(obj, "in_win_change") % note that this include any "events_changed" responce
            notify(obj,"detection_changed")
        end
    
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% Deal With Events %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function rmv_event(obj, time_point, event_type)
            if event_type == "Detected"
                error("IED_data2disp:rmv_event:no_removing_detected",...
                    "Can't remove detected events")
            else
                rmv_event@data2disp_model(obj, time_point, event_type);
                obj.match_manual_events_ied_obj();
            end
        end
    
        function add_event(obj, event_t, event_type)
            if event_type == "Detected"
                error("IED_data2disp:add_event:no_adding_detected",...
                    "Can't add to detected events, create a new event type instead")
            else
                add_event@data2disp_model(obj, event_t, event_type);
                obj.match_manual_events_ied_obj();
            end
        end
        
        function rmv_event_type(obj, event_type)
            if event_type == "Detected"
                error("IED_data2disp:rmv_event_type:no_removing_detected",...
                    "Can't remove detected events")
            else
                rmv_event_type@data2disp_model(obj, event_type);
                obj.match_manual_events_ied_obj();
            end
        end

        function add_event_type(obj, event_type)
            if event_type == "Detected"
                error("IED_data2disp:add_event_type:no_adding_detected",...
                    "detected is a reserved keyword for events detected by the IED algoritem")
            else
                add_event_type@data2disp_model(obj, event_type);
                obj.match_manual_events_ied_obj();
            end
        end
    
    end
    
    methods (Hidden)
        function match_manual_events_ied_obj(obj)
            obj.ied.manual_marked_events = obj.events_time.remove(["Detected"; obj.temp_event_types]);
        end

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%% Deal With IED Detections %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function change_detection_status(obj, new_status)
            % change current detection (ied.last_mark) and move to next
            obj.ied.accepted(obj.ied.last_mark) = new_status;
        end
        
        function move_IED_detection(obj, detect2move, only_init_accepted)
            % move to a specific IED detection.
            %
            % INPUTS:
            %   detect2move - number of event to move to, or
            %       "Next"/"Previous" to move relativly to ied.last_mark.
            % note that this function can move window to be outside data
            % limits.

            if ~exist("only_init_accepted","var") || isempty(only_init_accepted)
                only_init_accepted = obj.use_only_init_accepted;
            end

            % let user ask to move in a relative way
            if strcmp(detect2move, "Next") % will return false if detect2move is numeric
                if only_init_accepted
                    detect2move = obj.ied.last_mark + find(obj.init_accepted((obj.ied.last_mark+1):end),1);
                else
                    detect2move = obj.ied.last_mark + 1;
                end

                if isempty(detect2move) || detect2move > numel(obj.ied.pos)
                    error('IED_data2disp:move_IED_detection:out_of_bounds','Cant Move after last detection')
                end

            elseif strcmp(detect2move, "Previous")
                if only_init_accepted
                    detect2move = find(obj.init_accepted(1:(obj.ied.last_mark-1)),1,"last");
                else
                    detect2move = obj.ied.last_mark -1;
                end

                if isempty(detect2move) || detect2move < 1
                    error('IED_data2disp:move_IED_detection:out_of_bounds','Cant move before first detection')
                end
            end

            % find time of detection and move there
            detection_time = obj.collect_detection(detect2move);

            % move time
            obj.move_time(detection_time)

            % write down new marking location
            obj.ied.last_mark = detect2move;

            % report detection changed
            notify(obj,"detection_changed")
        end
    
        function thr_vals = collect_threshold(obj, win_limits)
            % collect the threshold values matching the signal that is
            % under use.
            %
            % INPUT:
            % win_limits - 
            %       numeric 2 element, start-stop values to collect
            %       threshold in. Must be sorted.

            
            
            if isempty(obj.thr_vals)
                if numel(obj.ied.thr) == 2
                    % assume that threshold is global - therefore, only need to
                    % replicate its mV value to match data.
                    obj.thr_vals = repelem(obj.ied.thr(2), length(obj.ied.sig));
                elseif numel(obj.ied.thr) == 1
                    % assume local z-score, over window of 5 sec.
                    sig_mu = movmean(obj.ied.sig, obj.ied.fs*5);
                    sig_sigma = movstd(obj.ied.sig, obj.ied.fs*5);
                    obj.thr_vals = obj.ied.thr.*sig_sigma + sig_mu;
                end
            end

            dataset_ind = obj.time2ind(obj.ied.sig2use, win_limits);
            thr_vals = obj.thr_vals(dataset_ind);
        end
    
        function [detection_time, detection_val] = collect_detection(obj, detection2collect)
            if ~exist("detection2collect","var") || isempty(detection2collect)
                detection2collect = obj.ied.last_mark;
            end

            % find time & value of detection
            detection_samp = obj.ied.pos(detection2collect);
            detection_time = detection_samp./obj.ied.fs;
            detection_val = obj.ied.sig(detection_samp);

        end
    
        function detect_samp = find_closest_detection(obj, time2approx, only_init_accepted)
            % find the closest sample that is in ied.pos to a given time.
            
            if ~exist("only_init_accepted","var") || isempty(only_init_accepted)
                only_init_accepted = obj.use_only_init_accepted;
            end

            % collect which samples to compare
            if only_init_accepted
                pos2comp = obj.ied.pos(obj.ied.accepted);
            else
                pos2comp = obj.ied.pos;
            end

            % find closest match in time
            t2comp = pos2comp./obj.ied.fs;
            [~,close_detect] = min(abs(time2approx-t2comp));
            detect_samp = pos2comp(close_detect);
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% Get Set Methods %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function res = get.data_sources(obj)
            % collect from ied, dup sig 2 display
            res = obj.ied.data_sources;
            res(obj.ied.sig2use + " "+ "Zoomed Out",:) = res(obj.ied.sig2use,:);
        end
        
    end

end