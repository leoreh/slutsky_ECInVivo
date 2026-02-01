classdef intervals
    %INTERVALS Class for efficient manipulation of time intervals [Start, Stop].
    %
    %   The intervals class provides a robust set of methods for managing and
    %   manipulating sets of time intervals. It supports consolidation, set
    %   operations (union, intersection, difference), and data restriction.
    %   The class ensures that intervals are always sorted by start time upon
    %   construction.
    %
    %   USAGE
    %       obj = intervals(intervalsMatrix)
    %       obj = intervals(startTimes, stopTimes)
    %
    %   PROPERTIES
    %       start   - (Nx1 double) Vector of interval start times.
    %       stop    - (Nx1 double) Vector of interval stop times.
    %       overlap - (logical) True if any intervals overlap.
    %       count   - (scalar) Number of intervals (Dependent).
    %       ints    - (Nx2 double) Matrix of [Start, Stop] (Dependent).
    %       dur     - (scalar) Total duration of intervals (Dependent).
    %
    %   METHODS
    %       consolidate - Merges overlapping or nearby intervals.
    %       intersect   - Calculates intersection with another interval set.
    %       union       - Calculates union with another interval set.
    %       minus       - Subtracts another interval set.
    %       contains    - Checks if data points lie within intervals.
    %       restrict    - Restricts data to the defined intervals.

    %% ====================================================================
    %  PROPERTIES
    %  ====================================================================
    properties (SetAccess = private)
        start   % (double) Nx1 Vector of Start times
        stop    % (double) Nx1 Vector of Stop times
        overlap % (logical) True if intervals overlap
    end

    properties (Dependent)
        count   % (int) Number of intervals
        ints    % (double) Nx2 Matrix [Start, Stop]
        dur     % (double) Total duration
    end

    %% ====================================================================
    %  CONSTRUCTOR
    %  ====================================================================
    methods
        function obj = intervals(varargin)
            %INTERVALS Constructor.
            %
            %   obj = intervals(matrix)
            %       matrix: Nx2 matrix of [Start, Stop] times.
            %
            %   obj = intervals(start, stop)
            %       start: Nx1 vector of start times.
            %       stop:  Nx1 vector of stop times.

            % -----------------------------------------------------------------
            % 1. Handle Empty Input
            % -----------------------------------------------------------------
            if nargin == 0
                obj.start   = [];
                obj.stop    = [];
                obj.overlap = false;
                return;
            end

            % -----------------------------------------------------------------
            % 2. Parse Input
            % -----------------------------------------------------------------
            if nargin == 1
                arg = varargin{1};
                if isa(arg, 'intervals')
                    obj = arg;
                    return;
                end
                if isempty(arg)
                    obj.start   = [];
                    obj.stop    = [];
                    obj.overlap = false;
                    return;
                end
                
                assert(size(arg, 2) == 2, 'Input must be an Nx2 matrix [Start, Stop].');
                startIn = arg(:,1);
                stopIn  = arg(:,2);
            else
                startIn = varargin{1};
                stopIn  = varargin{2};
            end

            % -----------------------------------------------------------------
            % 3. Validation & Formatting
            % -----------------------------------------------------------------
            startIn = double(startIn(:));
            stopIn  = double(stopIn(:));

            % Ensure Start <= Stop (Swap if necessary)
            if any(startIn > stopIn)
                invalid         = startIn > stopIn;
                tmp             = startIn(invalid);
                startIn(invalid)= stopIn(invalid);
                stopIn(invalid) = tmp;
            end

            % Sort strictly by Start time
            [obj.start, idx] = sort(startIn);
            obj.stop         = stopIn(idx);

            % Detect Overlaps
            % Since sorted, overlap exists if start(i) < stop(i-1)
            if length(obj.start) > 1
                % Note: Strict inequality allows touching intervals [0 10], [10 20] 
                % without setting overlap flag.
                obj.overlap = any(obj.start(2:end) < obj.stop(1:end-1));
            else
                obj.overlap = false;
            end
        end

    %% ====================================================================
    %  GETTERS
    %  ====================================================================
        function n = get.count(obj)
            n = length(obj.start);
        end

        function mat = get.ints(obj)
            mat = [obj.start, obj.stop];
        end

        function d = get.dur(obj)
            d = sum(obj.stop - obj.start);
        end

    %% ====================================================================
    %  SET OPERATIONS (CORE)
    %  ====================================================================
        
        % -----------------------------------------------------------------
        function objOut = consolidate(obj, varargin)
            %CONSOLIDATE Merge overlapping or close intervals.
            %
            %   objOut = obj.consolidate('maxGap', eps, 'strict', false)
            %
            %   INPUTS
            %       maxGap - (double) Merge if gap <= maxGap. Default 0.
            %       strict - (bool) If true, uses strict inequality (gap < maxGap).
            %                Default false.
            %
            %   OUTPUT
            %       objOut - (intervals) A new object with consolidated intervals.

            % Parse Inputs
            p = inputParser;
            addParameter(p, 'maxGap', eps, @isnumeric);
            addParameter(p, 'strict', false, @(x) islogical(x) || ischar(x));
            parse(p, varargin{:});
            
            maxGap = p.Results.maxGap;
            strict = p.Results.strict;
            if ischar(strict), strict = strcmp(strict, 'on'); end

            if obj.count <= 1
                objOut = obj;
                return;
            end

            startVec = obj.start;
            stopVec  = obj.stop;

            % -----------------------------------------------------------------
            % Vectorized Merge Logic
            % -----------------------------------------------------------------
            % Strategy:
            % 1. Compute the cumulative maximum of stops (maxPre). This 
            %    resolves complex nested overlaps by efficiently determining 
            %    the furthest end-point reached by any preceding interval.
            % 2. Calculate gaps between current start and previous max-stop.
            % 3. Identify breaks where gap exceeds threshold.
            
            maxPre = cummax(stopVec);
            gaps   = startVec(2:end) - maxPre(1:end-1);

            if strict
                isBreak = gaps >= maxGap;
            else
                isBreak = gaps > maxGap;
            end

            % Define Start/Stop indices for merged groups
            isStart = [true; isBreak];
            isStop  = [isBreak; true];

            % Extract
            newStart = startVec(isStart);
            newStop  = maxPre(isStop);

            objOut = intervals(newStart, newStop);
        end

        % -----------------------------------------------------------------
        function objOut = intersect(obj, query)
            %INTERSECT Intersection of two interval sets (A & B).
            %
            %   objOut = obj.intersect(query)
            %
            %   Computes independent intersection. Consolidates both inputs first.
            
            if ~isa(query, 'intervals'), query = intervals(query); end
            
            % Consolidate inputs to simplify logic
            A = obj.consolidate();
            B = query.consolidate();
            
            startA = A.start; stopA = A.stop;
            startB = B.start; stopB = B.stop;
            nA = length(startA); nB = length(startB);

            % -----------------------------------------------------------------
            % Linear Sweep O(N + M)
            % -----------------------------------------------------------------
            iInt = 1; iQ = 1;
            outStart = zeros(nA + nB, 1); % Preallocate with upper bound
            outStop  = zeros(nA + nB, 1);
            cnt = 0;

            while iInt <= nA && iQ <= nB
                % Intersection range
                start_ = max(startA(iInt), startB(iQ));
                stop_  = min(stopA(iInt), stopB(iQ));

                % Check validity (keep points if start == stop)
                if start_ <= stop_
                     cnt = cnt + 1;
                     outStart(cnt) = start_;
                     outStop(cnt)  = stop_;
                end
                
                % Advance the interval that ends first
                if stopA(iInt) < stopB(iQ)
                    iInt = iInt + 1;
                else
                    iQ = iQ + 1;
                end
            end
            
            % Trim and Create
            objOut = intervals(outStart(1:cnt), outStop(1:cnt));
        end

        % -----------------------------------------------------------------
        function objOut = union(obj, query)
            %UNION Union of two interval sets (A | B).
            %
            %   objOut = obj.union(query)
            
            if ~isa(query, 'intervals'), query = intervals(query); end
            
            combStart = [obj.start; query.start];
            combStop  = [obj.stop; query.stop];
            
            % Simply concatenate and consolidate
            temp = intervals(combStart, combStop);
            objOut = temp.consolidate();
        end

        % -----------------------------------------------------------------
        function [objOut, indices] = minus(obj, query, varargin)
            %MINUS Subtract intervals (A - B).
            %
            %   [objOut, indices] = obj.minus(query)
            %
            %   OUTPUT
            %       objOut  - Resulting intervals.
            %       indices - Original indices from A that survived (partial or full).
            
            if ~isa(query, 'intervals'), query = intervals(query); end
            
            % We need B Consolidated for efficient subtraction
            query = query.consolidate();
            
            startA = obj.start; stopA = obj.stop;
            startB = query.start; stopB = query.stop;
            
            nA = length(startA);
            nB = length(startB);
            
            if nB == 0
                objOut = obj;
                indices = (1:nA)';
                return;
            end
            
            % -----------------------------------------------------------------
            % Sweep Algorithm
            % -----------------------------------------------------------------
            outStart = nan(nA*2, 1);
            outStop  = nan(nA*2, 1);
            outIdx   = nan(nA*2, 1);
            cnt = 0;
            
            iQ = 1; 
            
            for iInt = 1:nA
                currStart = startA(iInt);
                currStop  = stopA(iInt);
                
                % Fast forward Query to relevant range
                while iQ <= nB && stopB(iQ) < currStart
                    iQ = iQ + 1;
                end
                
                % Check overlaps and subtract
                tempQ = iQ;
                while tempQ <= nB && startB(tempQ) < currStop
                    
                    % 1. Part of A before the overlap
                    if startB(tempQ) > currStart
                         cnt = cnt + 1;
                         % Grow buffer if needed
                         if cnt > length(outStart)
                             outStart = [outStart; nan(nA,1)]; %#ok<AGROW>
                             outStop  = [outStop; nan(nA,1)];  %#ok<AGROW>
                             outIdx   = [outIdx; nan(nA,1)];   %#ok<AGROW>
                         end
                         outStart(cnt) = currStart;
                         outStop(cnt)  = startB(tempQ);
                         outIdx(cnt)   = iInt;
                    end
                    
                    % 2. Advance Start of A (effectively cutting the overlap)
                    currStart = max(currStart, stopB(tempQ));
                    
                    tempQ = tempQ + 1;
                end
                
                % 3. Remaining part of A
                if currStart < currStop
                    cnt = cnt + 1;
                     if cnt > length(outStart)
                         outStart = [outStart; nan(nA,1)]; %#ok<AGROW>
                         outStop  = [outStop; nan(nA,1)];  %#ok<AGROW>
                         outIdx   = [outIdx; nan(nA,1)];   %#ok<AGROW>
                     end
                    outStart(cnt) = currStart;
                    outStop(cnt)  = currStop;
                    outIdx(cnt)   = iInt;
                end
            end
            
            % Trim Results
            if cnt > 0
                objOut = intervals(outStart(1:cnt), outStop(1:cnt));
                indices = outIdx(1:cnt);
            else
                objOut = intervals();
                indices = [];
            end
        end

    %% ====================================================================
    %  DATA OPERATIONS
    %  ====================================================================
        
        % -----------------------------------------------------------------
        function [inInts, intIdx] = contains(obj, points)
            %CONTAINS Check if points lie within intervals.
            %
            %   [inInts, intIdx] = obj.contains(points)
            %
            %   OUTPUT
            %       inInts - (logical) Nx1 Logic mask. True if point is inside.
            %       intIdx - (double) Nx2 Matrix. List of [PointIndex, IntervalIndex].
            
            points = double(points);
            nP = length(points);
            inInts = false(nP, 1);
            intIdx = zeros(0, 2);
            
            if obj.count == 0, return; end
            
            if obj.overlap
                % -------------------------------------------------------------
                % Path 1: Sweep (Supports Overlaps)
                % -------------------------------------------------------------
                n = obj.count;
                [sortedStart, sortIdx] = sort(obj.start);
                sortedStop = obj.stop(sortIdx);
                origIds = sortIdx;
                
                [sortedP, pSortIdx] = sort(points(:));
             
                % Preallocation
                estSize = nP; 
                outP = zeros(estSize, 1);
                outI = zeros(estSize, 1);
                cntOut = 0;
                
                tempMask = false(nP, 1);
                active = []; % Active Interval Indices
                iInt = 1;
                
                for iP = 1:nP
                    p = sortedP(iP);
                    origPIdx = pSortIdx(iP);
                    
                    % Add intervals that start before p
                    while iInt <= n && sortedStart(iInt) <= p
                        active(end+1) = iInt; %#ok<AGROW>
                        iInt = iInt + 1;
                    end
                    
                    if ~isempty(active)
                        % Prune expired intervals
                        validMask = sortedStop(active) >= p;
                        active = active(validMask);
                        
                        if ~isempty(active)
                            tempMask(iP) = true;
                            if nargout > 1
                                nNew = length(active);
                                if cntOut + nNew > length(outP)
                                    newSize = max(length(outP)*2, cntOut + nNew + 1000);
                                    outP(newSize, 1) = 0;
                                    outI(newSize, 1) = 0;
                                end
                                outP(cntOut+1 : cntOut+nNew) = origPIdx;
                                outI(cntOut+1 : cntOut+nNew) = origIds(active);
                                cntOut = cntOut + nNew;
                            end
                        end
                    end
                end
                
                inInts(pSortIdx) = tempMask;
                if nargout > 1
                    intIdx = [outP(1:cntOut), outI(1:cntOut)];
                    [~, sortRes] = sort(intIdx(:,1));
                    intIdx = intIdx(sortRes, :);
                end
                
            else
                % -------------------------------------------------------------
                % Path 2: Discretize (Fast, Requires No Overlaps)
                % -------------------------------------------------------------
                bins = discretize(points, [obj.start; inf]);
                candIdx = bins;
                
                valid = ~isnan(candIdx);
                if ~any(valid), return; end
                
                cIdx = candIdx(valid);
                pIndices = find(valid);
                pValid = points(valid);

                % Verify Stop Times
                isIn = pValid <= obj.stop(cIdx);
                
                finalP = pIndices(isIn);
                finalI = cIdx(isIn);
                
                inInts(finalP) = true;
                if nargout > 1
                    intIdx = [finalP, finalI];
                end
            end
        end

        % -----------------------------------------------------------------
        function [dataOut, shifts] = restrict(obj, data, varargin)
            %RESTRICT Restrict data to intervals.
            %
            %   [dataOut, shifts] = obj.restrict(data, 'flgShift', true)
            %
            %   INPUTS
            %       data     - (NxM double) Matrix where Col 1 is time.
            %       flgShift - (bool) If true, shifts segments to eliminate gaps.
            %                  Default false.
            %
            %   OUTPUTS
            %       dataOut  - (double) Filtered and shifted data.
            %       shifts   - (Nx4 double) [OldStart, OldStop, NewStart, NewStop]
            
            p = inputParser;
            addParameter(p, 'flgShift', false, @(x) islogical(x) || x==0 || x==1);
            parse(p, varargin{:});
            flgShift = p.Results.flgShift;
            
            shifts = [];
            if isempty(data), dataOut = []; return; end
            
            % 1. Filter Data
            mask = obj.contains(data(:,1));
            
            if ~any(mask)
                dataOut = [];
                return;
            end
            
            dataOut = data(mask, :);
            
            if ~flgShift
                return;
            end
            
            % 2. Compute Shifts (Time Dilation)
            % Consolidate to ensure clean timeline mapping
            cons = obj.consolidate();
            
            % Map filtered points to consolidated intervals
            [~, adjCons] = cons.contains(dataOut(:,1));
            
            % Sort by point index
            [~, sortOrd] = sort(adjCons(:,1));
            intIds = adjCons(sortOrd, 2);
            
            % Calculate New Starts
            durs = cons.stop - cons.start;
            newStarts = [0; cumsum(durs(1:end-1))];
            
            oldStarts = cons.start(intIds);
            targetStarts = newStarts(intIds);
            
            % Apply Shift
            offset = dataOut(:,1) - oldStarts;
            dataOut(:,1) = targetStarts + offset;
            
            shifts = [cons.start, cons.stop, newStarts, newStarts + durs];
        end

    end

    %% ====================================================================
    %  OVERLOADS (OPERATORS)
    %  ====================================================================
    methods
        function obj = and(obj, other); obj = intersect(obj, other); end
        function obj = or(obj, other);  obj = union(obj, other);     end
    end
end
