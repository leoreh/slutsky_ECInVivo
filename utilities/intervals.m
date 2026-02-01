classdef intervals
    %INTERVALS Class for efficient manipulation of intervals (start, stop).
    %
    %   USAGE:
    %       obj = intervals(intervals)
    %       obj = intervals(startTimes, endTimes)
    %
    %   METHODS:
    %       consolidate     - Merges overlapping intervals.
    %       intersect       - Intersection of two interval sets.
    %       union           - Union of two interval sets.
    %       minus           - Sutraction (Set difference).
    %       contains        - Checks if points are within intervals.
    %       restrict        - Restricts data to intervals (with shifting).
    %
    % ====================================================================
    %  PROPERTIES
    % ====================================================================
    properties (SetAccess = private)
        start   % (double) Nx1 Vector of Start times
        stop    % (double) Nx1 Vector of Stop times
        overlap % (logical) True if intervals overlap
    end
    
    properties (Dependent)
        count       % (int) Number of intervals
        ints        % (double) Nx2 Matrix [Start, Stop]
        dur         % (double) Total duration
    end
    
    % ====================================================================
    %  CONSTRUCTOR
    % ====================================================================
    methods
        function obj = intervals(varargin)
            %INTERVALS Constructor.
            %
            %   obj = intervals(mat)
            %   obj = intervals(start, stop)
            
            if nargin == 0
                obj.start = [];
                obj.stop = [];
                obj.overlap = false;
                return;
            end
            
            % Handle Input
            if nargin == 1
                arg = varargin{1};
                if isa(arg, 'intervals')
                    obj = arg;
                    return;
                end
                if isempty(arg)
                    obj.start = [];
                    obj.stop = [];
                    obj.overlap = false;
                    return;
                end
                
                % Assume Nx2 matrix
                assert(size(arg, 2) == 2, 'Input must be an Nx2 matrix [Start, Stop].');
                startIn = arg(:,1);
                stopIn = arg(:,2);
            else
                startIn = varargin{1};
                stopIn = varargin{2};
            end
            
            % Enforce Columns
            startIn = double(startIn(:));
            stopIn = double(stopIn(:));
            
            % Validation: Start <= Stop
            if any(startIn > stopIn)
                % Auto-fix swapped
                invalid = startIn > stopIn;
                tmp = startIn(invalid);
                startIn(invalid) = stopIn(invalid);
                stopIn(invalid) = tmp;
            end
            
            % Sort
            % We strictly sort to ensure O(N log N) ops later
            [obj.start, idx] = sort(startIn);
            obj.stop = stopIn(idx);
            
            % Check Overlap
            % Since sorted, overlap exists if start(i) < stop(i-1)
            if length(obj.start) > 1
                % Strict inequality allows touching [0 10] [10 20] without "overlap" flag
                % This is consistent with standard interval logic where touching != overlapping
                % for the purpose of ambiguity.
                obj.overlap = any(obj.start(2:end) < obj.stop(1:end-1));
            else
                obj.overlap = false;
            end
        end
        
        % ================================================================
        %  GETTERS
        % ================================================================
        function n = get.count(obj)
            n = length(obj.start);
        end
                
        function mat = get.ints(obj)
            mat = [obj.start, obj.stop];
        end
                
        function d = get.dur(obj)
            d = sum(obj.stop - obj.start);
        end
        
        % ================================================================
        %  CORE METHODS
        % ================================================================
        
        % -----------------------------------------------------------------
        function objOut = consolidate(obj, varargin)
            %CONSOLIDATE Merge overlapping or close intervals.
            %
            %   objOut = obj.consolidate('maxGap', 0, 'strict', false)
            %
            %   INPUTS:
            %       maxGap  - (double) Merge if gap <= maxGap. 
            %                 Default 0. 
            %                 Acts as a bridge for small gaps.
            %       strict  - (bool) If true, strictly less than is used for gap check.
            %                 If maxGap=0 and strict=true, touching intervals [0,1][1,2] 
            %                 will NOT merge (gap=0 is not < 0).
            %                 Default false (touching intervals merge).
            
            % Parse
            p = inputParser;
            addParameter(p, 'maxGap', 0, @isnumeric);
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
            stopVec = obj.stop;
            
            % Calculate gap between current start and previous stop
            % gap(i) corresponds to gap between interval (i-1) and (i)
            % Since start/stop are aligned: gap between i+1 and i is start(i+1) - stop(i)
            gaps = startVec(2:end) - stopVec(1:end-1);
            
            % Determine Merge Conditions
            if strict
                shouldMerge = gaps < maxGap;
            else
                shouldMerge = gaps <= maxGap;
            end
            
            % If no merges needed, return early (deep copy)
            if ~any(shouldMerge)
                objOut = obj;
                return;
            end
            
            % Vectorized Merge Logic
            % A "Start" of a merged block is one where the previous gap was NOT merged.
            % The first interval is always a start.
            % isStart(i) is true if gap(i-1) was not merged (i.e. gap > maxGap)
            % Note: gaps has length N-1.
            % index i of gaps corresponds to junction between i and i+1.
            
            % isStart(k) corresponds to interval k.
            % interval 1 is always start.
            % interval k (for k>1) is start if merge condition at k-1 (gaps(k-1)) is false.
            isStart = [true; ~shouldMerge];
            
            % A "Stop" of a merged block is one where the NEXT gap does NOT merge.
            % The last interval is always a stop.
            % interval k (for k<N) is stop if merge condition at k (gaps(k)) is false.
            isStop = [~shouldMerge; true];
            
            % Extract
            newStart = startVec(isStart);
            newStop = stopVec(isStop);
            
            % However, simple selection [start(isStart), stop(isStop)] works for simple gaps,
            % but does it work for nested/complex overlaps?
            % Ints class guarantees sorted starts, but does NOT guarantee sorted stops or proper nesting removal yet (overlap=true).
            % If we have [0 10] [2 5], start=[0 2], stop=[10 5].
            % Gap: 2-10 = -8. Merge=True.
            % isStart=[1 0]. isStop=[0 1].
            % newS = 0. newE = 5. -> Expect [0 10].
            % PROBLEM: Vectorized logic assumes Stop(i) extends far enough.
            % If we have nested intervals, we might pick a Stop that is smaller than a previous Stop.
            %
            % We need to resolve "effective max end" for the group.
            
            % If overlaps exist, we must handle the "max(stop)" propagation.
            % Since we are consolidating, we can pre-process stops to ensure they are monotonic-ish?
            % No, standard approach:
            
            if obj.overlap
                % If we have complex overlaps, the simple binary2bouts logic might fail on duration
                % extending backwards.
                % binary2bouts assumes sequential logic works or stops are somewhat aligned?
                % Actually binary2bouts does: iei = b(2:,1) - b(1:-1,2).
                % This relies on stop being the end of the previous BOUT.
                % In binary2bouts, input is sorted and non-overlapping (from binary vector).
                % HERE, inputs can overlap.
                
                % To use vectorized merge on *overlapping* intervals, we first need to 
                % resolve the max stop time for any chain of overlaps.
                % But we can do this simply by accumulating max stop.
                
                % Iterative approach is O(N), but Matlab loops are slow.
                % Cummax needed?
                % If we sort by start (already done), the merge group is contiguous.
                % For a merge group i..j, newStart = start(i), newStop = max(stop(i..j)).
                
                % Let's adjust the stops first to be "running max" within merge groups?
                % Can we determine merge groups?
                % A break in merge group is when gap > maxGap.
                % But gap calculation depends on the "Current Effective End".
                % If I have [0 20] [5 10] [25 30]. maxGap=0.
                % Gap 1: 5-20 = -15 (Merge).
                % Gap 2: 25-10 = 15 (No Merge? Wait. 25 SHOULD be compared to 20, not 10).
                
                % So strict vectorized logic "start(i+1) - stop(i)" fails if stop(i) < stop(i-1).
                % We need the cumulative max stop time relative to the *group start*?
                
                % Optimized loop for Matlab (JIT is good for simple loops).
                % The previous implementation used a loop. It was robust.
                % To make it "efficient but robust", we can keep the loop but optimize variables.
                % Or use `cummax`?
                
                % Let's try `cummax` on stops? 
                % [0 20], [5 10], [25 30] -> stops: [20 10 30]. cummax: [20 20 30].
                % gaps relative to cummax:
                % gap 1: 5 - 20 = -15.
                % gap 2: 25 - 20 = 5. -> Break.
                % This works for finding breaks!
                
                % Algorithm with cummax:
                % 1. maxPre = cummax(stop)
                % 2. gaps = start(2:end) - maxPre(1:end-1)
                % 3. isBreak = gaps > maxGap (or >= if strict)
                % 4. newStart = start([true; isBreak])
                % 5. Determine stops. Stop of a group is max(stop) within that group.
                %    The group ends at index i where isBreak(i) is true (start of next group).
                %    So we need max(stop) at the indices before breaks.
                %    Wait, maxPre IS the running max.
                %    At the end of a group (before a break), maxPre is exactly the max of that group.
                %    So we just need to grab values from maxPre at the end-of-group indices.
                
                % Let's trace [0 20] [5 10]. cummax=[20 20].
                % gap: 5 - 20 = -15. isBreak=False.
                % isBreak = [F].
                % isStart = [T, F]. -> newStart = start(1) = 0.
                % isStop = [F, T]. -> newStop = maxPre(2) = 20.
                % Result [0 20]. Correct.
                
                % Trace [0 10] [15 20]. cummax=[10 20].
                % gap: 15 - 10 = 5. maxGap=0. isBreak=True.
                % isStart = [T, T]. -> newStart = [0, 15].
                % isStop = [T, T]. -> newStop = [maxPre(1), maxPre(2)] = [10, 20].
                % Result [0 10] [15 20]. Correct.
                
                % Trace with strict. [0 10] [10 20]. cummax=[10 20].
                % gap: 10 - 10 = 0.
                % if strict=true, isBreak = (0 >= 0) -> True?
                % Strict logic: gap < maxGap (merge). So gap >= maxGap (break).
                % If maxGap=0, strict=true. 0 >= 0 -> True (Break).
                % Result [0 10] [10 20]. Correct.
                
                % If strict=false. isBreak = gap > maxGap.
                % 0 > 0 -> False (Merge).
                % Result [0 20]. Correct.
                
                % What if [0 10] [2 5] [2 20]. cummax=[10 10 20].
                % gap1: 2-10=-8. Merge.
                % gap2: 2-10=-8. Merge.
                % isBreak=[F F].
                % isStart=[T F F]. Start=0.
                % isStop=[F F T]. Stop=20.
                % Result [0 20]. Correct.
                
                % This looks like a solid robust vectorized algorithm using cummax.
                
                maxPre = cummax(stopVec);
                gaps = startVec(2:end) - maxPre(1:end-1);
                
                if strict
                    isBreak = gaps >= maxGap;
                else
                    isBreak = gaps > maxGap;
                end
                
                isStart = [true; isBreak];
                isStop = [isBreak; true];
                
                newS = startVec(isStart);
                newE = maxPre(isStop);
                
                objOut = intervals(newS, newE);
                
            else
                % If already no overlap, we just need to merge close ones.
                % But we still need to chain merges (A merges B, B merges C -> A merges C).
                % cummax helps here too, essentially same logic.
                % If sorted and no overlap, start(i) > stop(i-1).
                % stop(i) could perform a bridge.
                % cummax(stop) will just be stop (mostly) unless we merge?
                % No, if we have [0 10] [11 20] [12 22] (Wait, non-overlapping means start(i) >= stop(i-1)? 
                % The property `overlap` defines strict overlap. Touching is not "overlap" in current prop logic.
                % But consolidate might merge touching.
                % The cummax logic is universal. I will usage it for both cases.
                
                maxPre = cummax(stopVec);
                gaps = startVec(2:end) - maxPre(1:end-1);
                
                if strict
                    isBreak = gaps >= maxGap;
                else
                    isBreak = gaps > maxGap;
                end
                
                isStart = [true; isBreak];
                isStop = [isBreak; true];
                
                objOut = intervals(startVec(isStart), maxPre(isStop));
            end
        end
        
        % -----------------------------------------------------------------
        function [objOut, indices] = minus(obj, query, varargin)
            %MINUS Subtract intervals (A - B).
            %
            %   objOut = obj.minus(query)
            %
            %   Uses efficient O(N + M) verification on sorted lists.
            
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
            
            % Result Preallocation
            outStart = nan(nA*2, 1);
            outStop = nan(nA*2, 1);
            outIdx = nan(nA*2, 1);
            cnt = 0;
            
            iQ = 1; % Index for Query
            
            for iInt = 1:nA
                currStart = startA(iInt);
                currStop = stopA(iInt);
                
                % Find relevant B intervals
                while iQ <= nB && stopB(iQ) < currStart
                    iQ = iQ + 1;
                end
                
                % Process Overlaps
                tempQ = iQ;
                while tempQ <= nB && startB(tempQ) < currStop
                    % B(tempQ) overlaps with Curr
                    
                    % Part before B
                    if startB(tempQ) > currStart
                         cnt = cnt + 1;
                         % Grow if needed
                         if cnt > length(outStart)
                             outStart = [outStart; nan(nA,1)]; %#ok<AGROW>
                             outStop = [outStop; nan(nA,1)]; %#ok<AGROW>
                             outIdx = [outIdx; nan(nA,1)]; %#ok<AGROW>
                         end
                         outStart(cnt) = currStart;
                         outStop(cnt) = startB(tempQ);
                         outIdx(cnt) = iInt;
                    end
                    
                    % Advance Start
                    currStart = max(currStart, stopB(tempQ));
                    
                    tempQ = tempQ + 1;
                end
                
                % Remaining part
                if currStart < currStop
                    cnt = cnt + 1;
                     if cnt > length(outStart)
                         outStart = [outStart; nan(nA,1)]; %#ok<AGROW>
                         outStop = [outStop; nan(nA,1)]; %#ok<AGROW>
                         outIdx = [outIdx; nan(nA,1)]; %#ok<AGROW>
                     end
                    outStart(cnt) = currStart;
                    outStop(cnt) = currStop;
                    outIdx(cnt) = iInt;
                end
            end
            
            % Trim
            if cnt > 0
                outStart = outStart(1:cnt);
                outStop = outStop(1:cnt);
                indices = outIdx(1:cnt);
                objOut = intervals(outStart, outStop);
            else
                objOut = intervals();
                indices = [];
            end
        end
        
        % -----------------------------------------------------------------
        function objOut = intersect(obj, query)
            %INTERSECT Intersection of two interval sets (A & B).
            
            if ~isa(query, 'intervals'), query = intervals(query); end
            
            A = obj.consolidate();
            B = query.consolidate();
            
            startA = A.start; stopA = A.stop;
            startB = B.start; stopB = B.stop;
            
            nA = length(startA); nB = length(startB);
            
            % Sweep
            iInt = 1; iQ = 1;
            outStart = []; outStop = [];
            
            while iInt <= nA && iQ <= nB
                % Max of starts, Min of ends
                start_ = max(startA(iInt), startB(iQ));
                stop_ = min(stopA(iInt), stopB(iQ));
                
                if start_ < stop_ % Using < to avoid zero length if strict? 
                   % Original code used <=, but intersection of [0 10] [10 20] is 10 (point). 
                   % Usually intervals are [start, stop). If [start, stop] then point.
                   % ints class doesn't strictly define open/closed, but subtraction logic implies
                   % simple comparison.
                   % Let's keep consistent with valid intervals having duration > 0 or start <= stop?
                   % The constructor enforces start <= stop.
                   outStart(end+1) = start_; %#ok<AGROW>
                   outStop(end+1) = stop_; %#ok<AGROW>
                elseif start_ == stop_
                   % Point intersection?
                   % If we want to keep point intersections, we keep it.
                   outStart(end+1) = start_; %#ok<AGROW>
                   outStop(end+1) = stop_; %#ok<AGROW>
                end
                
                if stopA(iInt) < stopB(iQ)
                    iInt = iInt + 1;
                else
                    iQ = iQ + 1;
                end
            end
            
            objOut = intervals(outStart, outStop);
        end
        
        % -----------------------------------------------------------------
        function objOut = union(obj, query)
            %UNION Union of two interval sets (A | B).
            if ~isa(query, 'intervals'), query = intervals(query); end
            
            combStart = [obj.start; query.start];
            combStop = [obj.stop; query.stop];
            
            temp = intervals(combStart, combStop);
            objOut = temp.consolidate();
        end
        
        % ================================================================
        %  DATA OPERATIONS
        % ================================================================
        
        % -----------------------------------------------------------------
        function [inInts, intIdx] = contains(obj, points)
            %CONTAINS Check if points lie within intervals.
            %
            %   [inInts, intIdx] = contains(obj, points)
            %
            %   OUTPUT:
            %       inInts  - (logical) Nx1 Logic mask. True if point is inside any interval.
            %       intIdx  - (double) Nx2 Matrix. Adjacency List (Coordinate Format).
            %                 Column 1: Point Index
            %                 Column 2: Interval Index
            
            points = double(points);
            nP = length(points);
            inInts = false(nP, 1);
            intIdx = zeros(0, 2);
            
            if obj.count == 0, return; end
            
            % Sweep Path for Overlaps, Discretize for Fast Simple.            
            if obj.overlap
                % ---------------------------------------------------------
                % Path 1: Sweep (Supports Overlaps)
                % ---------------------------------------------------------
                
                % Sort intervals by start, track original IDs.
                n = obj.count;
                [sortedStart, sortIdx] = sort(obj.start);
                sortedStop = obj.stop(sortIdx);
                origIds = sortIdx;
                
                % Sort points
                [sortedP, pSortIdx] = sort(points(:));
             
                % Preallocate Result (Dynamic)
                estSize = nP; % Initial guess
                outP = zeros(estSize, 1);
                outI = zeros(estSize, 1);
                cntOut = 0;
                
                % Mask result
                tempMask = false(nP, 1);
                
                iInt = 1;
                active = []; % List of indices into sorted arrays
                
                for iP = 1:nP
                    p = sortedP(iP);
                    origPIdx = pSortIdx(iP);
                    
                    % Add new intervals that start before or at p
                    while iInt <= n && sortedStart(iInt) <= p
                        active(end+1) = iInt; %#ok<AGROW>
                        iInt = iInt + 1;
                    end
                    
                    % Check active intervals
                    if ~isempty(active)
                        % Remove expired: stop < p
                        validMask = sortedStop(active) >= p;
                        active = active(validMask);
                        
                        % Add matches
                        if ~isempty(active)
                            tempMask(iP) = true;
                            
                            if nargout > 1
                                nNew = length(active);
                                if cntOut + nNew > length(outP)
                                    % Grow - force column expansion
                                    newSize = max(length(outP)*2, cntOut + nNew + 1000);
                                    outP(newSize, 1) = 0;
                                    outI(newSize, 1) = 0;
                                end
                                
                                % Add
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
                    % Sort by Point Index for neatness
                    [~, sortRes] = sort(intIdx(:,1));
                    intIdx = intIdx(sortRes, :);
                end
                
            else
                % ---------------------------------------------------------
                % Path 2: Discretize (Fast, No Overlaps)
                % ---------------------------------------------------------
                
                % Edges: [start; inf]
                % Bin K means: start(K) <= p < start(K+1)
                bins = discretize(points, [obj.start; inf]);
                
                % Map bin directly to interval index (1-to-1)
                candIdx = bins;
                
                % Filter valid candidates (points < s1 are NaN)
                valid = ~isnan(candIdx);
                
                if ~any(valid)
                    return;
                end
                
                cIdx = candIdx(valid);
                pValid = points(valid);
                pIndices = 1:nP;
                pIndices = pIndices(valid)';
                
                % Check Stop
                isIn = pValid <= obj.stop(cIdx);
                
                % Final Selection
                finalP = pIndices(isIn);
                finalI = cIdx(isIn);
                
                % Create Mask
                inInts = false(nP, 1);
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
            %   [dataOut, shifts] = restrict(obj, data, 'flgShift', true)
            %
            %   INPUTS:
            %       data     - (double) NxM Matrix. Col 1 is time.
            %       flgShift - (bool) If true (default), shifts time to eliminate gaps.
            %
            %   OUTPUTS:
            %       dataOut  - (double) Restricted (and optionally shifted) data.
            %       shifts   - (double) Nx4 Matrix [OldStart, OldStop, NewStart, NewStop]
            %                  describing the mapping of intervals. Empty if flgShift=false.
            
             % Parse
            p = inputParser;
            addParameter(p, 'flgShift', true, @(x) islogical(x) || x==0 || x==1);
            parse(p, varargin{:});
            flgShift = p.Results.flgShift;
            
            shifts = [];
            if isempty(data), dataOut = []; return; end
            
            col1 = data(:,1);
            
            % Filter
            % Contains Output 1 is Mask
            mask = obj.contains(col1);
            
            if ~any(mask)
                dataOut = [];
                return;
            end
            
            dataOut = data(mask, :);
            
            if ~flgShift
                return;
            end
            
            % Shift Logic
            % We need to consolidate to determine the new timeline.
            cons = obj.consolidate();
            
            % Assign each point to a consolidated interval.
            % Since cons is non-overlapping, 'contains' (Output 2) returns Adjacency List.
            % Since filtered, every point is in some consolidated interval.
            % [~, adj]
            [~, adjCons] = cons.contains(dataOut(:,1));
            
            % adjCons is [PIdx, IntIdx]
            assert(size(adjCons, 1) == size(dataOut, 1), 'Consolidated contains mismatch');
            
            % Sort by point index to ensure matching with dataOut
            [~, sortOrd] = sort(adjCons(:,1));
            intIds = adjCons(sortOrd, 2);
            
            % Durations and New Starts
            durs = cons.stop - cons.start;
            newStarts = [0; cumsum(durs(1:end-1))];
            
            % Current Starts
            oldStarts = cons.start(intIds);
            targetStarts = newStarts(intIds);
            
            % Shift
            offset = dataOut(:,1) - oldStarts;
            dataOut(:,1) = targetStarts + offset;
            
            % Output Shifts Info
            shifts = [cons.start, cons.stop, newStarts, newStarts + durs];
        end
        
        % ================================================================
        %  OVERLOADS
        % ================================================================
        function obj = and(obj, other); obj = intersect(obj, other); end
        function obj = or(obj, other); obj = union(obj, other); end
    end
end
