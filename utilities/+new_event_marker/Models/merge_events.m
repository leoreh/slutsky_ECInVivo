function event_list = merge_events(event_list, overlap_matrix)
    % MERGE_EVENTS Merge overlapping events in a list.
    %
    %   event_list = merge_events(event_list, overlap_matrix)
    %   Merges overlapping events in the provided event list.
    %   If overlap_matrix is not provided, it will be calculated internally.
    %
    %   Inputs:
    %       event_list: An Nx2 numeric matrix where each row represents an event with the
    %                   first column containing the start time and the second column
    %                   containing the end time.
    %       overlap_matrix: (Optional) An NxN logical matrix where overlap_matrix(i, j)
    %                       is true if event i overlaps with event j, and false otherwise.
    %                       Diagonal elements (i.e., self-overlaps) must be false.
    %                       Default is calculated internally using check_overlap function.
    %
    %   Output:
    %       event_list: A modified event list with overlapping events
    %                   merged. The overlapping events will be removed, and
    %                   the merged events will be appended to the list end.
    %
    %   Example:
    %       event_list = [1, 5; 3, 7; 6, 9];
    %       event_list = merge_events(event_list);
    %
    % Author: Lior de Marcas (LdM: https://github.com/Liordemarcas)
    % Date: 2024-02-18

    if nargin == 1
        overlap_matrix = check_overlap(event_list);
    end

    nEvents = size(event_list,1);
    events_with_overlap = any(overlap_matrix,1);
    iIter = 0;
    
    while any(overlap_matrix(:))
        
        merged_events = double.empty(0,2);
        events_used = [];
        for iEv = 1:sum(events_with_overlap)
            % collect 1 event & what events it overlap with
            event2work = find(events_with_overlap,1);
            overlapping_events = overlap_matrix(:,event2work);
            overlapping_events(event2work) = true; % include current event
            
            % skip if we have an event we already used - 
            % cases of sequantail overlaps ([1 3; 2 5; 4 6]) will be taken
            % care of next iteration of the while loop
            if ismember(event2work,events_used)
                continue
            end

            % find what the merged event is
            overlapping_events_t = event_list(overlapping_events,:);
            merge_start = min(overlapping_events_t(:,1));
            merge_end = max(overlapping_events_t(:,2));

            % save events that were used, to skip them later
            events_used = [events_used;find(overlapping_events)]; %#ok<AGROW> this is unavoidable
            
            % save merged event
            merged_events(end+1,:) = [merge_start, merge_end]; %#ok<AGROW> this is unavoidable

            % move to next event using the "find" in the loop start
            events_with_overlap(event2work) = false;
        end

        % update event_list
        event_list(events_used,:) = [];
        event_list = [event_list;merged_events]; %#ok<AGROW> this is unavoidable

        % get new overlap_matrix
        overlap_matrix = check_overlap(event_list);
        events_with_overlap = any(overlap_matrix, 1);
        
        % make sure not to run forever - if this error is raised, something
        % probably went wrong
        iIter = iIter + 1;
        if iIter > 2*nEvents
            error("Didn't converge after %d iterations",iIter)
        end
    end

end
