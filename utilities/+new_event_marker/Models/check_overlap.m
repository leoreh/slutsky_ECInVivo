function overlap = check_overlap(events)
    % CHECK_OVERLAP - Determine if there are any overlaps between events
    %
    %   overlap = check_overlap(events)
    %
    %   Inputs:
    %       events: An Nx2 matrix where each row represents an event with the
    %               first column containing the start time and the second column
    %               containing the end time.
    %
    %   Output:
    %       overlap: An NxN logical matrix where overlap(i, j) is true if event i
    %                overlaps with event j, and false otherwise. Diagonal elements
    %                (i.e., self-overlaps) are always false.
    %
    %   Example:
    %       events = [1, 5; 3, 7; 6, 9];
    %       overlap = check_overlap(events);
    %
    % Author: Lior de Marcas (LdM: https://github.com/Liordemarcas)
    % Date: 2024-02-18

    start_times = events(:, 1);
    end_times = events(:, 2);

    % Create matrices of start and end times for efficient broadcasting
    start_matrix = repmat(start_times, 1, size(events, 1));
    end_matrix = repmat(end_times, 1, size(events, 1));

    % Check for overlaps using logical indexing
    contain_overlap_start   = (start_times' <= start_matrix) & (start_matrix <= end_times');
    contain_overlap_end     = (start_times' <= end_matrix)   & (end_matrix <= end_times');
    overlap = contain_overlap_start | contain_overlap_end;
    event_fully_contain = (start_matrix <= start_times') & (end_times' <= end_matrix);
    overlap = overlap | event_fully_contain;
    overlap = overlap & ~eye(size(overlap));  % Remove self-overlaps

    % overlap = sparse(overlap);
end