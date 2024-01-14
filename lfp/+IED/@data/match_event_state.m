function [pos_labels] = match_event_state(obj, labels, accepted)
    % MATCH_EVENT_STATE - Find the brain state that matches each Interictal
    % Epileptiform Discharge (IED) event based on provided labels.
    %
    %   [pos_labels] = MATCH_EVENT_STATE(obj, labels, accepted)
    %
    % INPUTS:
    %   obj         - An object of the IED.data class.
    %   labels      - An array of integers (1 to 8) representing the brain
    %                 states corresponding to different time intervals.
    %   accepted    - (Optional) A logical array indicating whether each IED
    %                 event is accepted. Default is 'obj.accepted'.
    %
    % OUTPUT:
    %   pos_labels  - An array containing the matched brain state for each
    %                 accepted IED event. Unmatched events are assigned the
    %                 value 8.
    %
    % NOTES:
    %   - The function ensures that the length of the signal and the number of
    %     provided labels differ by no more than 1 second.
    %   - The function assumes that labels correspond to a round number of
    %     seconds, where label 1 is for the interval [0, 1], label 2 for [1, 2],
    %     and so on.
    %   - Events occurring exactly on a second will be associated
    %     with the next second. For example, an event at 4.00
    %     seconds will be included in the label for second 4 to 5
    %     (not 3 to 4).
    %
    % EXAMPLES:
    %   % Using default 'accepted' values
    %   labels = [1 2 3 4 5 6 7 8];
    %   pos_labels = match_event_state(obj, labels);
    %
    %   % Providing custom 'accepted' values
    %   custom_accepted = logical([1 0 1 1 0 1 1 0]);
    %   pos_labels_custom = match_event_state(obj, labels, custom_accepted);
    %
    % See also:
    %   IED.data

    arguments
        obj IED.data
        labels (1,:) double {mustBeMember(labels,1:8)}
        accepted (1,:) logical {mustBeNumericOrLogical} = obj.accepted
    end
    % validate labels & signal are no more than 1 sec diffrent
    sig_len = length(obj.sig)./obj.fs;
    if ~ismembertol(sig_len, numel(labels), 1, "DataScale", 1)
        error('signal is diffrent than labels by more than 1 second')
    end

    % collect all accpted events
    true_pos = obj.pos(accepted);

    % find each pos matching secound. Note that this assume labels
    % are from round number of secounds - label 1 is [0 1], label 2
    % is [1 2], etc.
    true_pos_s = true_pos./obj.fs;
    matching_sec = floor(true_pos_s) + 1; % +1 due to matlab 1 indexing

    % find positions that exceed the labels, assign 8 to them, and
    % collect matching label to the rest
    pos_labels = 8.*ones(size(true_pos));
    pos_in_labels = matching_sec < numel(labels);
    pos_labels(pos_in_labels) = labels(matching_sec(pos_in_labels));

end