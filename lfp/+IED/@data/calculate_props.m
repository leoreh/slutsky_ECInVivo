function [events_width_half_amp, events_amp, event_width_10_amp] = calculate_props(obj, pos2calc, margs)
    % calculate_props - Calculate properties of interictal epileptiform discharges (IEDs)
    %
    % This function calculates various properties of interictal epileptiform
    % discharges (IEDs) detected in EEG signals. The properties include the
    % amplitude, width at half amplitude, and width at 10% of the amplitude.
    %
    % SYNTAX:
    % [events_width_half_amp, events_amp, event_width_10_amp] = calculate_ied_props(obj, margs)
    %
    % INPUT ARGUMENTS:
    %   obj     - An object of class IED.data containing EEG signal and related data.
    %  pos2calc - A logical vector, which detections in the ied object to
    %             calculate propeties on.
    %             Defaults to obj.accepted if not provided.
    %   margs   - (Optional) Double array representing additional parameters.
    %             Defaults to obj.marg if not provided.
    %
    % OUTPUT ARGUMENTS:
    %   events_width_half_amp   - Width at half amplitude for each detected IED.
    %   events_amp              - Amplitude of each detected IED.
    %   event_width_10_amp      - Width at 10% of the amplitude for each detected IED.
    %
    % NOTES:
    %
    % EXAMPLE USAGE:
    % % Assuming obj is an instance of IED.data
    % [width_half_amp, amp, width_10_amp] = calculate_ied_props(obj);
    %
    % SEE ALSO:
    % - IED.data
    %
    % Author: LdM
    

    arguments
        obj IED.data
        pos2calc logical = obj.accepted
        margs double = obj.marg
    end
    
    if ~isvector(pos2calc) || (numel(pos2calc) ~= numel(obj.pos))
        error("pos2calc must be a vector with same number of elements as in 'pos' in the IED.data object")
    end
    win_means = movmean(obj.sig(:),obj.fs);
    true_pos = obj.pos(pos2calc);
    local_segs = obj.extract_discharges(margs, pos2calc);
    detec_p = ceil(size(local_segs,2)/2);

    for iEvt = numel(true_pos):-1:1

        % calculate amplitude - marked peak (window center) to window mean
        event_peak = local_segs(iEvt,detec_p);
        event_mean = win_means(true_pos(iEvt));
        events_amp(iEvt) = abs(event_peak - event_mean);

        % calculate width
        half_amp = events_amp(iEvt)/2;
        tnth_amp = events_amp(iEvt)/10;
        if event_peak > event_mean
            % detection is above signal mean
            events_width_half_amp(iEvt) = width_from_thrshold(event_mean + half_amp, @gt);
            event_width_10_amp(iEvt) = width_from_thrshold(event_mean + tnth_amp, @gt);
        else
            % detection is below signal mean
            events_width_half_amp(iEvt) = width_from_thrshold(event_mean - half_amp, @lt);
            event_width_10_amp(iEvt) = width_from_thrshold(event_mean - tnth_amp, @lt);
        end

    end


    function [spike_width] = width_from_thrshold(cross_val, cross_fun)
        % Nested function! calculate spike width at a given threshold.
        % Inputs:
        %   cross_val - what value width is calculated at.
        %   cross_fun - either @lt ( < ) or @gt ( > ), what is the
        %               directions of the spike compared to cross value.
        % Outputs:
        %   spike_width - time diff between points of calculations

        % find when crossing value is passed around spike center (detect_p)
        right_thr2_pass = cross_fun(local_segs(iEvt,detec_p:end), cross_val);
        left_thr2_pass = cross_fun(local_segs(iEvt,1:detec_p), cross_val);

        % find the last sequantial threshold cross around the spike center
        left_thr_mark = find(diff(left_thr2_pass),1,"last") + 1;
        right_thr_mark = find(diff(right_thr2_pass),1,"first") + detec_p - 1;

        % calculate spike width in [s]
        spike_width = (right_thr_mark-left_thr_mark)./obj.fs;
        if isempty(spike_width)
            spike_width = nan;
        end
    end

end


% graphics - display
%{
clearvars ans
[~,ans{1}] = ied.extract_discharges(margs);
ans{1} = ans{1}*1000;
figure('WindowState','maximized');
plot(ans{1},local_segs(iEvt,:));
%ylim([-10 10])
hold on
scatter(ans{1}(detec_p),local_segs(iEvt,detec_p),'kO','filled')
ans{2} = [left_thr2_pass(1:(end-1)),true,right_thr2_pass(2:end)];
plot(ans{1}(ans{2}),local_segs(iEvt,ans{2}),'r','LineWidth',2)
scatter(ans{1}([left_thr_mark,right_thr_mark]),local_segs(iEvt,[left_thr_mark,right_thr_mark]),'Or','filled','MarkerEdgeColor','k')
yline(cross_val,'--g');
%}

% EOF