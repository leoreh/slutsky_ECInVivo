function [wv] = align_waveforms(wv,options)
% align all waveforms to the waveform with the biggest sink.
% Assume each row is a diffrent waveform, each col is a diffrent sample.

arguments
    wv (:,:) double {mustBeNumeric}
    options.ref_type  (1,:) {mustBeTextScalar,mustBeMember(options.ref_type,["biggest_sink","latest_sink","biggest_amp"])} = "biggest_sink"
    options.use_idx   (:,1) {mustBeScalarOrEmpty,mustBePositive, mustBeInteger} = []
end

if ~isempty(options.use_idx)
    % if user specified requested ref_idx, use it
    if options.use_idx > size(wv,1)
        error("Invalid name-value argument 'use_idx'. Value must be a wv row index.") 
    else
        ref_idx = options.use_idx;
    end
else
    switch options.ref_type
        case "biggest_sink"
            % find the waveform with the biggest sink
            sink_size = min(wv,[],2);
            [~,ref_idx] = max(abs(sink_size));
        case "latest_sink"
            % find the waveform whoes biggest sink is the latest
            [~,sink_loc] = min(wv,[],2);
            [~,ref_idx] = max(sink_loc);
        case "biggest_amp"
            % find the waveform that have the biggest amplitude (max-min diff)
            amp_size = max(wv,[],2) - min(wv,[],2);
            [~,ref_idx] = max(abs(amp_size));
    end     
end

% align all other waveforms to it
for iWave = size(wv,1):-1:1
    % using constant reference, move all to match. based on:
    % [~,wv(iWave,:)] = alignsignals(wv(ref_idx,:),wv(iWave,:),[],'truncate');
    
    % find delay - how many steps to move current to the right, to match ref
    d = finddelay(wv(iWave,:),wv(ref_idx,:));
    
    % move non-refernce signal (using nan)
    if d > 0
        % move to the right
        wv(iWave,:) = [nan(1,d), wv(iWave,1:(end-d))];
    else
        % move to the left
        d = -d;
        wv(iWave,:) = [wv(iWave,(d+1):end), nan(1,d)];
    end  
    
end

end