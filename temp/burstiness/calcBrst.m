
function [burst_data]=single_burst_detect(spktimes, isiThr, minSpksBrst)

isi=diff(spktimes);

% get indices of first and last spike in each burst
brstIdx = binary2epochs('vec', isi <= isiThr);

% single_burst_stats;
nspksBrst = diff(brstIdx, 1, 2) + 1;

[r] = find(nspksBrst >= minSpksBrst);
short_bursts = (numel(nspksBrst)-numel(r))/numel(nspksBrst);
filt_burst_spike_num = brstIdx(r,:);
filt_spikes_in_burst = nspksBrst(r,:);





burst_times(:,1) = spktimes(filt_burst_spike_num(:,1));
burst_times(:,2) = spktimes(filt_burst_spike_num(:,2));
burst_duration = diff(burst_times,1,2)*1000;    % [ms]
in_burst_freq = filt_spikes_in_burst./(burst_duration/1000);
percent_spikes_in_bursts = (sum(filt_spikes_in_burst)/numel(spktimes));

burst_freq = numel(burst_duration)/3600; %me: again 1200 because of the 20 min
if isempty(burst_times) || size(burst_times,1)==1
    IBI=[];
else
    for b=1:size(burst_times,1)-1
        IBI(b,1) = 1000*(burst_times(b+1,1)- burst_times(b,2));
    end

end
burst_data.burst_spike_num = filt_burst_spike_num;
burst_data.burst_times = burst_times;
burst_data.spikes_in_burst = filt_spikes_in_burst;
burst_data.burst_duration = burst_duration;
burst_data.in_burst_freq = in_burst_freq;
burst_data.percent_spikes_in_bursts = percent_spikes_in_bursts;
burst_data.IBI = IBI;

burst_data.short_bursts=short_bursts;
burst_data.burst_freq = burst_freq;
burst_data.norm_burst_freq= burst_freq/(numel(spktimes)/3600);


end