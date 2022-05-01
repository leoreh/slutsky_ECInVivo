function [short_bursts, filt_burst_spike_num, filt_spikes_in_burst]= find_short_bursts(min_spikes_per_burst, burst_spike_num,  spikes_in_burst); 

[r] = find(spikes_in_burst >= min_spikes_per_burst);
short_bursts = (numel(spikes_in_burst)-numel(r))/numel(spikes_in_burst);
filt_burst_spike_num = burst_spike_num(r,:);
filt_spikes_in_burst = spikes_in_burst(r,:);

end