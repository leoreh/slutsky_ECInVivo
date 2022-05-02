
%%
function [burst_data]=single_burst_detect(spikes, min_ISI, min_spikes_per_burst)
ISI=diff(spikes);
[r]= find(ISI<=min_ISI);
bs=[];
for i=1:numel(r)
   
    bs(i)=r(i);
    if ~isequal(i,1) % start from the seconed ISI
        while isequal(r(i-1), (r(i)-1)) %%in bs-nullyfies the indices of the ISIs that belong to the same burst , except for the first ISI of each burst.
            bs(i)=0;
            i=i+1;
            if i>=numel(r)
                break
            end
        end
    end
end
if ~isempty(bs);
    bs_r=find(bs);
    bs=bs(bs_r);%bs is now the indices of the first ISI of each burst
    be=zeros(numel(bs),1);%
    for i=1:numel(bs)
        if isempty(find(ISI(bs(i):numel(ISI))>min_ISI)) %if there are no too-long ISIs since the i's burst begining, so be(i) is all the ISIs (one long burst). Important if the ISI ends with a burst, for example.
            be(i)=numel(ISI);
        else
            be(i)=(find(ISI(bs(i):numel(ISI))>min_ISI, 1, 'first'))+(bs(i))-1;
        end
    end
end
%%
%%single_burst_stats;
if ~isempty(bs)
    burst_spike_num = [bs' be];
    spikes_in_burst = diff(burst_spike_num,1,2)+1;%number of ISI in burst
    [short_bursts, filt_burst_spike_num, filt_spikes_in_burst]= find_short_bursts(min_spikes_per_burst, burst_spike_num,  spikes_in_burst);
    burst_times(:,1) = spikes(filt_burst_spike_num(:,1));
    burst_times(:,2) = spikes(filt_burst_spike_num(:,2));
    burst_duration = diff(burst_times,1,2)*1000; 
    in_burst_freq = filt_spikes_in_burst./(burst_duration/1000); 
    percent_spikes_in_bursts = (sum(filt_spikes_in_burst)/numel(spikes));
    
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
    burst_data.norm_burst_freq= burst_freq/(numel(spikes)/3600);
else
    burst_data.burst_spike_num = [];
    burst_data.burst_times = [];
    burst_data.spikes_in_burst = [];
    burst_data.burst_duration = [];
    burst_data.in_burst_freq = [];
    burst_data.percent_spikes_in_bursts = [];
    burst_data.IBI = [];
    burst_data.short_bursts = [];
    burst_data.burst_freq = 0;
    burst_data.norm_burst_freq= 0;
end
%%
%%
end