
function   all_su_bursts=single_unit_bursts(data, isiThr, minSpksBrst) % usually (data, 20 (ms), 2) %I added data_h and data_n
%%%%All burst functions,except for "find_short_bursts" are here, sequently
%data - basic data struct including information regarding recording time, number of units, spike times etc.
%min_ISI - min inter-spike interval in ms for spikes to be considered a burst
%min_spikes_per_burst - min number of spikes for a burst to be accepted
%% preallocations
spktimes = data; % change to seconds
recLen = max(cellfun(@max, spktimes, 'uni', true));


isiThr = 0.02;  % units must be as spktimes
binsize = 3600;  
% bins
bins = n2chunks('n', recLen, 'chunksize', binsize, 'overlap', [1 0]);
bins(1) = 0;
nbins = size(bins, 1);
nunits = length(spktimes); % get from spktimes

% initialize. change recLen according to binsize
sep_unit_spikes = cell(nbins, nunits); % I am thinking theat maxunit fuu hours gives the numcet og hours the uni "exsits" in the recording
all_burst_data = cell(nbins, nunits);
mean_spikes_in_burst = zeros(nbins, nunits);

mean_burst_duration = zeros(nbins, nunits);
mean_in_burst_freq = zeros(nbins, nunits);
mean_percent_spikes_in_bursts = zeros(nbins, nunits);
mean_IBI = zeros(nbins, nunits);
mean_short_bursts = zeros(nbins, nunits);
mean_burst_freq = zeros(nbins, nunits);
mean_norm_burst_freq= zeros(nbins, nunits);

%% make cell array of spikes per unit per time
for ibin=1:nbins
    for u=1:nunits
        [r] = find(spktimes{u}(:)>3600*(ibin-1) & spktimes{u}(:)<3600*ibin); % I think 1200 is becuse the calculation is only over the first 20 min of each hour
        sep_unit_spikes{ibin,u} = spktimes{u}(r);
    end
end
%%
for ibin=1:nbins
    
    for u=1:nunits
        spikes = sep_unit_spikes{ibin,u}(:);
        [burst_data] = single_burst_detect(spikes, isiThr, minSpksBrst);
        all_burst_data{ibin,u} = burst_data;
        mean_spikes_in_burst(ibin,u) = mean(burst_data.spikes_in_burst);
%         spikes_in_bursts.(my_field)= [spikes_in_bursts.(my_field);burst_data.spikes_in_burst];
        mean_burst_duration(ibin,u) = mean(burst_data.burst_duration);
        mean_in_burst_freq(ibin,u) = mean(burst_data.in_burst_freq);
        mean_percent_spikes_in_bursts(ibin,u) = mean(burst_data.percent_spikes_in_bursts);
        mean_IBI(ibin,u) = mean(burst_data.IBI)/1000; %%%% divided by 1000 so that the answer will be in seconds
        mean_short_bursts(ibin,u) = mean(burst_data.short_bursts);
        mean_burst_freq(ibin,u) = mean(burst_data.burst_freq);
        mean_norm_burst_freq(ibin,u)= mean(burst_data.norm_burst_freq);
    end
end
%% remove non-bursting units
for u=1:nunits
    if all(mean_burst_freq(:,u)<0.004)
        mean_spikes_in_burst(:,u) = NaN;
        mean_burst_duration(:,u) = NaN;
        mean_in_burst_freq(:,u) = NaN;
        mean_percent_spikes_in_bursts(:,u) = 0; %%chaned NaN to 0
        mean_IBI(:,u) = NaN;
        mean_short_bursts(:,u) = mean(burst_data.short_bursts);
        mean_burst_freq(:,u) = NaN;
        mean_norm_burst_freq(:,u)= NaN;
        
    end
end
all_su_bursts.mean_spikes_in_burst = mean_spikes_in_burst;

all_su_bursts.mean_burst_duration = mean_burst_duration;
all_su_bursts.mean_in_burst_freq = mean_in_burst_freq;
all_su_bursts.mean_percent_spikes_in_bursts = mean_percent_spikes_in_bursts;
all_su_bursts.mean_percent_spikes_in_bursts( isnan(all_su_bursts.mean_percent_spikes_in_bursts)) =0; %%Added on 21/1/19 by Max  
all_su_bursts.mean_IBI = mean_IBI;
all_su_bursts.mean_burst_freq = mean_burst_freq;
all_su_bursts.mean_norm_burst_freq = mean_norm_burst_freq;
end

