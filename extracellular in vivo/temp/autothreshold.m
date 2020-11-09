function [spikes, thresholds, STDthresholds] = autothreshold(signal,smplfreq,varargin)
%
% Function to detect neuronal action potentials in extracellular neuronal
% field potential measurement data. The function analyses the histogram of
% the input data, and sets detection amplitude tresholds automatically.
% Negative and positive thresholds are set independently. 
% [SPIKES,THRESHOLDS,STDTHRESHOLDS] = AUTOTHRESHOLD(signalvector,samplingfrequency)
% 
% AUTHOR
%
% Dr. Jarno M. A. Tanskanen
% Computational Biophysics and Imaging Group
% Department of Electronics and Communications Engineering, and BioMediTech
% Tampere University of Technology
% Finland
% E-mail: tanskanen //at// ieee.org
%
% VERSION
%
% Version 0.9, July 3, 2016
%
% The function is complete but not exhaustively tested. Comments may also
% be incomplete.
%
% Also checking the histogram gradient's properties is insufficient for 
% general use, i.e., histograms with no usable features are not screened 
% for; this has been planned for a future version. Currently, the user must 
% observe the gradient images to check that the gradient possesed features 
% resulting in useful thresholds.
%
% (Capitalization for parameters and the function name are for highlighting
% only. Please, write all using lower case when used.)
%
% USAGE AND PARAMETERS
%
% AUTOTHRESHOLD is a function to detect action potentials spikes by
% objectively thresholding the input signal. The input signal is assumed to
% be a field potential signal exhibiting neuronal (or other) action
% potential signal, with the desired action potentiathrs reaching
% sufficiently much above the background noise level.
%
% Usage:
%
% [SPIKES,THRESHOLDS,STDTHRESHOLDS] = AUTOTHRESHOLD(S,F) threshold the data
% in vector S using the default parameter values. F is the sampling
% frequency of S in samples per second. The output vector SPIKES carries
% the spike time stamps in samples, and THRESHOLDS is a vector of length two
% with the two thresholds. STDTHRESHOLDS is a 2-element vector with
% the negative and postive thresholds in input signal standard deviations.
% S and F must be supplied.
%
% [...] = AUTOTHRESHOLD(S,F,STRT,STP) finds the threholds
% based on the analyses of only the section S(STRT:STP) and thereafter
% thresholds the signal in vector S accordingly.
%
% [...] = AUTOTHRESHOLD(S,F,STRT,STP,VARARGIN) analyzes S
% using NTHRESHOLDS threshold values between the minimun and maximun of S.
%
% Optional Parameters (to be written to VARARGIN)
%
% Parameter name        Values and description
%
% 'nthresholds'          The total number of threshold values to use in the
%                       analysis. Default is 100.
%
% 'start'               The sample number of the first sample to analyze. 
%                       Default is 1.
%
% 'stop'                The sample number of the last sample to analyze.
%                       Default is lenth(s).
%
% 'graphs'              1 - plot graphs (default)
%                       0 - do not plot graphs
%
% 'refract'             Refractory period; time in ms from the time the
%                       signal exceeds the positive threshold or falls
%                       below the negative threshold until the next spike
%                       may be detected. Default is 1.5.
%
% 'stdlimits'           2-element vector with the limits for the
%                       automatically found thresholds in standard
%                       deviations of the input signal. Default is [3 10].
%
% For example:
% [spikes,thresholds,stdthresholds] = autothreshold(s,f,'Naverager',10,'strt',10001,'stp',20000)
%
% CITING
%
% When describing any use or modification of the program or the method,
% whether in its original or derivative form, please, cite [1]. The method
% originally proposed in [2] and a version of the function was also
% presented at [3], which may optionally be sited in addition to [1].
%
% COPYRIGHT
%
% Copyright (c) 2016, Jarno M. A. Tanskanen
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution
%     * Neither the name of the Uni nor the names
%       of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written
%       permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% (Copyright text is a copy of the text appearing in the Matlab Central
% File Exchange.)
%
% REFERENCES
%
% [1] Jarno M. A. Tanskanen, Fikret E. Kapucu, Inkeri Vornanen, and Jari A.
% K. Hyttinen, "Automatic objective thresholding to Detect Neuronal Action
% Potentials," accepted for publication in Proc. 24th European Signal 
% Processing Conference (EUSIPCO 2016), Aug. 29 - Sept. 2, 2016, Budapest, 
% Hungary.
%
% [2] J. M. A. Tanskanen, F. E. Kapucu, and J. A. K. Hyttinen, "On the
% threshold based neuronal spike detection, and an objective criterion for
% setting the threshold," in Proc. 7th Annual International IEEE EMBS
% Conference on Neural Engineering, Montpellier, France, Apr. 2015, pp.
% 1016?1019. http://dx.doi.org/10.1109/NER.2015.7146799
%
% [3] J. M. A. Tanskanen, F. E. Kapucu, I. Vornanen, K. Lenk, J. Hyttinen,
% "Objective thresholding of MEA data for action potential detection, "
% Front. Neurosci. Conference Abstract: MEA Meeting 2016 | 10th
% International Meeting on Substrate-Integrated Electrode Arrays. 
% http://dx.doi.org/10.3389/conf.fnins.2016.93.00028
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform some basic checing of the input data validity.
% This is not nearly exhaustive at the moment.
if nargin < 2,
    error('Input error: Too few input arguments.');
end
if length(size(signal)) > 2,
    error('Input data error: Input data must be a vector.');
end
if min(size(signal)) > 1,
    error('Input data error: Input data must be a vector.');
end
if any(any(isnan(signal))),
    error('Input data error: Input data contains NaN''s.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Default values for optional parameters
graphs = 1;           % Graphs on
strt = 1;             % Start analysing from the beginning of the input
% signal
Ls = length(signal)
stp = smplfreq*60     % By default, analyze 1 min
if stp > Ls           % If signal shorter than 1 min, analyze until the end
    stp = Ls          % of the input signal
end;
Nthresholds = 500;    % The number of threshold values for making the
% histogram
Naverager = 10;       % The '10' will be made apdaptive in the
% following versions, since the features
% of the smoothed gradient depends on
% the ratio of the smoothing averager
% length (Naverager) and the number of
% thresholds (Nthresholds)
refract = 1.5;        % Refractory period in milliseconds
STDlimits = [3 10];   % Lower and upper limits for acceptable threholds in
% standard deviatons of the entire signal to analyze
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get possible user input values of optional parameters
if (rem(length(varargin),2)==1)
    error('Input error: optional parameters must occur in pairs.');
else
    
    for i=1:2:(length(varargin)-1)
        
        if ~ischar(varargin{i}),
            error('Input error: Optional parameter names must be strings).');
        end
        
        switch lower(varargin{i})
            
            case 'nthresholds'
                Nthresholds = varargin{i+1};
                
            case 'graphs'
                graphs = varargin{i+1};
                
            case 'stp'
                if strt < 1
                    error('Input error: Analysis end sample number is too large.');
                end;
                stp = varargin{i+1};
                
            case 'strt'
                if strt > length(signal)
                    error('Input error: Analysis start sample number is too large.');
                end;
                strt = varargin{i+1};
                
            case 'refract'
                refract = varargin{i+1};
        end;
    end;
end;
d = signal(strt:stp); % Data vector to analyze
STDsignal = std(d);   % Standard deviation of the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find threshold values equispaced between the min and max of the input
% signal.
thrs = linspace(min(d),max(d),Nthresholds);
% Initialize histogram data vector (spkcnts for all one-sided
% thresholds)
spkcnt = zeros(size(thrs));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Form histograms
% Find spike counts for all negative thresholds
m = 1;
for i = thrs(thrs < 0);
    sn = signal < i; % Logical '1' for all samples smaller than the
    % threshold
    spkcnt(m) = length(find(diff(sn)==1)); % Each continuous string of 1s
    % is considered to be one spike.
    % Count only rising edges.
    m = m + 1;
end;
% Find spike counts for zero and all positive thresholds, write them to the
% same vector as the negative thresholds, after the negative thresholds
for i = thrs(thrs >= 0);
    sp = signal > i; % Logical '1' for all samples greater than the
    % threshold
    spkcnt(m) = length(find(diff(sp)==1)); % Each continuous string of 1s
    % is considered to be one spike.
    % Count only rising edges.
    m = m + 1;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the gradient of the spike count histogram and smooth with
% a moving averager, filtering forthward and backward to cancel phase
% errors due to filtering
gr = gradient(spkcnt);
hf = ones(1,Naverager)./Naverager;
grf = filtfilt(hf,1,gr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analyze the gradient of the histogram (grf)
[~, locsp] = findpeaks(grf); % find positive peaks
[~, locsn] = findpeaks(-grf); % find negative peaks
% Locations of the maximun and minimun of the smoothed gradient
mxpos = find(grf==max(grf));
mxneg = find(grf==min(grf));
% Negative threshold is at the negative peak closest to the maximum of the
% smoothed gradient
negth = max(locsn(locsn<mxpos));
% Positive threshold is at the positive peak closest to the maximum of the
% smoothed gradient
posth = min(locsp(locsp>mxneg));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot histogram and its gradient, if desired.
if (graphs == 1)
    % Plot, in a new figure, the histogram with the left vertical axis, and
    % smoothed gradient with the right vertical axis, and plot the
    % thresholds on them
    figure;
    [h, l1, l2]=plotyy(thrs,spkcnt,thrs,grf);
    hold on;
    
    % Set the axes nicely
    minthrs = min(thrs); % smallest threshold
    maxthrs = max(thrs); % largest threshold
    y1 = ylim(h(1)) % current vertical axis limit for the left axis
    y2 = ylim(h(2)); % current vertical axis limit for the right axis
    axis(h(1),[minthrs maxthrs -y1(2) y1(2)]);
    axis(h(2),[minthrs maxthrs y2]);
    
    % Set the linewidth and colors nicely
    set(l1,'LineWidth',2,'Color','k')
    set(l2,'LineWidth',2,'Color','r')
    set(h(1),'YColor','k')
    set(h(2),'YColor','r')
    % Set the axis labels
    ylabel(h(1),'Spike count','Color','k'); % label left y-axis
    ylabel(h(2),'Gradient','Color','r'); % label right y-axis
    xlabel(h(2),'Amplitude') % label x-axis
    
    % Plot vertica lines at the thresholds on the histogram
    hl1 = line([thrs(negth) thrs(negth)],[y1(2) -y1(2)]);
    hl2 = line([thrs(posth) thrs(posth)],[y1(2) -y1(2)]);
    
    % Set the line properties nicely
    set(hl1,'LineWidth',2,'Color','m')
    set(hl2,'LineWidth',2,'Color','g')
    
    hold off;
    
    
    
    % Plot, in a new figure, the data and the thresholds
    figure;
    t = (1:length(d))/smplfreq*1000;
    plot(t,d,'k');
    
    hold on;
    plot(t,ones(size(t))*thrs(negth),'m','LineWidth',2);
    plot(t,ones(size(t))*thrs(posth),'g','LineWidth',2);
    ylabel('Amplitude');
    xlabel('Time (ms)')
    x = xlim;
    y = ylim;
    axis([x(1) max(t) y]) % The limits will be automatized
    
    hold off;
    
end;
%%
% Find times in samples when s is greater than or equal to the positive
% threshold, or smaller than or equal to the negative threshold
spsz = zeros(size(signal));
spsz((signal >= posth) | (signal <= negth)) = 1; % binary string with 1s
% where a threhold is
% satisfied
% spsz = spsz;
% Find all peaks (local minima and maxima) in data
[~, ppeaklocs] = findpeaks(signal);
[~, npeaklocs] = findpeaks(-signal);
peaklocs = sort([npeaklocs ppeaklocs]); % all peak locations in samples
peakz = zeros(size(signal));
peakz(peaklocs) = 1; % binary string with 1s at local maxima
%%
% Spike time points are where there is a peak in the signal and the signal
% satisfies a threshold
spks = spsz & peakz;
%%
% Apply refractory period (this part could be made prettier)
rfper = floor(smplfreq/1000*refract); % number of samples in
% a refractory period rounded down to
% a full sample (next sample will
% already be outside the refractory
% period)
l = length(spks); % length of data for cutting out extras possibly
% produced by the following procedure
spklocs = find(spks == 1); % spike time stamps in samples
% Go through all the spike time stamps one by one and set all samples to
% zero within the refractory period, if the spike was not removed by an
% earlier refractory period application
for m = spklocs
    if (spks(m) == 1) % check that the spike was not
        % removed in a previous
        % application of the refractory period
        spks(m + 1 : m + rfper) = 0; % set all samples to
        % zero starting from with the sample
        % right after a spike and ending at
        % the end of the refractory period
    end;
end;
spks = spks(1:l); % remove extra samples possibly created by the
% application of the refractory period
%%
% Form outputs and check that the found thresholds are within
% reasonable limits
thresholds = [thrs(negth) thrs(posth)];
limits = STDlimits*STDsignal;
if (thresholds(1) > -limits(1) ||  thresholds(1) < -limits(2) || ...
        thresholds(2) < limits(1) || thresholds(2) > limits(2))
    disp('************************************************');
    disp(' ');
    disp('WARNING! Found threshold beyond validity limits.');
    disp(' ');
    disp('Verify the validity of the results and/or apply');
    disp('user-defined thresholds validity limits.');
    disp(' ');
    disp('************************************************');
    
end;
spikes = spks;
STDthresholds = thresholds/STDsignal;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%