function brst = brst_maxInt(spktimes, varargin)
% spktimes_brstMaxInt - Detects bursts using the Max Interval method
%
% This function implements the Max Interval burst detection algorithm,
% which is defined by fixed thresholds for inter-spike intervals (ISIs) and
% burst properties. It follows the implementation described in Cotterill et
% al. (2016) and the NeuroExplorer manual.
%
% INPUT
%   spktimes    cell array of spike times per unit (e.g., {unit1_times, unit2_times}).
%               Times should be in seconds.
%   varargin    Name-value pairs for parameters:
%       'basepath'      char. fullpath to recording folder {pwd}
%       'maxISI_start'  numeric. Max ISI allowed at start of burst {0.1 [s]}
%       'maxISI_end'    numeric. Max ISI allowed within a burst {0.2 [s]}
%       'minIBI'        numeric. Minimum Inter-Burst Interval {0.2 [s]}
%       'minDur'        numeric. Minimum burst duration {0.05 [s]}
%       'minSpks'       numeric. Minimum number of spikes in a burst {5}
%       'flgSave'       logical. Save struct to file {true}
%       'flgForce'      logical. Analyze even if file exists {false}
%
% OUTPUT
%   brst        struct containing burst statistics and locations.
%               Fields include:
%               .detect     - Number of bursts detected
%               .nspks      - Mean number of spikes per burst
%               .brstDur    - Mean burst duration
%               .freq       - Mean intra-burst frequency
%               .bspks      - Fraction of spikes in bursts
%               .ibi        - Mean inter-burst interval
%               .all        - Struct array with details of every burst
%                             (times, nspks, etc.) for each unit.
%
% ALGORITHM
%   1. Core Burst Detection:
%      - Identify "core" bursts segments where all ISIs are < maxISI_start.
%   2. Extension:
%      - Extend these segments by including subsequent spikes if their ISI
%        is <= maxISI_end.
%   3. Merging:
%      - Merge any two bursts if the interval between them (end of first to
%        start of second) is < minIBI.
%   4. Filtering:
%      - Discard bursts with duration < minDur. - Discard bursts with
%      number of spikes < minSpks.
%
% REFERENCES
%      - Cotterill, E., et al. (2016). A comparison of computational
%        methods for detecting bursts in neuronal spike trains and their
%        application to human stem cell-derived neuronal networks. J.
%        Neurophysiol.
%      - NeuroExplorer Manual, section 2.24.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'basepath', pwd, @ischar)
addParameter(p, 'maxISI_start', 0.1, @isnumeric)
addParameter(p, 'maxISI_end', 0.2, @isnumeric)
addParameter(p, 'minIBI', 0.2, @isnumeric)
addParameter(p, 'minDur', 0.05, @isnumeric)
addParameter(p, 'minSpks', 5, @isnumeric)
addParameter(p, 'flgSave', true, @islogical)
addParameter(p, 'flgForce', false, @islogical)

parse(p, varargin{:})
basepath        = p.Results.basepath;
maxISI_start    = p.Results.maxISI_start;
maxISI_end      = p.Results.maxISI_end;
minIBI          = p.Results.minIBI;
minDur          = p.Results.minDur;
minSpks         = p.Results.minSpks;
flgSave         = p.Results.flgSave;
flgForce        = p.Results.flgForce;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load if exists
[~, basename] = fileparts(basepath);
brstfile = fullfile(basepath, [basename, '.mi_brst.mat']);
if exist(brstfile, 'file') && ~flgForce
    load(brstfile, 'brst')
    return
end

nunits = length(spktimes);

% initialize output structure
brst.info.runtime       = datetime("now");
brst.info.algorithm     = 'MaxInterval';
brst.info.input         = p.Results;

brst.detect             = zeros(1, nunits);
brst.nspks              = nan(1, nunits);
brst.brstDur            = nan(1, nunits);
brst.freq               = nan(1, nunits);
brst.bspks              = zeros(1, nunits);
brst.ibi                = nan(1, nunits);

% Initialize 'all' field as a cell array or struct array to hold individual burst data
% We will use a cell array for simplicity if units have different numbers of bursts,
% or a struct array if acceptable. Let's use a cell array of structs for flexibility.
brst.all = cell(1, nunits);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc burstiness per unit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iunit = 1 : nunits

    st = spktimes{iunit};
    if isempty(st) || length(st) < minSpks
        continue;
    end

    % Ensure column vector
    st = st(:);

    % --- Phase 1: Identify potential bursts based on maxISI_start ---
    isi = diff(st);

    % We need to find sequences of spikes.
    % A burst must start with an ISI < maxISI_start (Strict!).
    % It continues as long as ISIs are <= maxISI_end.

    % However, NeuroExplorer manual says:
    % "The algorithm finds short ISIs. A short ISI is an interval <= Max Interval.
    %  The algorithm finds the start and end of the burst using Max Interval.
    %  ... A burst consists of at least Min Number of Spikes ...
    %  If the interval between two bursts is less than Min Inter-Burst Interval, the bursts are merged."
    %
    % Wait, the Max Interval method usually has TWO interval thresholds:
    % 1. Max Begin ISI (start of burst)
    % 2. Max End ISI (end of burst / continuation)
    %
    % Cotterill et al. 2016 says:
    % "The MI method ... defines a burst as a series of spikes with ISIs less than a threshold...
    %  NeuroExplorer implementation uses ... Max ISI to start a burst ... Max ISI to end a burst"

    bursts = []; % [start_time, end_time, n_spikes]

    nspikes = length(st);
    inBurst = false;
    currentBurstStart = -1;
    currentBurstSpikes = 0;

    % Indices of spikes in the current burst
    curr_spk_idxs = [];

    % Pre-calculate ISIs
    all_isis = diff(st);

    % We will identify bursts by indices of spikes [start_idx, end_idx]
    raw_bursts_idx = [];

    k = 1;
    while k < nspikes
        if ~inBurst
            % Look for start of a burst
            % Note: sjemea (Cotterill et al.) uses strict inequality (<)
            if all_isis(k) < maxISI_start
                inBurst = true;
                currentBurstStartIdx = k;
                % k+1 is the second spike in this interval
                % so the burst definitely includes k and k+1
                k = k + 1;
            else
                k = k + 1;
            end
        else
            % We are in a burst. Check if we continue.
            % Note: sjemea allows continuation if ISI <= maxISI_end
            % (Implied by ending if > maxISI_end)
            if all_isis(k) <= maxISI_end
                % Continue burst
                k = k + 1;
            else
                % Burst ends.
                % The last spike in the burst was k. (The interval k->k+1 was too long).
                currentBurstEndIdx = k;

                raw_bursts_idx = [raw_bursts_idx; currentBurstStartIdx, currentBurstEndIdx]; %#ok<AGROW>

                inBurst = false;
                k = k + 1;
            end
        end
    end

    % Close last burst if active
    if inBurst
        raw_bursts_idx = [raw_bursts_idx; currentBurstStartIdx, nspikes]; %#ok<AGROW>
    end

    if isempty(raw_bursts_idx)
        continue;
    end

    % --- Phase 2: Merge bursts closer than minIBI ---
    % IBI is defined as Time(Start of Burst i+1) - Time(End of Burst i)
    % Actually, Cotterill and NeuroExplorer usually define IBI as the interval between the bursts.
    % NeuroExplorer Manual 2.24: "If the interval between two bursts (Time_Start_B(i+1) - Time_End_B(i)) is less than Min Inter-Burst Interval, ... merged."

    merged_bursts_idx = [];
    if size(raw_bursts_idx, 1) > 0
        curr_start = raw_bursts_idx(1, 1);
        curr_end   = raw_bursts_idx(1, 2);

        for b = 2:size(raw_bursts_idx, 1)
            next_start = raw_bursts_idx(b, 1);
            next_end   = raw_bursts_idx(b, 2);

            % Valid IBI check
            % Time of end spike of current burst
            t_end = st(curr_end);
            % Time of start spike of next burst
            t_start = st(next_start);

            ibi = t_start - t_end;

            if ibi < minIBI
                % Merge
                curr_end = next_end;
            else
                % Commit current burst
                merged_bursts_idx = [merged_bursts_idx; curr_start, curr_end]; %#ok<AGROW>
                % Start new
                curr_start = next_start;
                curr_end   = next_end;
            end
        end
        % Commit last
        merged_bursts_idx = [merged_bursts_idx; curr_start, curr_end]; %#ok<AGROW>
    end

    % --- Phase 3: Filter by minDur and minSpks ---
    final_bursts_idx = [];
    for b = 1:size(merged_bursts_idx, 1)
        s_idx = merged_bursts_idx(b, 1);
        e_idx = merged_bursts_idx(b, 2);

        % Check spike count
        % Number of spikes = (e_idx - s_idx) + 1
        n_spks = e_idx - s_idx + 1;

        if n_spks < minSpks
            continue;
        end

        % Check duration
        dur = st(e_idx) - st(s_idx);

        if dur < minDur
            continue;
        end

        final_bursts_idx = [final_bursts_idx; s_idx, e_idx]; %#ok<AGROW>
    end

    if isempty(final_bursts_idx)
        continue;
    end

    % --- Collect Stats ---

    nb = size(final_bursts_idx, 1);

    b_struct.times = zeros(nb, 2);
    b_struct.nspks = zeros(nb, 1);
    b_struct.dur   = zeros(nb, 1);
    b_struct.freq  = zeros(nb, 1);
    b_struct.ibi   = nan(nb, 1); % First IBI is usually NaN or maybe pre

    first_last_times = zeros(nb, 2);

    for b = 1:nb
        s_idx = final_bursts_idx(b, 1);
        e_idx = final_bursts_idx(b, 2);

        b_struct.times(b, :) = [st(s_idx), st(e_idx)];
        b_struct.nspks(b)    = e_idx - s_idx + 1;
        b_struct.dur(b)      = st(e_idx) - st(s_idx);
        b_struct.freq(b)     = b_struct.nspks(b) / b_struct.dur(b);

        first_last_times(b, :) = [st(s_idx), st(e_idx)];
    end

    % Calculate IBIs (Time between end of prev and start of curr)
    if nb > 1
        b_struct.ibi(2:end) = first_last_times(2:end, 1) - first_last_times(1:end-1, 2);
    end

    % Fraction of spikes in bursts
    total_spikes_in_bursts = sum(b_struct.nspks);

    % Store in output structure
    brst.detect(iunit)  = nb;
    brst.nspks(iunit)   = mean(b_struct.nspks);
    brst.brstDur(iunit) = mean(b_struct.dur);
    brst.freq(iunit)    = mean(b_struct.freq);
    brst.bspks(iunit)   = total_spikes_in_bursts / length(st);
    brst.ibi(iunit)     = mean(b_struct.ibi, 'omitnan');

    brst.all{iunit} = b_struct;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flgSave
    save(brstfile, 'brst')
end

end
