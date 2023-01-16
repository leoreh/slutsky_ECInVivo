function [snips, stamps] = snipFromBinary(varargin)

% maps a binary file to memory and snips segments sorrounding specific
% samples. the file can be .dat, .temp_wh, .lfp etc. user defines the
% window length sorrounding each snippet (does not have to be symmetrical).
% stamps can be a single vector or a cell of vectors. snips will be
% organized in a cell of 3d mats accordingly (ch x sampels x stamps). can
% detrend each snip based on specific segments, perfrom L2 noramlization,
% and re-snip such that the peak is at the center of the snip. 
%
% INPUT:
%   fname       string. full path, name and extension of binary file to snip from.
%               exists in basepath or if fname can be extracted from basepath
%   stamps      vec / cell. pointers to dat files from which snipping will occur
%               [samples]. if cell than snips will also be organized in a cell.
%               this is to avoid the need to call the function (and memmap)
%               multiple times when, e.g. organizing clusters.
%   raw         data field of memory map to binary file. this can be helpful
%               when calling the function several times to keep snips relatively
%               small
%   win         vec of 2 elements. determines length of snip. for example,
%               win = [5 405] each snip will be 401 samples, starting
%               5 samples after the corresponding stamp and ending 405
%               samples after stamp. if win = [-16 16] than snip will be of
%               33 samples symmetrically centered around stamp.
%   nchans      numeric. number of channels in dat file {35}.
%   ch          vec. channels to load from dat file {[]}. if empty than
%               all channels will be loaded. can be a cell of equal length to
%               stamps such that each cell specifies the channels to load for
%               the corresponding stamps.
%   align_peak  char. align snip such that the peak in at the center. can
%               be 'min' or 'max' which will determine the polarity of the
%               peak
%   precision   char. sample precision of binary file {'int16'}
%   rmv_trend   numeric. used to detrend each snippet.
%               the trend can be determined based on specific samples. for
%               example, if rmv_trend = [8 10], the trend will be calculated
%               based on the first 8 and last 10 samples. if
%               rmv_trend = 0 than no trending will occur. if rmv_trend >
%               diff(win) than the entire snip will be used to find the
%               linear fit.
%               than the trend will be calculated on the entire snippet.
%   l2norm      logical. L2 normalize snippets {false}.
%   saveVar     logical. save snips or not. name of file will be taken from
%               path of fname.
%
% OUTPUT
%   snips       array of mats ch x sampels x stamps. class double
%   stamps      same as input but with correction due to alignment
% 
% CALLS:
%   class2bytes
%
% TO DO LIST:
%
% 10 apr 20 LH  updates:
% 13 may 20 LH      separate detrend and normalize
% 12 dec 20 LH      detrend based on segment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'fname', '', @ischar);
addOptional(p, 'stamps', []);
addOptional(p, 'raw', []);
addOptional(p, 'win', [-16 16], @isnumeric);
addOptional(p, 'nchans', 35, @isnumeric);
addOptional(p, 'ch', []);
addOptional(p, 'align_peak', 'min', @ischar);
addOptional(p, 'precision', 'int16', @ischar);
addOptional(p, 'rmv_trend', [6], @isnumeric);
addOptional(p, 'l2norm', false, @islogical);
addOptional(p, 'saveVar', false, @islogical);

parse(p, varargin{:})
fname       = p.Results.fname;
stamps      = p.Results.stamps;
raw         = p.Results.raw;
win         = p.Results.win;
nchans      = p.Results.nchans;
ch          = p.Results.ch;
align_peak  = p.Results.align_peak;
precision   = p.Results.precision;
rmv_trend   = p.Results.rmv_trend;
l2norm      = p.Results.l2norm;
saveVar     = p.Results.saveVar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~iscell(stamps)
    stamps = {stamps};
end
ncells = length(stamps);

% handle channels
if isempty(ch)
    cell(ncells, 1);
    for icell = 1 : ncells
        ch{icell} = 1 : nchans;
    end
end
if ~iscell(ch)
    ch = {ch};
end

% size of one data point in bytes
nbytes = class2bytes(precision);

% -------------------------------------------------------------------------
% build regressor to find the trend of the snip using the first and last
% samples as specified in rmv_trend. this is used, for example, so a spike
% in the center of a snip does not influence the linear fit. Based on
% matrix form of polynomial regression.
sniplength = diff(win) + 1;
psamp = floor(sniplength / 2);

% find segements of snip to use when finding the trend
if length(rmv_trend) > 2
    error('rmv_trend should be a 2-element vector')
elseif length(rmv_trend) == 1
    rmv_trend = repmat(rmv_trend, 1, 2);
end
if sum(rmv_trend) > sniplength 
    use_idx = [1 : sniplength];
elseif sum(rmv_trend) > 0
    use_idx = [1 : rmv_trend(1), sniplength : -1 : sniplength - rmv_trend(2) + 1];
else
    use_idx = [];
end

% build regressor
reg_len = length(use_idx);
w = [0 : reg_len - 1] / (reg_len - 1);
w = [reshape(w, reg_len, []), ones(reg_len, 1)];
[q, r] = qr(w, 0);
w_full = [0 : sniplength - 1] / (sniplength - 1);
w_full = [reshape(w_full, sniplength, []), ones(sniplength, 1)];

% -------------------------------------------------------------------------
% map binary file
minput = false;
if isempty(raw)
    [basepath, basename] = fileparts(fname);
    info = dir(fname);
    if isempty(info)
        error('%s does not exist', fname)
    end
    nsamps = info.bytes / nbytes / nchans;
    m = memmapfile(fname, 'Format', {precision, [nchans, nsamps] 'mapped'});
    raw = m.Data;
else
    minput = true;
end
nsamps = size(raw.mapped, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% snip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nextracting from %s\n', fname)

% go over stamps and snip data
for icell = 1 : length(stamps)
    
    fprintf('\nworking on cell #%d: ', icell)
    
    % initialize
    rmvsnip = [];
    nsnips = length(stamps{icell});
    snips{icell} = nan(length(ch{icell}), sniplength, nsnips);
    
    for isnip = 1 : nsnips
        
         % print progress
        if mod(isnip, 10000) == 0
            if isnip ~= 10000
                fprintf(repmat('\b', 1, length(txt)))
            end
            txt = ['Extracted ', num2str(isnip), ' / ', num2str(nsnips), ' snips'];
            fprintf('%s', txt)     
        end     
        
        % fix special case where snip is at begining / end of recording
        if stamps{icell}(isnip) + win(1) < 1 ||...
                stamps{icell}(isnip) + win(2) > nsamps
            warning('\nskipping stamp %d because snip incomplete', isnip)
            rmvsnip = [rmvsnip, isnip];
            continue
        end

        % snip

        v = double(raw.mapped(ch{icell}, stamps{icell}(isnip) + win(1) :...
            stamps{icell}(isnip) + win(2)));
        
        % find channel with max amplitude 
        [~, ampMaxCh] = max(range(v, 2));
        
        % realign according to min / max
        if strcmp(align_peak, {'min', 'max'})
            switch align_peak
                case 'min'
                    [~, ia] = min(v, [], 2);
                case 'max'
                    [~, ia] = max(v, [], 2);
                otherwise
            end
            wvpeak = ia(ampMaxCh) ;
            snipshift = wvpeak - psamp;
            if snipshift ~= 0
                stamps{icell}(isnip) = stamps{icell}(isnip) + snipshift;
                v = double(raw.mapped(ch{icell}, stamps{icell}(isnip) + win(1) :...
                    stamps{icell}(isnip) + win(2)));
            end
        else
            ia = psamp;
        end

        
        % find and remove trend
        if sum(rmv_trend) > 0
            trend = w_full * (r \ q' * v(:, use_idx)');
            v = v - trend';
        end
        
        % L2 normalize
        if l2norm
            v = v ./ vecnorm(v, 2, 2);
        end
               
        snips{icell}(:, :, isnip) = v;
                        
    end

    % remove incomplete snips
    for isnip = 1 : length(rmvsnip)
        snips{icell}(:, :, rmvsnip(isnip)) = [];
        stamps{icell}(rmvsnip(isnip)) = []
    end


end

% return outputs as vec / mat instead of cell 
if length(stamps) == 1
    stamps = cell2mat(stamps);
    snips = cell2mat(snips);
end

clear raw
if minput
    clear m
end

% save
if saveVar
    savename = fullfile(basepath, [basename, '.snips']);
end

end

% EOF