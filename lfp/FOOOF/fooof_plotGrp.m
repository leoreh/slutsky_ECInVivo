function fh = fooof_plotGrp(f1f, dimAvg, dimGrp, dimRow)
% Plots group-level FOOOF results with flexible dimension handling.
% Creates a tiled layout: 3 columns (PSD original, aperiodic, residual).
% Rows are determined by 'dimRow', groups within tiles by 'dimGrp'.
% Data for plot_stdShade is averaged over 'dimAvg'.
%
% fh = fooof_plotGrp(f1f, dimAvg, dimGrp, dimRow)
%
% INPUTS:
%   f1f      - struct containing FOOOF results. Expected fields:
%              f1f.freqs    - [1 x nFreqs] vector of frequencies.
%              f1f.psd_orig - [..., nFreqs] matrix of original PSDs.
%              f1f.psd_ap   - [..., nFreqs] matrix of aperiodic fits.
%                            The dimensions before nFreqs can be arbitrary.
%   dimAvg   - Scalar, specifies the dimension index to average over for plot_stdShade.
%   dimGrp   - Scalar, specifies the dimension index for grouping plots within each tile.
%   dimRow   - Scalar, specifies the dimension index for creating rows in the tiled layout.
%
% OUTPUTS:
%   fh - handle to the generated figure.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataSize = size(f1f.psd_orig);
nDims = length(dataSize);

if nDims < max([dimAvg, dimGrp, dimRow, nDims])
    error('Input data has fewer dimensions than specified by dimAvg, dimGrp, or dimRow combined with the frequency dimension.');
end

dFreq = nDims; % Frequency dimension is assumed to be the last one

% Robustly extract the frequency vector
freqs_raw = f1f.freqs;
freqs_size = size(freqs_raw);
freqs_nDims = length(freqs_size);
freq_slicer = cell(1, freqs_nDims);
for kDim = 1:(freqs_nDims-1)
    freq_slicer{kDim} = 1;
end
freq_slicer{freqs_nDims} = ':';
freqs = squeeze(freqs_raw(freq_slicer{:}));

% If freqs happens to be a column vector after slicing and squeezing, ensure it's a row vector
if iscolumn(freqs)
    freqs = freqs';
end

nFreqs = dataSize(dFreq);

if length(freqs) ~= nFreqs
    error('Length of f1f.freqs does not match the size of the frequency dimension in psd_orig.');
end

nRows = dataSize(dimRow);
nGrps = dataSize(dimGrp);

% Calculate psd_res: log10(original PSD) - log10(aperiodic component)
f1f.psd_res = log10(f1f.psd_orig) - log10(f1f.psd_ap);
f1f.psd_res(isinf(f1f.psd_res)) = NaN; 

psdFields = {'psd_orig', 'psd_ap', 'psd_res'};
colTitles = {'Original PSD', 'Aperiodic Fit', 'Flattened Spectrum (Residuals)'};
yScales = {'log', 'log', 'linear'}; 
yLabels = {'Power', 'Power', 'Power (log10)'};

% Define a set of distinguishable colors for groups
grpClrs = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], ...
           [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], ...
           [0.6350 0.0780 0.1840]}; 
if nGrps > length(grpClrs)
    grpClrs = repmat(grpClrs, 1, ceil(nGrps / length(grpClrs)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fh = figure;
set(fh, 'WindowState', 'maximized', 'Name', 'FOOOF Generalized Group Analysis');
set(fh, 'DefaultAxesFontSize', 12);

tileDims = [nRows, 3]; 
th = tiledlayout(tileDims(1), tileDims(2));
th.TileSpacing = 'compact';
th.Padding = 'compact';

for iRow = 1 : nRows
    for iCol = 1 : 3 % Iterate through psd_orig, psd_ap, psd_res
        
        axh = nexttile(th); 
        cla(axh); 
        hold(axh, 'on');
        
        psdFieldName = psdFields{iCol};
        psdData = f1f.(psdFieldName);
        
        yScale = yScales{iCol};
        yLabel = yLabels{iCol};
        
        lgdHndls = [];
        lgdEntries = {};
        
        for iGrp = 1 : nGrps
            slicer = cell(1, nDims);
            for kDim = 1:nDims
                slicer{kDim} = ':'; 
            end
            
            slicer{dimRow} = iRow;
            slicer{dimGrp} = iGrp;
            
            otherDimsIdx = setdiff(1:nDims, [dFreq, dimAvg, dimRow, dimGrp]);
            for kOther = otherDimsIdx
                slicer{kOther} = 1; 
            end
            
            slicedCube = psdData(slicer{:});
            
            permOrder = 1:nDims;
            idxDFreq = find(permOrder == dFreq); permOrder(idxDFreq) = [];
            idxDimAvg = find(permOrder == dimAvg); permOrder(idxDimAvg) = [];
            finalPermOrder = [dFreq, dimAvg, permOrder];
            
            permutedData = permute(slicedCube, finalPermOrder);
            
            stdShadeData = squeeze(permutedData);
            
            nObsAvg = dataSize(dimAvg);
            if size(stdShadeData, 1) == nObsAvg && size(stdShadeData, 2) == nFreqs
                 if nFreqs == 1 && nObsAvg ~=1 
                 elseif nObsAvg ~= nFreqs 
                    stdShadeData = stdShadeData';
                 end
            end
            stdShadeData = reshape(stdShadeData, [nFreqs, nObsAvg]);

            if ~isempty(stdShadeData)
                ph = plot_stdShade('dataMat', stdShadeData', 'xVal', freqs, 'axh', axh, ...
                                   'clr', grpClrs{iGrp}, 'alpha', 0.3);
                if ~isempty(ph) 
                    lgdHndls(end+1) = ph;
                    lgdEntries{end+1} = sprintf('Group %d', iGrp);
                end
            end
        end
        
        set(axh, 'xscale', 'log', 'yscale', yScale);
        xlabel(axh, 'Frequency (Hz)');
        ylabel(axh, yLabel);
        
        tileTitle = sprintf('Row %d: %s', iRow, colTitles{iCol});
        title(axh, tileTitle, 'Interpreter', 'none');
        
        if ~isempty(lgdHndls)
            legend(axh, lgdHndls, lgdEntries, 'Location', 'northeast', 'Interpreter', 'none');
        end
        grid(axh, 'on');
    end
end

figTitle = 'FOOOF Generalized Group Analysis';
title(th, figTitle, 'interpreter', 'none', 'FontSize', 16);

end
% EOF 