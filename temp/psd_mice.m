% psd_mice

bandNames = ["broad", "swa", "delta", "theta", "alpha", "beta", "gamma"];
mice = ["lh96"; "lh107"; "lh105"];
datapath = 'F:\Data\';

for imouse = 1 : length(mice)
       
    mname = mice{imouse};
    
    % analyze
    bands = sessions_psd(mname, 'flgNormBand', false, 'flgAnalyze', false,...
        'flgNormTime', true, 'flgEmg', true, 'idxBsl', [1], 'flgDb', false);
    
    % load
    bandsfile = fullfile(datapath, mname, [mname, '_psdBands.mat']);
    v(imouse).bands = load(bandsfile, 'bands');
end

% cat bands across mice. array is freqBand x state x session x mouse
bands = catfields([v.bands], 'catdef', 'addim');

% organize to prism
sbands = [2, 3, 4, 7];
nstates = size(bands.bands, 2);
for istate = 1 : nstates
    for iband = 1 : length(sbands)
        bcell{istate} = squeeze(bands.bands(sbands(iband), istate, :, :));
    end
end




mname = 'lh111';
frTiles = sessions_frStates(mname, 'flgAnalyze', true,...
    'flgEmg', false, 'ntiles', 2);




