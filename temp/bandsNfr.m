


mname = 'lh107';
[bands, tbands] = sessions_catVarTime('mname', mname,...
    'dataPreset', 'bands', 'graphics', false,...
    'basepaths', {}, 'xTicksBinsize', 6, 'markRecTrans', true);

[fr, tfr] = sessions_catVarTime('mname', mname,...
    'dataPreset', 'fr', 'graphics', false,...
    'basepaths', {}, 'xTicksBinsize', 6, 'markRecTrans', true);

fs = 20 * 60;
fsRat = fs / unique(diff(tbands));
bandsN = downsample(bands, fsRat);
tbandsN = downsample(tbands, fsRat) / 60 / 60;

fsRat = fs / unique(diff(tfr));
frN = downsample(fr, fsRat);
tfrN = downsample(tfr, fsRat) / 60 / 60;

fh = figure;
th = tiledlayout(2, 1, 'TileSpacing', 'Compact');
axh = nexttile;
plot(tfrN, frN)

axh = nexttile;
plot(tbandsN, bandsN)