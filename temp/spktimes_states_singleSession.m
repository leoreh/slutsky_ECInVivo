% spike timing metrics
st = spktimes_metrics('spikes', v.spikes, 'sunits', [],...
    'bins', bins, 'forceA', true, 'saveVar', true, 'fullA', false);

% brst (mea)
brst = spktimes_meaBrst(v.spikes.times, 'binsize', [], 'isiThr', 0.02,...
    'minSpks', 3, 'saveVar', true, 'force', true, 'bins', bins);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate gain ratio for various brst params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

minSpks = [2 : 2 : 10];
isiThr = [0.1];

clear b 
for ithr = 1 : length(minSpks)
    b(ithr) = spktimes_meaBrst(v.spikes.times, 'binsize', [], 'isiThr', isiThr,...
        'minSpks', minSpks(ithr), 'saveVar', false, 'force', true, 'bins', bins);
end

clear gfactor
for ithr = 1 : length(minSpks)
    for ivar = 1 : length(brstVar)
        vec1 = b(ithr).(brstVar{ivar})(1, :)';
        vec2 = b(ithr).(brstVar{ivar})(2, :)';
        gfactor.(brstVar{ivar})(:, ithr) = (vec2 - vec1) ./ sum([vec1, vec2]')';
    end
end

% -------------------------------------------------------------------------
% plot gain factor of vars as a function of brst params

fh = figure;
th = tiledlayout(1, length(brstVar), 'TileSpacing', 'Compact');
title(th, num2str(isiThr))
for ivar = 1 : length(brstVar)
    axh = nexttile;
    dataMat = gfactor.(brstVar{ivar})(unitIdx, :);
    plot_boxMean(dataMat, 'clr', 'k', 'allPnts', false)
    title(brstVar{ivar})
    ylabel('GainFactor')
    xticklabels(split(num2str(minSpks)))
    xlabel('minSpks')
end


% -------------------------------------------------------------------------
% plot mean value of var per state, for each brst param

fh = figure;
th = tiledlayout(1, length(minSpks) + 1, 'TileSpacing', 'Compact');

% mfr
axh = nexttile;
plot_boxMean(mfrStates(unitIdx, :), 'clr', 'k', 'allPnts', false)
set(gca, 'yscale', 'log')
ylabel('MFR [Hz]')
xticklabels(stateNames(sstates))

% brst var 
ivar = 2;
for ithr = 1 : length(minSpks)
    axh = nexttile;
    dataMat = b(ithr).(brstVar{ivar})(:, unitIdx)';
    plot_boxMean(dataMat, 'clr', 'k', 'allPnts', false)
    set(gca, 'yscale', 'log')
    ylabel(brstVar{ivar})
    title(num2str(minSpks(ithr)))
    xticklabels(stateNames(sstates))
end


%%%
fh = figure;
% plot(mfrStates(:, 1), mfrStates(:, 2), '.', 'MarkerSize', 20)
plot(mfrStates(:, 2), v.fr.states.gain(4, :), '.', 'MarkerSize', 20)


% -------------------------------------------------------------------------
% predicition: prcnt spks in brsts should be correlated w/ mfr gain factor
fh = figure;
th = tiledlayout(2, 2, 'TileSpacing', 'Compact');

xval = b(1).spkprct(2, unitIdx);

axh = nexttile;
yval = mfrStates(unitIdx, 1);
plot(xval, yval, '.', 'MarkerSize', 20)
xlim([0 100])
xlabel('Spks In Brsts [%]')
ylim([0 2])
ylabel('MFR in AW')
p = polyfit(xval, yval, 1);
yfit = polyval(p, xval);
hold on
plot(xval, yfit, '--r')

axh = nexttile;
yval = mfrStates(unitIdx, 2);
plot(xval, yval, '.', 'MarkerSize', 20)
xlim([0 100])
xlabel('Spks In Brsts [%]')
ylim([0 2])
ylabel('MFR in NREM')
p = polyfit(xval, yval, 1);
yfit = polyval(p, xval);
hold on
plot(xval, yfit, '--r')

axh = nexttile;
axh.Layout.TileSpan = [1, 2];
yval = v.fr.states.gain(4, unitIdx);
plot(xval, yval, '.', 'MarkerSize', 20)
xlim([0 100])
xlabel('Spks In Brsts [%]')
ylabel('MFR gain')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spktimes metrics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sstates = [1, 4, 5];
stateNames = v.ss.info.names;
bins = v.ss.stateEpochs(sstates);

stVar = ["royer"; "doublets"; "lidor"; "mizuseki"; "cv2"];

clear gfactor
for ivar = 1 : length(stVar)
    vec1 = st.(stVar{ivar})(1, :)';
    vec2 = st.(stVar{ivar})(2, :)';
    gfactor.(stVar{ivar}) = (vec2 - vec1) ./ sum([vec1, vec2]')';
end


fh = figure;
th = tiledlayout(1, 2, 'TileSpacing', 'Compact');
axh = nexttile;
plot_boxMean(v.fr.states.gain(4, unitIdx)', 'clr', 'k', 'allPnts', true)
ylabel('GainFactor')

axh = nexttile;
dataMat = [cell2nanmat(struct2cell(gfactor), 2)];
plot_boxMean(dataMat(unitIdx, :), 'clr', 'k', 'allPnts', true)
ylabel('GainFactor')
xticklabels(stVar)


ivar = 3;
fh = figure;
plot(st.(stVar{1})(1, unitIdx)', v.fr.states.gain(4, unitIdx)', '.', 'MarkerSize', 20)







% -------------------------------------------------------------------------

% mfr in states
mfrStates = cellfun(@(x) mean(x, 2), v.fr.states.fr, 'uni', false);
mfrStates = cell2mat(mfrStates(sstates));


% get mfr in states
mfr = cellfun(@(x) mean(x, 2), fr.states.fr, 'uni', false);
clear mfrStates
for istate = 1 : length(sstates)
    mfrStates(:, istate) = vertcat(mfr{:, sstates(istate)});
end

% select units
brsty = mean(st.lidor, 1, 'omitnan');
unitBrsty = brsty > prctile(brsty(units.rs), 50);
unitMfr = fr.mfr > prctile(fr.mfr(units.rs), 50);
unitIdx = unitBrsty & units.rs;

% plot
fh = figure;
th = tiledlayout(1, 2, 'TileSpacing', 'Compact');

% mfr of all rs units
axh = nexttile;
plot_boxMean(mfrStates(units.rs, :), 'clr', 'k', 'allPnts', false)
set(gca, 'yscale', 'log')
ylabel('MFR [Hz]')
xticklabels(stateNames(sstates))

% mfr of only bursty units
axh = nexttile;
ivar = 1;
plot_boxMean(mfrStates(unitIdx, :), 'clr', 'k', 'allPnts', false)
set(gca, 'yscale', 'log')
ylabel('MFR [Hz]')
xticklabels(stateNames(sstates))

