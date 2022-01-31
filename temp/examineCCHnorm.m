ccg_bnsz = 0.001;
ccg_dur = 0.05;
[ccg50, ccg50_tstamps] = CCG(mea.spktimes,...
    [], 'binSize', ccg_bnsz,...
    'duration', ccg_dur, 'Fs', 1 / mea.info.fs);
[ccg50norm, ~] = CCG(mea.spktimes,...
    [], 'binSize', ccg_bnsz,...
    'duration', ccg_dur, 'Fs', 1 / mea.info.fs, 'norm', 'rate');

fh = figure;
setMatlabGraphics(false)
targets = [12, 105, 20, 32];
ncells = length(targets);
reference = 19;
sbidx = reshape([1 : ncells * 4], ncells, 4);
cnt = 1;
for icell = 1 : ncells
    subplot(ncells, 4, sbidx(cnt))
    bh = bar(ccg50_tstamps,...
        squeeze(ccg50(:, reference, targets(icell))), 'BarWidth', 1);
    hold on
    plot([0, 0], ylim, '--k')
    ylabel('Counts')
    xlabel('Time [ms]')
    bh.FaceColor = 'k';
    bh.EdgeColor = 'none';
    box off
    axis tight
    title(sprintf('%d => %d', reference, targets(icell)))
    cnt = cnt + 1;
    
    subplot(ncells, 4, sbidx(cnt))
    bh = bar(ccg50_tstamps,...
        squeeze(ccg50norm(:, reference, targets(icell))), 'BarWidth', 1);
    hold on
    plot([0, 0], ylim, '--k')
    ylabel('Rate [Hz]')
    xlabel('Time [ms]')
    bh.FaceColor = 'k';
    bh.EdgeColor = 'none';
    box off
    axis tight
    title(sprintf('%d => %d - norm', reference, targets(icell)))
    cnt = cnt + 1;
    ylim([0 20])
    
    subplot(ncells, 4, sbidx(cnt))
    bh = bar(ccg50_tstamps,...
        squeeze(ccg50(:, targets(icell), reference)), 'BarWidth', 1);
    hold on
    plot([0, 0], ylim, '--k')
    ylabel('Counts')
    xlabel('Time [ms]')
    bh.FaceColor = 'k';
    bh.EdgeColor = 'none';
    box off
    axis tight
    title(sprintf('%d => %d', targets(icell), reference))
    cnt = cnt + 1;
    
    subplot(ncells, 4, sbidx(cnt))
    bh = bar(ccg50_tstamps,...
        squeeze(ccg50norm(:, targets(icell), reference)), 'BarWidth', 1);
    hold on
    plot([0, 0], ylim, '--k')
    ylabel('Rate [Hz]')
    xlabel('Time [ms]')
    bh.FaceColor = 'k';
    bh.EdgeColor = 'none';
    box off
    axis tight
    title(sprintf('%d => %d - norm', targets(icell), reference))
    cnt = cnt + 1;
    ylim([0 20])
end
