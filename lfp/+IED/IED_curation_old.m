function IED_curation_old(varargin)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'sig', [], @isnumeric)
addParameter(p, 'fs', 1250, @isnumeric)
addParameter(p, 'thr', [10 0], @isnumeric)
addParameter(p, 'basepath', pwd, @isstr);
addParameter(p, 'basename', [], @isstr);
addParameter(p, 'saveVar', true, @islogical);

parse(p, varargin{:})
sig = p.Results.sig;
fs = p.Results.fs;
thr = p.Results.thr;
basepath = p.Results.basepath;
basename = p.Results.basename;
saveVar = p.Results.saveVar;

% params
tstamps = [1 : length(sig)]' / fs;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% manual curation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fh = figure;
    set(fh, 'units','normalized','outerposition',[0 0 1 1]);
    suptitle({'Inspect IIS: left = accept; right = decline; middle = previous'...
        'Wait for curser before pressing'});
    sp1 = subplot(2, 1, 1);
    ylabel('Voltage [mV]')
    set(gca, 'TickLength', [0 0])
    box off
    sp2 = subplot(2, 1, 2);
    i = 1;
    marg1 = fs * 5;
    marg2 = fs * 0.5;
    yrange = [min(sig) max(sig)];
    while i <= nspks
        idx1(1) = max([1 pos(i) - marg1]);
        idx1(2) = min([length(sig) pos(i) + marg1]);
        idx2(1) = max([1 pos(i) - marg2]);
        idx2(2) = min([length(sig) pos(i) + marg2]);
        hold(sp1, 'off')
        plot(sp1, tstamps(idx1(1) : idx1(2)), sig(idx1(1) : idx1(2)))
        hold(sp1, 'on')
        axis tight
        plot(sp1, [idx1(1) idx1(2)] / fs, [thr(2) thr(2)], '--r')
        plot(sp1, [idx1(1) idx1(2)] / fs, -[thr(2) thr(2)], '--r')
        scatter(sp1, pos(i) / fs, peak(i), '*');
        axis(sp1, 'tight')
        ylim(sp1, yrange)
        ylabel(sp1, 'Voltage [mV]')
        set(sp1, 'TickLength', [0 0])
        box off
        
        hold(sp2, 'off')
        plot(sp2, tstamps(idx2(1) : idx2(2)), sig(idx2(1) : idx2(2)))
        hold(sp2, 'on')
        plot(sp2, [idx2(1) idx2(2)] / fs, [thr(2) thr(2)], '--r')
        plot(sp2, [idx1(1) idx1(2)] / fs, -[thr(2) thr(2)], '--r')
        scatter(sp2, pos(i) / fs, peak(i), '*');
        ylim(sp2, [min(sig(idx2(1) : idx2(2))) max(sig(idx2(1) : idx2(2)))])
        xlim(sp2, [idx2(1) idx2(2)] / fs)
        ylabel(sp2, 'Voltage [mV]')
        xlabel(sp2, 'Time [s]')
        set(sp2, 'TickLength', [0 0])
        box off
        title(sp2, sprintf('spk %d/%d', i, nspks))
        while 1
            [~, ~, button] = ginput(1);
            if button == 1
                iis.accepted(i) = 1;
                break
            elseif button == 3
                iis.accepted(i) = 0;
                break
            elseif button == 2
                i = i - 2;
                break
            end
            prevButton = button;
        end
        i = i + 1;
    end
    close(fh)
        
    seg(~iis.accepted, :) = [];
    peak(~iis.accepted) = [];
    pos(~iis.accepted) = [];
    nspks = length(pos);
    fprintf('after manual curation: %d IIS\n', nspks);    
    if nspks == 0
        fprintf('no spikes detected\n');
        return
    end   

end