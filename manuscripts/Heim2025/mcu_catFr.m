

grp = 'wt';
mice = [mcu_sessions(grp)];
vars = {'fr'; 'units'};

flg_plot = true;
pertWin = [1 : 3000];

for iMouse = 1 : nMice

    basepaths = mcu_sessions(mice{iMouse});
    nfiles = length(basepaths);
    
    % load vars
    v = basepaths2vars('basepaths', basepaths, 'vars', vars);
    
    % cat units
    units = catfields([v(:).units], 'addim');
    rs = squeeze(units.clean(1, :, :));
    fs = squeeze(units.clean(1, :, :));

    % cat fr
    fr = catfields([v(:).fr], 'cell');
    frt = cell2padmat(fr.strd, 2);
    t = [1 : length(frt)] * 60;         % hr

    % find perturbation onset from mfr
    mfr = mean(frt, 1, 'omitnan');
    mfr = mea_frDenoise(mfr, t, 'flgPlot', false, 'frameLenSec', 600);
    [pertOnset, ~] = findchangepts(mfr(pertWin), ...
        'Statistic', 'mean', 'minDistance', 5,...
        'MaxNumChanges', 1);
    
    % Visualization to verify
    if flg_plot
        figure;
        plot(mfr);
        hold on;
        xline(pertOnset, 'r-', 'LineWidth', 2);
        title('Detected Perturbation Onset');
    end
    
    % Store data
    fr_mouse{iMouse} = frt;
    t_mouse{iMouse} = t - pertOnset;

end

