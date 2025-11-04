
mname = {'lh122'; 'lh123'; 'lh126'; 'lh129'; 'lh130'};
nmice = length(mname);

iUnit = 1;

% Perturbation time from manual inspection
pertHW = 4 * 60;        % Half window for clipping recording
pertTime(1) = 1680;
pertTime(2) = 1850;
pertTime(3) = 1730;
pertTime(4) = 1950;
pertTime(5) = 2390;

for iMouse = 1 : nmice

    queryStr = [mname{iMouse}, '_mk801'];
    basepaths = mk801_sessions(queryStr);

    vars = {'fr', 'units'};
    v = basepaths2vars('basepaths', basepaths, 'vars', vars);
    nPaths = length(basepaths);
    
    clear yData
    for iPath = 1 : nPaths
        fr = v(iPath).fr;
        units = v(iPath).units;
        uIdx = units.clean(iUnit, :);
        yData{iPath} = fr.strd(uIdx, :);                             
    end
    yData = cell2padmat(yData, 2);
    
    % manually inspect pert time [sample]
    % xData = [1 : length(yData)];
    % hFig = figure;
    % plot(xData, mean(yData, 1, 'omitnan'))
    % hAx = gca;
    % hAx.XAxis.Exponent = 0;

    recWin = pertTime(iMouse) - 6 * 60 : pertTime(iMouse) + 48 * 60;
    frCell{iMouse} = yData(:, recWin);

end

frMat = cell2padmat(frCell, 1)';
xData = ([1 : length(frMat)] / 60) - 6;

hFig = figure;
plot(mean(yData, 1, 'omitnan'))


pertWin = [24, 72] * 60;
min(mean(mData(pertWin(1) : pertWin(2)), 1, 'omitnan'))