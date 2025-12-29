try
    % Mock Data
    % Create 2 files, each with 2 units
    v = struct();

    % File 1
    spk1 = sort(rand(100, 1) * 1000); % Unit 1
    spk2 = sort(rand(150, 1) * 1000); % Unit 2
    v(1).mea.spktimes = {spk1, spk2};
    v(1).fr.info.group = 'Control';

    % File 2
    spk3 = sort(rand(50, 1) * 1000);
    spk4 = sort(rand(500, 1) * 1000); % High rate
    v(2).mea.spktimes = {spk3, spk4};
    v(2).fr.info.group = 'KO';

    % Parameters
    winLim = [0, 1000];
    winBsl = [0, 200];
    winSs  = [800, 1000];

    fprintf('Running spk2ca_rcv with mock data...\n');

    % Call Function
    [tbl, hFig] = spk2ca_rcv(v, ...
        'winLim', winLim, ...
        'winBsl', winBsl, ...
        'winSs', winSs, ...
        'flgPlot', true);

    % Checks
    assert(height(tbl) == 4, 'Table height %d (expected 4)', height(tbl));
    assert(ismember('Mito_LogRatio', tbl.Properties.VariableNames), 'Missing Mito_LogRatio');
    assert(ismember('Rate_LogRatio', tbl.Properties.VariableNames), 'Missing Rate_LogRatio');
    assert(isa(hFig, 'matlab.ui.Figure'), 'hFig is not a figure handle');

    fprintf('VERIFICATION PASSED.\n');

    % Optional: Inspect table
    disp(tbl(1, {'Group', 'Rate_Bsl', 'Rate_Ss', 'Rate_LogRatio'}));

    close(hFig);

catch ME
    fprintf('VERIFICATION FAILED: %s\n', ME.message);
    fprintf('Stack:\n');
    disp(ME.stack(1));
end
