% mcu_stg

mname = 'mcu_bsl';
vars = ["spikes"; "units"; "session"; "swv_metrics"];
[basepaths, v] = mcu_sessions(mname, vars);
nfiles = length(basepaths);

% mono synaptic interactions (spike transmission gain)
for ifile = 1 : nfiles
    
    % recording params   
    basepath = basepaths{ifile};
    [~, basename] = fileparts(basepath);
    cd(basepath)
    fs = v(ifile).session.extracellular.sr;

    stg = calc_stgs('spktimes', v(ifile).spikes.times, 'basepath', pwd,...
        'winCalc', [0, Inf], 'saveVar', true, 'mancur', false,...
        'forceA', true, 'fs', fs);
    
    mancur_stg('stg', stg)
    
    % cell explorer
    mono_res = ce_MonoSynConvClick(v(ifile).spikes,...
        'includeInhibitoryConnections', true, 'bout', [0 4 * 60 * 60]);


    monofile = fullfile(basepath, [basename, '.mono_res.cellinfo.mat']);
    save(monofile, 'mono_res');
    gui_MonoSyn(monofile);
    load(monofile)
    mono_res

    cpair = [93, 54];
    
end