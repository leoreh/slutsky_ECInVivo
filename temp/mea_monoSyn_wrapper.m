
mono_res = ce_MonoSynConvClick(spikes,'includeInhibitoryConnections',true/false); % detects the monosynaptic connections

spikes.times = mea.spktimes;
spikes.shankID = mea.ch;
spikes.cluID = [1 : length(mea.spktimes)];


mono_res = ce_MonoSynConvClick(spikes,...
    'includeInhibitoryConnections', true, 'epoch', [0 Inf]);

basepath = pwd;
[~, basename] = fileparts(basepath);
monofile = fullfile(basepath, [basename, '.mono_res.cellinfo.mat']);
save(monofile, 'mono_res');
setMatlabGraphics(true)
gui_MonoSyn(monofile);
load(monofile)
mono_res

units = [30, 45];

clr = 'rb';
plot_monoSyn('spktimes', mea.spktimes, 'swv', swv, 'units', units,...
    'clr', clr, 'fs', 10000)
