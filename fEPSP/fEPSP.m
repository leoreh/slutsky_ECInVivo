
% this is a wrapper for fEPSP analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basepath = 'E:\data\fEPSP\Ortal';
graphics = true;
saveFig = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nsessions = 2;
amp = stability('nsessions', nsessions, 'inspect', true, 'basepath', basepath,...
    'graphics', graphics, 'saveFig', saveFig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% arrange data to load.
% make sure order of files is according to stimulus intensity
intensity = [0.2 : 0.1 : 0.8, 1];
filename{1} = 'IO_1AP_0.06Hz_0.5msduration_0.2mA';
filename{2} = 'IO_1AP_0.06Hz_0.5msduration_0.3mA';
filename{3} = 'IO_1AP_0.06Hz_0.5msduration_0.4mA';
filename{4} = 'IO_1AP_0.06Hz_0.5msduration_0.5mA';
filename{5} = 'IO_1AP_0.06Hz_0.5msduration_0.6mA';
filename{6} = 'IO_1AP_0.06Hz_0.5msduration_0.7mA';
filename{7} = 'IO_1AP_0.06Hz_0.5msduration_0.8mA';
filename{8} = 'IO_1AP_0.06Hz_0.5msduration_1mA';

pk_io = io('filename', filename, 'intensity', intensity, 'inspect', false,...
 'basepath', basepath, 'graphics', graphics, 'saveFig', saveFig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
intensity = [0.4 0.8];
clear filename
filename{1} = 'STP2_5AP_50Hz_0.5msduration_0.4mA';
filename{2} = 'STP2_5AP_50Hz_0.5msduration_0.8mA';

pk_stp = stp('filename', filename, 'intensity', intensity, 'inspect', false,...
 'basepath', basepath, 'graphics', graphics, 'nstim', 5, 'saveFig', saveFig);
