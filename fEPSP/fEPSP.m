
% this is a wrapper for fEPSP analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basepath = 'E:\Data\fEPSP\Daniel\18nov19';
graphics = true;
saveFig = false;
saveVar = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE: load file via GUI, remove DC, and plot specific traces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = import_wcp();

% find artifact onset from derirative and remove DC
[~, art] = max(diff(data.S(:, 1)));
data.S = rmDC(data.S, [1, art - 0.003 * data.fs]);

trace_select = [1 : 5];

[~, ~] = rmTraces(data.S(:, trace_select), 'x', data.T');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nsessions = 4;
amp = stability('nsessions', nsessions, 'inspect', false, 'basepath', basepath,...
    'graphics', graphics, 'saveFig', saveFig, 'pkMet', 'avgMin');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% arrange data to load.
% make sure order of files is according to stimulus intensity
intensity = [0.2 : 0.1 : 0.8];

% alt 1 to arrange filename if named sequentially
% firstfile = 5; 
% prefix = '230319';
% for i =  1 : length(intensity)
    % if firstfile - 1 + i < 10
      %   filename{i} = sprintf('%s_00%d', prefix, firstfile - 1 + i);
    % else
      %   filename{i} = sprintf('%s_0%d', prefix, firstfile - 1 + i);
    % end
% end
% alt 2 to arrange filename manually
filename{1} = 'IO_1AP_0.06Hz_0.5msduration_0.2mA';
filename{2} = 'IO_1AP_0.06Hz_0.5msduration_0.3mA';
filename{3} = 'IO_1AP_0.06Hz_0.5msduration_0.4mA';
filename{4} = 'IO_1AP_0.06Hz_0.5msduration_0.5mA';
filename{5} = 'IO_1AP_0.06Hz_0.5msduration_0.6mA';
filename{6} = 'IO_1AP_0.06Hz_0.5msduration_0.7mA';
filename{7} = 'IO_1AP_0.06Hz_0.5msduration_0.8mA';
% filename{8} = '230319_012';
% filename{9} = '230319_013';
% filename{10} = '230319_014';
% filename{11} = '230319_015';

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
