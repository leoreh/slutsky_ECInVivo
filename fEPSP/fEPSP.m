
% this is a wrapper for fEPSP analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basepath = 'D:\Data\fEPSP\Ortal\Acute12';
graphics = true;
saveFig = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nsessions = 1;
amp = stability('nsessions', nsessions, 'basepath', basepath,...
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


for i = 1 : length(filename)
    data = import_wcp(fullfile(basepath, [filename{i} '.wcp']));
    
    % find artifact onset from derirative For each protocol the artifact is at
    % the same location, hence only calculated for the 1st trace
    [~, art] = max(diff(data.S(:, 1)));
    
    S{i} = rmDC(data.S, art, fs);
    
    pk = abs(min(S{i}(art + 0.004 * fs : end, :)));
    p(i) = mean(pk);
end

figure
plot(intensity, p, '*')

figure
for i = 1 : length(filename)
    plot(mean(S{i}, 2))
    hold on
end
legend



