% this is a wrapper for fEPSP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% import data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = import_wcp();

% parameters
fs = 1 / diff(data.T(1:2));

% plot data and examine manually that all traces are OK
figure
plot(data.T, data.S)
axis tight
box off
xlabel('Time [s]')
ylabel('Voltage [mV]')
title('Traces')

data.S = rmDC(data.S, fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find peak amplitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% amplitude relative to baseline. index can be used to validate response
% position
p = min(data.S(art + 0.002 * data.fs : end, :));

% time axis in minutes
dt = 15;
x = [1 : dt : length(p) * dt] / 60;

figure
plot(x, p, '*')
axis tight
box off
xlabel('Time [m]')
ylabel('Amplitude [mV]')
title('Amplitude')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IO protocol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load traces. make sure order of filename is according to stimulus intensity
intensity = [0.2 : 0.1 : 0.8, 1];
basepath = 'E:\Data\fEPSP\Ortal\Acute12';
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
    
    % find artifact onset from derirative. For each protocol the artifact is at
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



