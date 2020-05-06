%%% user input
basepath = 'E:\Data\Others\DZ\IIS\WT';
mouse = 2;

%%% load data
cd(basepath)
filename = dir('*.lfp.*');
files = natsort({filename.name});
[~, basename] = fileparts(files{nfiles(mouse)});
[~, basename] = fileparts(basename);
load([basename '.lfp.mat'])
load([basename '.bs.mat'])

sig = lfp.data;

%%% manually marks burst-suppresion
ep = markEp(lfp.timestamps / 60, lfp.data(:, 1));

%%% convert episodes of burst to binary vector
x = floor(ep * fs) * 60;
epbin = zeros(length(sig), 1);
for i = 1 : size(x, 1)
    epbin(x(i, 1) : x(i, 2)) = 1;
end

%%% percent accuracy
accuracy = 1 - sum(abs(bs.binary - epbin)) / length(sig);

%%% graphics
figure
subplot(2, 1, 1)
plot(lfp.timestamps / 60, sig, 'k')
Y = ylim;
hold on
fill([bs.stamps fliplr(bs.stamps)]' / fs / 60, [Y(1) Y(1) Y(2) Y(2)],...
    'r', 'FaceAlpha', 0.25,  'EdgeAlpha', 0);
fill([ep fliplr(ep)]', [Y(1) Y(1) Y(2) Y(2)],...
    'g', 'FaceAlpha', 0.25,  'EdgeAlpha', 0);
xlabel('Time [m]')
ylabel('Voltage [mV]')

subplot(2, 1, 2)
plot(lfp.timestamps / 60, sig, 'k')
Y = ylim;
hold on
fill([bs.stamps fliplr(bs.stamps)]' / fs / 60, [Y(1) Y(1) Y(2) Y(2)],...
    'r', 'FaceAlpha', 0.25,  'EdgeAlpha', 0);
fill([ep fliplr(ep)]', [Y(1) Y(1) Y(2) Y(2)],...
    'g', 'FaceAlpha', 0.25,  'EdgeAlpha', 0);
xlim([12 15])
xlabel('Time [m]')
ylabel('Voltage [mV]')
legend({'Raw', 'Algorithm', 'Manual (green)'})