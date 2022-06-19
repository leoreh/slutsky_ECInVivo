
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepaths = [
    {'F:\Data\lh93\lh93_210811_102035'},...         % acsf
    {'F:\Data\lh95\lh95_210824_083300'},...
    {'F:\Data\lh96\lh96_211201_070100'},...
    {'F:\Data\lh99\lh99_211218_090630'},...
    {'F:\Data\lh100\lh100_220405_100406'},...
    {'F:\Data\lh107\lh107_220509_095738'},...
    {'F:\Data\lh93\lh93_210813_110609'},...         % ket
    {'F:\Data\lh95\lh95_210825_080400'},...
    {'F:\Data\lh96\lh96_211204_084200'},...
    {'F:\Data\lh100\lh100_220403_100052'},...
    {'F:\Data\lh107\lh107_220501_102641'},...
    ];

nsessions = length(basepaths);
varsFile = ["fr"; "spikes"; "datInfo"; "session"; "cell_metrics.cellinfo";...
    "units"; "swv_metrics"; "st_metrics"];
varsName = ["fr"; "spikes"; "datInfo"; "session"; "cm";...
    "units"; "swv"; "st"];
v = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);

% cat units
units = catfields([v.units], 'catdef', 'long', 'force', false);
rs = find(units.clean(1, :));
fs = find(units.clean(2, :));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize waveforms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% resample waveforms. this is only because tdt and oe have
% different sampling rates
fs = 20000;
wv_len = 32;        % [samples]
x_resamp = linspace(0, 1, wv_len);
x_time = [1 : wv_len] / fs * 1000; 
x_time = x_time - x_time(ceil(wv_len / 2));

% initialize
wv = nan(length(units.clean), wv_len);
cnt = 1;

% go over sessions, resample and concat
for isession = 1 : nsessions
    wv_tmp = v(isession).swv.wv;
    if size(wv_tmp, 2) ~= wv_len
        x_orig = linspace(0, 1, size(wv_tmp, 2));
        wv_tmp = [interp1(x_orig, wv_tmp', x_resamp, 'spline', nan)]';
    end
    wv(cnt : cnt + size(wv_tmp, 1) - 1, :) = wv_tmp;
    cnt = cnt + size(wv_tmp, 1);
end

% normalize waveforms 
wv_norm = wv ./ vecnorm(wv, 2, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize classification metrics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

st = catfields([v.st], 'catdef', 'long', 'force', false);
swv = catfields([v.swv], 'catdef', 'long', 'force', false);
fr = catfields([v.fr], 'catdef', 'long', 'force', false);

xdata = swv.tp;
ydata = st.royer;
zdata = bz_NormToRange(fr.mfr, [10 50]);
ccbins = v(1).st.info.acg_narrow_tstamps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fig params
nrows = 2;
fh = figure;
th = tiledlayout(nrows, 2);
th.TileSpacing = 'compact';
th.Padding = 'none';

% plot rs and fs waveforms
axh(1) = nexttile(1);
hold on
plot(x_time, wv_norm(units.clean(2, :), :), 'b')
plot(x_time, wv_norm(units.clean(1, :), :), 'r')

% burst vs. tp separation, size is mfr
axh(2) = nexttile(2);
hold on
scatter(xdata(units.clean(2, :)), ydata(units.clean(2, :)), zdata(units.clean(2, :)), 'b', 'filled')
scatter(xdata(units.clean(1, :)), ydata(units.clean(1, :)), zdata(units.clean(1, :)), 'r', 'filled')
set(gca, 'yscale', 'log')
xlabel('Trough To Peak (ms)')
ylabel('Burst Index')

% example acg of rs
axh(3) = nexttile(3);
cc = squeeze(st.acg_narrow(:, :, rs(1)));
plot_ccg(cc, ccbins, 'clr', 'r',...
    'pred', [], 'sigbins1', [], 'sigbins2', [])

% example acg of fs
axh(4) = nexttile(4);
cc = squeeze(st.acg_narrow(:, :, fs(1)));
plot_ccg(cc, ccbins, 'clr', 'b',...
    'pred', [], 'sigbins1', [], 'sigbins2', [])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reaclculate stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for isession = 1 : 2
    basepath = basepaths{isession};
    cd(basepath)
    st = spktimes_metrics('bins', [], 'forceA', true);
end






cell_metrics = CellExplorer('basepaths', basepaths);





cell_metrics = CellExplorer('basepath', pwd);
v(isession).cm.putativeCellType
fh = gcf;
x = fh.Children(3).Children.XData;
y = fh.Children(4).Children.YData;
figure
plot(x, y)







