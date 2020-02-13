function plotIsolation(basepath, spikes, saveFig)

% INPUT
%   basepath    path to recording
%   spikes      must include isolation fields (L, iDist, ISI). See cluValid
%   saveFig     save figure {true} or not (false)
%
% 04 dec 18 LH. 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % arguments
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nargs = nargin;
% if nargs < 1 || isempty(basepath)
%     basepath = pwd;
% end
% if nargs < 2 || isempty(spikes)
%     warning('spikes will be loaded from %s', basepath)
%     spikes = getSpikes('basepath', basepath);
% end
% nunits = length(spikes.UID);
% if nargs < 3 || isempty(clu)
%     clu = 1 : nunits;
% elseif length(clu) > nunits
%     error('specified more units to plot than are available in spikes')
% end
% if nargs < 4 || isempty(saveFig)
%     saveFig = true;
% end
% if ~isfield(spikes, 'L')
%     spikes.L = nan(nunits, 1);
% end
% if ~isfield(spikes, 'iDist')
%     spikes.iDist = nan(nunits, 1);
% end
% if ~isfield(spikes, 'isi')
%     spikes.isi = nan(nunits, 1);
% end
% if ~isfield(spikes, 'su')
%     spikes.su = nan(nunits, 1);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grpcolor = ['k', 'b', 'r', 'm', 'g', 'y'];
Lval = log10([0 0.05]);
iDistval = log10([20 300]);
ISIval = log10([0 1]);

f = figure;
subplot(2, 3, 1)
plot(log(spikes.isi), log10(spikes.lRat), '*')
ax = gca;
line(([ISIval(2) ISIval(2)]), ax.YLim)
line(ax.XLim, ([Lval(2) Lval(2)]))
xlabel('log(ISI ratio)')
ylabel('log(L ratio)')
axis tight
ax.TickLength = [0 0];
ax.Box = 0;
p = patch([ax.XLim(1) ISIval(2) ISIval(2) ax.XLim(1)], [ax.YLim(1) ax.YLim(1) Lval(2) Lval(2)], 'k');
p.FaceAlpha = 0.1;

subplot(2, 3, 2)
plot(log10(spikes.isi), log10(spikes.iDist), '*')
ax = gca;
line(([ISIval(2) ISIval(2)]), ax.YLim)
line(ax.XLim, ([iDistval(1) iDistval(1)]))
xlabel('log(ISI ratio)')
ylabel('log(Isolation Distance)')
axis tight
ax.TickLength = [0 0];
ax.Box = 0;
p = patch([ax.XLim(1) ISIval(2) ISIval(2) ax.XLim(1)], [iDistval(1) iDistval(1) ax.YLim(2) ax.YLim(2)], 'k');
p.FaceAlpha = 0.1;

subplot(2, 3, 3)
plot(log10(spikes.lRat), log10(spikes.iDist), '*')
ax = gca;
line(([Lval(2) Lval(2)]), ax.YLim)
line(ax.XLim, ([iDistval(1) iDistval(1)]))
p = patch([ax.XLim(1) Lval(2) Lval(2) ax.XLim(1)], [iDistval(1) iDistval(1) ax.YLim(2) ax.YLim(2)], 'k');
p.FaceAlpha = 0.1;
axis tight
ax.TickLength = [0 0];
ax.Box = 0;
xlabel('log(L ratio)')
ylabel('log(Isolation Distance)')

%%
subplot(2, 3, [4 : 6])
ax = gca;
scatter3(log10(spikes.lRat), log10(spikes.iDist), log10(spikes.isi), '*')
xlabel('log(L ratio)')
ylabel('log(Isolation Distance)')
zlabel('log(ISI ratio)')
hold on
ax = gca;
v = [Lval(2), iDistval(1), ax.ZLim(1); ax.XLim(1) iDistval(1) ax.ZLim(1); ax.XLim(1) ax.YLim(2) ax.ZLim(1); Lval(2) ax.YLim(2) ax.ZLim(1);...
    Lval(2), iDistval(1), ISIval(2); ax.XLim(1) iDistval(1) ISIval(2); ax.XLim(1) ax.YLim(2) ISIval(2); Lval(2) ax.YLim(2) ISIval(2)];
faces = [1 2 3 4; 1 4 8 5; 1 2 6 5; 2 3 7 6; 3 4 8 7; 5 6 7 8];
p = patch('Faces', faces, 'Vertices', v);
p.FaceAlpha = 0.1;
ax.TickLength = [0 0];
ax.Box = 0;

if saveFig
    savepdf('isolation distance', basepath, f)
end

%% bias of separation measurments to high-rate neurons
for i = 1 : length(spikes.UID)
    nspikes(i) = length(spikes.times{i});
end
nspikes = nspikes';

faces = figure;

subplot(1, 3, 1)
plot(log10(nspikes), log10(spikes.iDist), '*')
hold on
p = polyfit(log10(nspikes), log10(spikes.iDist), 1);
plot(log10(nspikes), polyval(p, log10(nspikes)), 'k', 'lineWidth', 2)
xlabel('log(nspikes)')
ylabel('log(Isolation Distanace)')
ax = gca;
ax.TickLength = [0 0];
ax.Box = 0;
title('Isolation Distanace')

subplot(1, 3, 2)
plot(log10(nspikes), log10(spikes.lRat), '*')
hold on
p = polyfit(log10(nspikes), log10(spikes.lRat), 1);
plot(log10(nspikes), polyval(p, log10(nspikes)), 'k', 'lineWidth', 2)
xlabel('log(nspikes)')
ylabel('log(L ratio)')
ax = gca;
ax.TickLength = [0 0];
ax.Box = 0;
title('L ratio')

subplot(1, 3, 3)
plot(log10(nspikes), spikes.isi, '*')
hold on
p = polyfit(log10(nspikes), spikes.isi, 1);
plot(log10(nspikes), polyval(p, log10(nspikes)), 'k', 'lineWidth', 2)
xlabel('log(nspikes)')
ylabel('ISI ratio')
ax = gca;
ax.TickLength = [0 0];
ax.Box = 0;
title('ISI ratio')
suptitle('bias of separation measurments to high-rate neurons')

if saveFig
    savepdf('separation bias', basepath, faces)
end


end

% EOF 
