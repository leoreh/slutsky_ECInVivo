setMatlabGraphics(true)


% -------------------------------------------------------------------------
% mainly from https://www.mathworks.com/help/stats/mvncdf.html

% close all
% 
% % plot params
% n = 30;
% cLimit = [-0.1, 1.2];
% 
% % dist params
% m = [1 -1];
% s = [1, 0.2; 0.2, 0.6];
% 
% % create grid
% [x1, y1] = meshgrid(linspace(-1, 3, n)',linspace(-3, 1, n)');
% x = [x1(:) y1(:)];
% 
% % evaluate cdf
% p = mvncdf(x, m, s);
% z = reshape(p, n, n);
% 
% [x2, y2] = meshgrid(linspace(-4, -1, n)',linspace(-3, 3, n)');
% x = [x2(:) y2(:)];
% % x1 = x2; y1 = y2;
% % x = [y1(:) x1(:)];
% p2 = mvncdf(x, m, s);
% z2 = reshape(-p2, n, n);
% 
% % z = z + z2;
% % z = z2;
% z = bz_NormToRange(z, [0, 1]);
% 
% % plot
% 
% fh = figure;
% surf(x1, y1, z,...
%     'EdgeColor', [0.5, 0.5, 0.5], 'LineStyle', '-', 'FaceColor', 'flat')
% xlabel('Relevance')
% ylabel('Cost')
% zlabel('Interest')
% xticks([-1, 3])
% xticklabels([-1 1])
% yticks([-3, 1])
% yticklabels([1 0])
% zticks([0, 1])
% zticklabels([-1 1])
% colormap(flipud(bone))
% c = colorbar('TicksMode', 'manual', 'Ticks', cLimit, 'TickLabels', [0, 1]);
% c.Label.String = 'Worthwhile';
% c.Label.FontSize = 12;
% clim(cLimit)
% 
% cd('C:\Users\Leore\Downloads')
% set(fh,'Renderer','Painter')
% hgexport(fh,'2D.eps');
% 






% 
% 
% setMatlabGraphics(true)
% 
% 
% % -------------------------------------------------------------------------
% % mainly from https://www.mathworks.com/help/stats/mvncdf.html
% 
% % close all
% 
% % plot params
% n = 30;
% cLimit = [-0.1, 1.2];
% 
% % dist params
% m = [1.4 -2];
% s = [0.5, 0.2; 0.2, 0.6];
% 
% % create grid
% [x1, y1] = meshgrid(linspace(-1, 3, n)',linspace(-3, 1, n)');
% x = [x1(:) y1(:)];
% 
% % evaluate cdf
% p = mvncdf(x, m, s);
% 
% % plot
% z = reshape(p, n, n);
% z = bz_NormToRange(z, [0, 1]);
% 
% fh = figure;
% surf(x1, y1, z,...
%     'EdgeColor', [0.5, 0.5, 0.5], 'LineStyle', '-', 'FaceColor', 'flat')
% xlabel('Relevance')
% ylabel('Cost')
% zlabel('Interest')
% xticks([-1, 3])
% xticklabels([-1 1])
% yticks([-3, 1])
% yticklabels([1 0])
% zticks([0, 1])
% zticklabels([-1 1])
% colormap(flipud(bone))
% c = colorbar('TicksMode', 'manual', 'Ticks', cLimit, 'TickLabels', [0, 1]);
% c.Label.String = 'Worthwhile';
% c.Label.FontSize = 12;
% clim(cLimit)
% 
% cd('C:\Users\Leore\Downloads')
% set(fh,'Renderer','Painter')
% hgexport(fh,'2D.eps');


% -------------------------------------------------------------------------
scatter 

n = 10000;
s = 0.1;
m = 0;
distName = 'hn';
x = -random(distName, m, s, 1, n);       % relevant
y = random(distName, m, s, 1, n);        % time;
z = -random(distName, m, s, 1, n);       % interesting

c = mvncdf([x; y; z]', [1, 1, 1], scov);


scov = eye(3) .* [1, 1, 1] * 1000;
c = mvnpdf([x; y; z]', [1, 1, 1], scov);
% c = bz_NormToRange(c, [-1, 0]);

x = bz_NormToRange(x, [-1 1]);
y = bz_NormToRange(y, [0 1]);
z = bz_NormToRange(z, [-1 1]);


Limit = [0 20];
fh = figure;
scatter3(x, y, z, 8, c, 'filled')
colormap(gca, "bone")
colorbar('Direction', 'normal')
xlabel('Relevant')
ylabel('Time Consuming')
zlabel('Interesting')
% xlim(sort(-Limit))
% ylim(Limit)
% zlim(sort(-Limit))


fh = figure;
th = tiledlayout(1, 3);
axh = nexttile;
plot(histcounts(x))
axh = nexttile;
plot(histcounts(y))
axh = nexttile;
plot(histcounts(z))














% 
% -------------------------------------------------------------------------
% scatter 

n = 10000;
s = 0.1;
m = 0;
distName = 'hn';
x = -random(distName, m, s, 1, n);       % relevant
y = random(distName, m, s, 1, n);        % time;
z = -random(distName, m, s, 1, n);       % interesting

% c = mvncdf([x; y; z]', [1, 1, 1], scov);


scov = eye(3) .* [1, 1, 1] * 1000;
c = mvnpdf([x; y; z]', [1, 1, 1], scov);
% c = bz_NormToRange(c, [-1, 0]);

x = bz_NormToRange(x, [-1 1]);
y = bz_NormToRange(y, [0 1]);
z = bz_NormToRange(z, [-1 1]);


n = 25;
m = [1 -1, 1];
scov = [1, 0.6, 0.3; 0.6, 1, 0.3; 0.3, 0.3, 1];
% scov = eye(3) .* [1, 1, 1] * 100;

% create grid
xq = linspace(-1, 3, n)';
yq = linspace(-3, 1, n)';
zq = linspace(-3, 1, n)';
[x1, y1, z1] = meshgrid(xq, yq, zq);
x = [x1(:), y1(:), z1(:)];

% evaluate cdf
p = mvncdf(x, m, scov);
% z = reshape(p, n, n);




Limit = [0 20];
fh = figure;
scatter3(x1(:), y1(:), z1(:), 30, p, 'filled')
colormap(gca, "bone")
colorbar('Direction', 'normal')
xlabel('Relevant')
ylabel('Time Consuming')
zlabel('Interesting')
% xlim(sort(-Limit))
% ylim(Limit)
% zlim(sort(-Limit))


fh = figure;
th = tiledlayout(1, 3);
axh = nexttile;
plot(histcounts(x))
axh = nexttile;
plot(histcounts(y))
axh = nexttile;
plot(histcounts(z))













