    
% params
basepath = 'C:\Users\LeoreHeim\Google Drive\PhD\Slutsky\Equipment\Running wheel\Results\30jun19';
filename = fullfile(basepath, 'Data.xlsx');
sheet = 'Data';
miceOrder = [5 1 2 3 4];
% miceOrder = [1 2 3 4];  % mouse 1 wheel not working

% load data
[num, ~, raw] = xlsread(filename, sheet);
data = num(8 : end, 1 : 5);
data = data(:, miceOrder);

% convert timestamps
t = datetime({raw{12 : 51, 1}});
t = [t, datetime({raw{52 : end, 1}})];
[h, m, ~] = hms(t);
[~, d, ~] = ymd(t);

t = strcat(num2str(h'), num2str(m'));
t(isspace(t)) = '0';
t = str2num(t);

idx = find(abs(h - 7) == 0);
idx = idx(logical([1 (diff(idx) > 3)]));
% idx2 = find(abs(h - 19) == 0);
% idx2 = idx2(logical([1 (diff(idx2) > 3)]));
% idx = sort([idx idx2])

% x = cumsum(m) / 60 / 24;

surgery = find(d == 11 & h == 11, 1);

% graphics
f = figure;
plot(movmean(data, 20))
hold on
y = ylim;
plot([surgery, surgery], y, '--k', 'LineWidth', 2)
xticks(idx)
xticklabels(repmat(['07:00'], length(idx), 1))
xlabel('Time [clock]')
ylabel('Running [turns]')
axis tight
box off
set(gca, 'TickLength', [0 0])
legend({'m1', 'm2', 'm3', 'm4', 'm5'})
