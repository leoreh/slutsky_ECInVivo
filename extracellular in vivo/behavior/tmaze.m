

% basepath = 'C:\Users\LeoreHeim\AppData\Local\Temp\OneNote';
% cd(basepath);
% filelist = dir('*.xlsx')

%%% load data
filename{1} = 'C:\Users\LeoreHeim\AppData\Local\Temp\OneNote\16.0\Exported\{5CA3CD89-AB3C-433F-BC60-B91343326B0C}\NT\9\May 22 - Spreadsheet10.xlsx';
filename{2} = 'C:\Users\LeoreHeim\AppData\Local\Temp\OneNote\16.0\Exported\{5CA3CD89-AB3C-433F-BC60-B91343326B0C}\NT\10\May 23 - Spreadsheet11.xlsx';
filename{3} = 'C:\Users\LeoreHeim\AppData\Local\Temp\OneNote\16.0\Exported\{5CA3CD89-AB3C-433F-BC60-B91343326B0C}\NT\11\May 24 - Spreadsheet12.xlsx';
filename{4} = 'C:\Users\LeoreHeim\AppData\Local\Temp\OneNote\16.0\Exported\{5CA3CD89-AB3C-433F-BC60-B91343326B0C}\NT\12\May 25 - Spreadsheet13.xlsx';
weightname = 'C:\Users\LeoreHeim\AppData\Local\Temp\OneNote\16.0\Exported\{5CA3CD89-AB3C-433F-BC60-B91343326B0C}\NT\8\Cohort 1  - Spreadsheet9.xlsx';
sheet = 'Sheet1';

ndays = length(filename);

durMat = [];
trialVec = [];
trialsuccess = [];
for i = 1 : ndays
    [num, txt, raw] = xlsread(filename{i}, sheet);
    session{i}.name = raw{1, 1};
    session{i}.weight =  num(2, 2 : 2 : end);
    session{i}.trialdur =  num(4 : end, 2 : 2 : end);
    session{i}.trialdir = num(5 : end, 3 : 2 : end);
    trialsuccess = [trialsuccess; session{i}.trialdir(:, 5 : 9)];
    session{i}.success = nanmean(session{i}.trialdir);
    successMat(i, :) = session{i}.success(~isnan(session{i}.success));
    durMat = [durMat; session{i}.trialdur];
    trialVec = [trialVec, 1 : length(session{i}.trialdur)];
end

%%% graphics
% success rate
figure
for i = 1 : 5
    plot([1 : ndays], successMat(:, i))
    hold on
end
plot([1 : ndays], mean(successMat, 2), 'LineWidth', 3, 'Color', 'k')
ylim([0.5 1])
xticks(1 : ndays)
xticklabels(num2cell([1 : ndays]));
yticks([0.5 0.75 1])
yticklabels(num2cell([0.5 0.75 1]));
box off
set(gca, 'TickLength', [0 0])
ylabel('Success Rate')
xlabel('Days')

% trialDur
figure
for i = 1 : 9
    plot(durMat(:, i), '*')
    hold on
end
plot(nanmean(durMat, 2), 'LineWidth', 3, 'Color', 'k')
axis tight
ylabel('Trial Duration [s]')
xlabel('Trial')
xticks(1 : length(trialVec))
xticklabels(num2cell(trialVec))
box off
set(gca, 'TickLength', [0 0])

% weight
[num, txt, raw] = xlsread(weightname, sheet);
weightMat = num([2, 5 : end], 2 : end);
figure
subplot(2, 1, 1)
plot(weightMat)
box off
set(gca, 'TickLength', [0 0])
ylabel('Weight')
xlabel('Days')
    
subplot(2, 1, 2)
plot(weightMat ./ weightMat(1, :))
box off
set(gca, 'TickLength', [0 0])
ylabel('Weight')
xlabel('Days')

% trial success
figure
subplot(2, 1, 1)
plot(trialsuccess)
xticks(1 : length(trialVec))
xticklabels(num2cell(trialVec))
box off
set(gca, 'TickLength', [0 0])
ylabel('Wrong / Right')
xlabel('Trial')
subplot(2, 1, 2)
plot(cumsum(trialsuccess, 'omitnan'))
xticks(1 : length(trialVec))
xticklabels(num2cell(trialVec))
box off
set(gca, 'TickLength', [0 0])
ylabel('Cummulative Success')
xlabel('Trial')

