% ED_CHTEST  Does ed_pickCh land on a good channel?
%
% Scored against ed_chSweep.mat, which ran the real detector on EVERY channel
% of raMCU3/4/5 and measured the AUC with which fastZ separates the curated
% discharges. Success = the pick's AUC is within ~0.01 of the best channel's,
% and its candidate count is not the pathological one.

DEVDIR = fileparts(mfilename('fullpath'));
load(fullfile(DEVDIR, 'ed_chSweep.mat'), 'res');

fprintf('%-8s %5s %7s %7s %7s %9s %9s\n', 'mouse', 'pick', 'AUC', ...
    'best', 'worst', 'nCand', 'bestCand');
for iB = 1 : numel(res)
    r = res{iB};
    basepath = fullfile('D:\Data\RA', r.basename(1 : 6), r.basename);

    ch = ed_pickCh(basepath, 'flgForce', true);
    iPick = find(r.ch == ch, 1);
    [aucBest, iBest] = max(r.auc);

    fprintf('%-8s %5d %7.3f %7.3f %7.3f %9d %9d\n', r.basename(1 : 6), ch, ...
        r.auc(iPick), aucBest, min(r.auc), r.nCand(iPick), r.nCand(iBest));
end
