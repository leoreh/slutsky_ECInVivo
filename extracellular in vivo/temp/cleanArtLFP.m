% remove last minutes from recording

basepath{1} = 'E:\Data\Others\DZ\Field\Data\APPPS1_short';
basepath{2} = 'E:\Data\Others\DZ\Field\Data\APPPS1_long';
basepath{3} = 'E:\Data\Others\DZ\Field\Data\WT_long';
basepath{4} = 'E:\Data\Others\DZ\Field\Data\WT_short';

binsize = (2 ^ nextpow2(30 * 1250));
marg = round(fs * 0.05);

for i = 2 : length(basepath)
    cd(basepath{i})
    filename = dir('*.abf');
    files = natsort({filename.name});
    nfiles = 1 : length(files);
    
    for j = 1 : length(nfiles)
        [~, basename] = fileparts(files{nfiles(j)});
        
        lfp = getLFP('basepath', basepath{i}, 'ch', 1, 'chavg', {},...
            'fs', 1250, 'interval', [0 inf], 'extension', 'abf', 'pli', true,...
            'savevar', true, 'force', false, 'basename', basename);
        
        thr = [0 mean(lfp.data) + 5 * std(lfp.data)];
        iis.out = 0;
        iis = getIIS('sig', lfp.data, 'fs', 1250, 'basepath', basepath{i},...
            'graphics', false, 'saveVar', false, 'binsize', binsize,...
            'marg', 0.1, 'basename', basename, 'thr', thr, 'smf', 6,...
            'saveFig', false, 'forceA', true, 'spkw', true);
        
        art(i, j) = length(iis.out);
        idx = cell(1, length(iis.out));
        for k = 1 : length(iis.out)
            if iis.peakPos(iis.out(k)) + marg > length(lfp.data)
                idx{k} = iis.peakPos(iis.out(k)) - marg : length(lfp.data);
            else
                idx{k} = iis.peakPos(iis.out(k)) - marg : iis.peakPos(iis.out(k)) + marg;
            end
        end
        lfp.data([idx{:}]) = [];
        
        % remove first x minutes
        if i == 3
            if j == 4
                lfp.data(1 : fs * 60 * 21) = [];
            elseif j == 9
                lfp.data(1 : fs * 60 * 5) = [];
            end
        elseif i == 1
            if j == 4
                lfp.data(1 : fs * 60 * 1) = [];
            end
        end
        
        % inspect data
        figure; plot(lfp.timestamps, lfp.data)
        
        % remove max artifacts
        [~, x] = min(lfp.data);
        lfp.data(x - marg : x + marg) = [];

        % remove last 2 minutes
        lfp.data(end : -1 : end - 2 * 60 * lfp.fs) = [];

        % correct timestamps
        lfp.timestamps = [1 : length(lfp.data)] / fs;
        
        % save
        filename = [basename '.lfp.mat'];
        save([basepath{i}, filesep, filename], 'lfp')
    end
    
end



% basepath = 'E:\Data\Others\DZ\Field\Acute recordings\2h-3h\WT';
% cd(basepath)
% filename = dir('*.abf');
% files = {filename.name};
%
% i = 6;
%
% [~, basename] = fileparts(files{i});
%
% lfp = getLFP('basepath', basepath, 'chans', 1, 'chavg', {},...
%     'fs', 1250, 'interval', [0 inf], 'extension', 'abf', 'pli', true,...
%     'savevar', true, 'force', false, 'basename', basename);
%
% % remove last x minutes
% x = 2;
% lfp.data(end : -1 : end - x * 60 * lfp.fs) = [];
%
% % remove first x minutes
% x = 2;
% lfp.data(1 : fs * 60 * x) = [];
%
% % remove specific artifacts
% [~, x] = max(lfp.data);
% marg = fs * 0.2;
% lfp.data(x - marg : x + marg) = [];
%
% % correct timestamps
% lfp.timestamps = [1 : length(lfp.data)] / fs;
%
% % inspect data
% figure; plot(lfp.timestamps, lfp.data)
%
%
% filename = [basename '.lfp.mat'];
% save([basepath, filesep, filename], 'lfp')
%
%
%
%
%
