
% this is a wrapper for fEPSP analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basepath = 'C:\Users\LeoreHeim\Downloads\Data-20200317T100406Z-001\Data\hDlx-Gq\after CNO';
files = [2  : 5];

cd(basepath)
filename = dir('*.wcp');
filename = natsort({filename.name});
ch = [1];

% typical params for stp
start = [0.01 : 0.02 : 0.09];
nstim = 5;
% typical params for IO / stability
start = 0.03;
nstim = 1;

% go over files, calc response and save data as .mat
for i = files
    raw = import_wcp(fullfile(basepath, filename{i}));
    fprintf('\nloading %s\n', fullfile(basepath, filename{i}));
    [~, basename] = fileparts(filename{i});
    for j = ch
        [data.amp{j}, data.rm{j}] = getFieldAmp('sig', raw.S{j},...
            'fs', raw.fs, 'start', start, 'stop', [], 'inspect', false,...
            'basepath', basepath, 'nstim', nstim,...
            'graphics', true, 'saveVar', false, 'saveFig', true,...
            'filename', basename);
        data.sig{j} = raw.S{j};
    end   
    % arrange and save
    data.t = raw.T;
    data.fs = raw.fs;
    save([basename '.mat'], 'data')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% concatenate traces from different files
amp = cell(ch);
k = 1;
for i = files
    [~, basename] = fileparts(filename{i});
    load(basename)
    for j = ch
        amp{j} = [amp{j} data.amp{j}];
        trace{k, j} = data.sig{j};
    end
    ntraces(k) = size(amp{j}, 2);
    k = k + 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% vector of intensities (mA). The order must correspond to the order of
% files.
intensity = [0.2 : 0.1 : 0.4];
files = [2 : 4];

% arrange data
k = 1;
for i = files
    [~, basename] = fileparts(filename{i});
    load(basename)
    for j = ch
        amp{k, j} = (data.amp{j});
        trace{j}(k, :) = mean(data.sig{j}');
    end
    k = k + 1;
end
