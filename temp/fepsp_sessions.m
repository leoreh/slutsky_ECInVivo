
% 
% % experiment folder
% exppath = 'G:\Data\hb_a2';
% cd(exppath)
% basepaths = dir('*_stp*');
% basepaths = fullfile({basepaths.folder}, {basepaths.name});
% nfiles = length(basepaths);
% 
% for ifile = 1 : nfiles
%     
%     ifile = 6;
%     cd(basepaths{ifile})
%     [~, basename] = fileparts(basepaths{ifile});
%     load([basename, '.lfp.mat'])
%     [~, wcpfiles] = fileparts(lfp.files)
%     intens = [50, 70];
% 
%     if ~iscell(wcpfiles)
%         wcpfiles = {wcpfiles};
%     end
% 
%     fepsp_wcpPipeline('basepath', exppath, 'wcpfiles', wcpfiles,...
%         'intens', intens, 'fsOut', [], 'recname', basename,...
%         'fepsp_protocol', lfp.fepsp_protocol)
% 
% end

% organizes recordings from wcp. assumes all .wcp files are in basepath and
% are named with unique numbers as suffixes (e.g. "xxx_034"). user selects
% specific files by specifying the last 2-3 digits of the filename. the
% function load the data, organizes it in a separate folder and calls the
% slutsky_fepsp functions.
% 
% INPUT
%   basepath        char. dir with .wcp files
%   slctfiles       numeric. files to select
%   fepsp_protocol  char. can be 'freerun', 'io', or 'stp'
%   recname         char. name of output folder
%   fsOut           numeric. requested sampling frequency. if empty will
%                   not downsample
%   intens          numeric. intensity values of stimulations [uA]
% 


% fepsp_sessions

% experiment folder
exppath = 'G:\Data\hb_a2';
cd(exppath)
basepaths = dir('*_io*');
basepaths = fullfile({basepaths.folder}, {basepaths.name});
nfiles = length(basepaths);

% stim session names
[~, basenames] = fileparts(basepaths);
ids = cellfun(@(x) x(end - 2 : end), basenames, 'uni', false);

% load data
varsFile = ["fepsp_traces"; "fepsp_results"; "lfp"];
varsName = ["traces"; "results"; "lfp"];
v = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);

% graphics
fh = figure;
hold on
for ifile = 1 : nfiles
    yval = cell2nanmat(v(ifile).results.avg_traces.Amp);
    xval = v(ifile).lfp.intens;
    plot(xval, yval)
end
legend(ids, 'Location', 'best')


