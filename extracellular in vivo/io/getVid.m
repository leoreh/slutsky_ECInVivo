function vid = getVid(filename, varargin)

% gets video struct from .xls of ToxTrac (Rodriguez et al., 2018)
%  
% INPUT
%   filename    string. filename of .xls file
%   basepath    string. path to load filename and save output {pwd}
%   sheet       string. sheet to extract data from {'Seq_0_Arena_1_1'}
%   graphics    plot figure {1}.
%   saveFig     save figure {1}.
%   saveVar     save variable {1}.
% 
% OUTPUT
%   vid         structure with the following fields:
%       fps         frames per second of video file
%       nframes     number of frames analyzed
%       drate       detection rate of object (mouse)
%       t           frame timestamps
%       x           x coordinates (pixels)
%       y           y coordinates (pixels)
%       spd         instantaneous speed. starts at third frame
%       acc         instantaneous acceleration. starts at fifth frame
% 
% TO DO LIST
%     understand how spd and acc is calculated in ToxTrac
%     calibrate camera and arena in ToxTrac
%     convert pixels to mm
% 
% 21 mar 19 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'sheet', 'Seq_0_Arena_1_1', @ischar);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'saveFig', false, @islogical);
addOptional(p, 'saveVar', true, @islogical);

parse(p,varargin{:})
basepath = p.Results.basepath;
sheet = p.Results.sheet;
graphics = p.Results.graphics;
saveFig = p.Results.saveFig;
saveVar = p.Results.saveVar;

cd(basepath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get and arrange data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load xls
[~, txt, raw] = xlsread(filename, sheet);

vid.fps = xlsread(filename, sheet, 'B34');
vid.nframes = xlsread(filename, sheet, 'B35');
vid.drate = xlsread(filename, sheet, 'F15');
vid.res = xlsread(filename, sheet, 'B32:C32');

% position
idx = find(max(contains(txt, 'POSITIONS')));
vid.t = cell2mat(raw(6 : end, idx));
vid.x = cell2mat(raw(6 : end, idx + 2));
vid.y = vid.res(2) - cell2mat(raw(6 : end, idx + 3));

% speed
idx = find(max(contains(txt, 'INSTANT SPEED')));
vid.spd = cell2mat(raw(6 : end, idx + 2));
vid.spd(3 : end) = vid.spd(1 : end - 2);
vid.spd(1 : 2) = nan; 

% acceleration
idx = find(max(contains(txt, 'INSTANT ACCELERATION')));
vid.acc = cell2mat(raw(6 : end, idx + 2));
vid.acc(5 : end) = vid.acc(1 : end - 4);
vid.acc(1 : 4) = nan; 

% save variable
if saveVar   
    [~, filename] = fileparts(basepath);
    save([basepath, '\', filename, '.vid.mat'], 'vid')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    figure
    subplot(1, 2, 1)
    plot(vid.x, vid.y, 'Color', 'k', 'LineWidth', 1.5)
    xlim([0 vid.res(1)])
    ylim([0 vid.res(2)])
    axis tight
    xlabel('X Coordinates [px]')
    ylabel('Y Coordinates [px]')
    title('ToxTrac')
    if saveFig
        filename = 'VideoTrajectory';
        savePdf(filename, basepath, f)
    end
end

end