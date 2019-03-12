function pk = stp(varargin)

% loads data specified in filename. Each file should contain traces with
% the same stimulus intensity. allows user to remove unwanted traces.
% plots amplitude as a function of stimulus intensity and superimposed
% traces.
%
% INPUT
%   filename    cell array of strings describing filename number of stability sessions to concatenate
%   intensity   vector of intensities (mA). The order must correspond to
%               the order of filename
%   inspect     logical. inspect traces {1} or not (0).
%   nstim       number of stimulations in single trace.
%   basepath    recording session path {pwd} to save figure and variables
%   graphics    logical. plot graphics {1} or not.
%   saveFig     logical. saveFig to current path {1} or not (0).
%
% OUTPUT
%   pk          vector with peak amplitude for each intensity
%
% 09 mar 19 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'filename', {}, @iscell);
addOptional(p, 'intensity', [0.1 : 0.1 : 1]);
addOptional(p, 'inspect', true, @islogical);
addOptional(p, 'nstim', 5);
addOptional(p, 'basepath', pwd);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'saveFig', false, @islogical);

parse(p,varargin{:})
filename = p.Results.filename;
inspect = p.Results.inspect;
intensity = p.Results.intensity;
nstim = p.Results.nstim;
basepath = p.Results.basepath;
graphics = p.Results.graphics;
saveFig = p.Results.saveFig;

if length(intensity) ~= length(filename)
    if isempty(filename)
        for i = 1 : length(intensity)
            [filename{i}, ~, ~] = uigetfile('*.wcp');
        end
    else
        error('intensity and filename must be of equal length')
    end
end

% initialize vars
nfiles = length(filename);
d = zeros(1, nstim);
pk = zeros(nfiles, nstim);
lg = cell(1, nfiles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get and analyse data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : nfiles
    data = import_wcp(fullfile(basepath, [filename{i} '.wcp']));
    art_t = 0.003 * data.fs;

    % find artifact onset\s from derirative
    [~, idx] = sort(diff(data.S(:, 1)));
    didx = nan;
    while ~isempty(didx)
        art = idx(end : -1 : end - nstim + 1);
        % make sure maximums are from different stimuli
        for j = 1 : length(art) - 1
            d(j) = min(abs(art(j) - art(j + 1 : end)));
        end
        didx = find(d < 0.001 * data.fs & d > 0);
        for j = 1 : length(didx)
            idx(idx == art(didx(j))) = [];
        end
    end
    art = sort(art);
    
    % remove DC
    data.S = rmDC(data.S, [1, art(1) - art_t]);
    
    % manually inspect and remove unwanted traces
    if inspect
        [data.S, data.rm_idx] = rmTraces(data.S);
    end
    
    bl = mean(abs(min(data.S(art(1) + art_t : art(2) - art_t, :))));
    for j = 1 : nstim
        if j == nstim
            pk(i, j) = mean(abs(min(data.S(art(j) + art_t : end, :)))) / bl;
        else
            pk(i, j) = mean(abs(min(data.S(art(j) + art_t : art(j + 1) - art_t, :)))) / bl;
        end
    end
    
    tr(i, :) = mean(data.S, 2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if graphics
    
    f = figure;   
    subplot(2, 1, 1)
    for i = 1 :  nfiles
        plot([1 : nstim], pk(i, :), '-o')
        hold on
    end
    xlabel('Stimulus [#]')
    ylabel('Amplitude Change [%]')
    title('STP')   
    axis tight
    box off
    set(gca,'TickLength',[0, 0])   
    
    for i = 1 : nfiles
        lg{i} = sprintf('%.2f mA', intensity(i));
    end

    legend([lg])
    subplot(2, 1, 2)
    plot(data.T, tr)
    xlabel('Time [s]')
    ylabel('Amplitude [mV]')
    axis tight
    box off
    set(gca,'TickLength',[0, 0])
    
    if saveFig
        filename = 'Stability';
        savePdf(filename, basepath, f)
    end
end

end
