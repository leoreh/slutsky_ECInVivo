function pk = io(varargin)

% loads data specified in filename. Each file should contain traces with
% the same stimulus intensity. allows user to remove unwanted traces.
% plots amplitude as a function of stimulus intensity and superimposed
% traces.
%
% INPUT
%   filename    cell array of strings describing filename number of stability sessions to concatenate
%   intensity   vector of intensities (mA). The order must correspond to
%               the order of filename or a GUI will be used.
%   inspect     logical. inspect traces {1} or not (0).
%   basepath    recording session path {pwd} to save figure and variables
%   graphics    logical. plot graphics {1} or not.
%   saveFig     logical. saveFig to current path {1} or not (0).
%   saveVar     logical. save output to current path {1} or not (0).
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
addOptional(p, 'basepath', pwd);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'saveFig', false, @islogical);
addOptional(p, 'saveVar', false, @islogical);

parse(p,varargin{:})
filename = p.Results.filename;
inspect = p.Results.inspect;
intensity = p.Results.intensity;
basepath = p.Results.basepath;
graphics = p.Results.graphics;
saveFig = p.Results.saveFig;
saveVar = p.Results.saveVar;

if length(intensity) ~= length(filename)
    if isempty(filename)
        for i = 1 : length(intensity)
            [filename{i}, ~, ~] = uigetfile('*.wcp');
        end
    else
        error('intensity and filename must be of equal length')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get and analyse data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : length(filename)
    data = import_wcp(fullfile(basepath, [filename{i} '.wcp']));
    
    % find artifact onset from derirative
    [~, art] = max(diff(data.S(:, 1)));
    % set minimum artifact if derirative fails
    art = max(art, 0.01 * data.fs);
    
    % remove DC
    data.S = rmDC(data.S, [1, art - 0.003 * data.fs]);
    
    % manually inspect and remove unwanted traces
    if inspect
        [data.S, data.rm_idx] = rmTraces(data.S);
    end
    
    pk(i) = mean(abs(min(data.S(art + 0.004 * data.fs : end, :))));
    tr(i, :) = mean(data.S, 2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if graphics
    
    f = figure;
    title('Input-Output')
    
    subplot(2, 1, 1)
    plot(intensity, pk, '-o', 'Color', 'k')
    xlabel('Intensity [mA]')
    ylabel('Amplitude [mV]')
    axis tight
    box off
    set(gca,'TickLength',[0, 0])
    
    subplot(2, 1, 2)
    plot(data.T, tr)
    xlabel('Time [s]')
    ylabel('Amplitude [mV]')
    axis tight
    box off
    set(gca,'TickLength',[0, 0])
    
    for i = 1 : length(filename)
        lg{i} = sprintf('%.2f mA', intensity(i));
    end
    legend([lg])
    
    if saveFig
        filename = 'Stability';
        savePdf(filename, basepath, f)
    end
end

if saveVar
    cd(basepath)
    save(pk_io, pk);
end

end
