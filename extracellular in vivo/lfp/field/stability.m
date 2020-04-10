function amp = stability(varargin)

% loads data through GUI
% allows user to remove unwanted traces
% plots amplitude as a function of time
% repeats for specified number of sessions
% 
% INPUT
%   nsessions   number of stability sessions to concatenate
%   inspect     logical. inspect traces {1} or not (0).
%   basepath    recording session path {pwd} to save figure and variables
%   graphics    logical. plot graphics {1} or not.
%   saveFig     logical. saveFig to current path {1} or not (0).
%   saveVar     logical. save output to current path {1} or not (0).
%   pkMet       method to calculate peak amplitude
%               'absMin'    for each trace finds absolute minimum
%               'avgMin'    calculates peak of all traces at the same
%                           position, determined by the average minimum 
%                           excluding outliers
%               timestamp   calculates peak of all traces at the number
%                           specified (in sec)
% 
% OUTPUT
%   amp         vector of amplitude (min - baseline) for each trace  
% 
% 09 mar 19 LH  updates:
% 24 mar 19 LH  added pkMet

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'nsessions', 1);
addOptional(p, 'inspect', true, @islogical);
addOptional(p, 'basepath', pwd);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'saveFig', false, @islogical);
addOptional(p, 'saveVar', false, @islogical);
addOptional(p, 'pkMet', 'avgMin');

parse(p,varargin{:})
nsessions = p.Results.nsessions;
inspect = p.Results.inspect;
basepath = p.Results.basepath;
graphics = p.Results.graphics;
saveFig = p.Results.saveFig;
saveVar = p.Results.saveVar;
pkMet = p.Results.pkMet;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get and analyse data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1 : nsessions
    data = import_wcp();
    
    % find artifact onset from derirative
    [~, art] = max(diff(data.S(:, 1)));
    
    % remove DC
    data.S = rmDC(data.S, [1, art - 0.003 * data.fs]);
    
    % manually inspect and remove unwanted traces
    if inspect
        [data.S, data.rm_idx] = rmTraces(data.S);
    end
    
    % find amplitude
    if strcmp(pkMet, 'absMin')
        a{i} = abs(min(data.S(art + 0.003 * data.fs : end, :)));
    elseif strcmp(pkMet, 'avgMin')
        [~, idx] = min(data.S(art + 0.003 * data.fs : end, :));
        idx = round(trimmean(idx, 5));
        a{i} = abs(data.S(idx, :));
    elseif isnumeric(pkMet)
        idx = find(pkMet == data.T);
        a{i} = abs(data.S(idx, :));
    else
        error('pkMet unrecognized')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if graphics    
    amp = [a{:}];

    dt = 15;                                    % inter-trace interval [s]
    x = [1 : dt : length(amp) * dt] / 60;       % time axis in minutes

    f = figure;
    plot(x, amp, '*')
    axis tight
    box off
    xlabel('Time [m]')
    ylabel('Amplitude [mV]')
    title('Stability')
    set(gca,'TickLength',[0, 0])
    
    ntraces = 0;
    for i = 1 : length(a) - 1
        ntraces = ntraces + length(a{i});
        line([x(ntraces) x(ntraces)], get(gca, 'YLim'), 'Color', 'k', 'LineStyle', '--')
        % breakxaxis([x(ntraces), x(ntraces + 1)]); % works only for 2
        % sessions
    end
    
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
