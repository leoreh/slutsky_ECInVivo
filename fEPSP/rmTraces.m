function [mat, rm_idx] = rmTraces(mat, varargin)

% plots data and allows the user to remove unwanted traces
% 
% INPUT
%   mat         matrix of traces (columns) x samples (rows)
%   x           vector for x-axis, can be empty
%   basepath    recording session path {pwd} to save figure
%   saveFig     logical. saveFig to current path or not.
% 
% OUTPUT
%   mat         same as input without unwanted traces   
%   rm_idx      index of removed channels
% 
% TO DO LIST
%   replace while loop with interactive (linked) plot
% 
% 09 mar 19 LH


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'x', []);
addOptional(p, 'basepath', pwd);
addOptional(p, 'saveFig', false, @islogical);

parse(p,varargin{:})
x = p.Results.x;
basepath = p.Results.basepath;
saveFig = p.Results.saveFig;

if isempty(x)
    x = 1 : length(mat);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm_idx = nan;
while ~isempty(rm_idx)

    f = figure;
    plot(x, mat)
    axis tight
    box off
    xlabel('Time [s]')
    ylabel('Voltage [mV]')
    title('Traces')
    legend(string([1 : size(mat, 2)]))
    legend('off')
    set(gca,'TickLength',[0, 0])
    
    datacursormode on;
    dcm = datacursormode(gcf);
    set(dcm, 'UpdateFcn', @curserLegend)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get traces to remove
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    prompt = sprintf('\nWhich traces would you like to remove? type as vector, then press return.\n\n');
    rm_idx = input(prompt);
    mat(:, rm_idx) = [];
    
    if ~isempty(rm_idx)
        close
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveFig
    filename = 'Traces';
    savePdf(filename, basepath, f)
end

end