function [mat, rm_idx] = rmTraces(mat, varargin)

% plots data and allows the user to remove unwanted traces
% 
% INPUT
%   mat         matrix of samples (rows) x traces (columns 
%   x           vector for x-axis, can be empty
%   basepath    recording session path {pwd} to save figure
%   saveFig     logical. saveFig to current path or not.
% 
% OUTPUT
%   mat         same as input but unwanted traces are nan
%   rm_idx      index of removed traces
% 
% TO DO LIST
%   replace while loop with interactive (linked) plot
% 
% 09 mar 19 LH  updates:
% 12 feb 20 LH  nan instead of removing traces



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

rm = nan;
rm_idx = [];
while ~isempty(rm)

    f = figure;
    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
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
    rm = input(prompt);
    mat(:, rm) = nan;
    
    if ~isempty(rm)
        close
    end
    
    rm_idx = [rm_idx; rm];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveFig
    filename = 'Traces';
    savePdf(filename, basepath, f)
end

end