function output_txt = curserLegend(~, event_obj)

% adds legend entry to data cursor mode text
% the legend must be explicitly defined
% 
% EXAMPLE: 
%           datacursormode on;
%           dcm = datacursormode(gcf);
%           set(dcm, 'UpdateFcn', @curserLegend)
% 
% 09 mar 19 LH

pos = get(event_obj, 'Position');
output_txt = {...
    ['X: ', num2str(pos(1), 4)]...
    ['Y: ', num2str(pos(2), 4)] ...
    ['Trace:', event_obj.Target.DisplayName]};

end