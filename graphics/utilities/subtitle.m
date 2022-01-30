function hh = subtitle(varargin)
%SUBTITLE Graph subtitle.
%   SUBTITLE('txt') adds the specified subtitle to the axes or chart 
%   returned by the gca command. Reissuing the subtitle command causes the 
%   new subtitle to replace the old subtitle.
%
%   SUBTITLE(...,'Property1',PropertyValue1,'Property2',PropertyValue2,...)
%   sets the values of the specified properties of the subtitle.
%
%   SUBTITLE(target,...) adds the subtitle to the specified target object.
%
%   H = SUBTITLE(...) returns the handle to the text object used as the 
%   subtitle.
%
%   See also XLABEL, YLABEL, ZLABEL, TEXT, TITLE.

%   Copyright 1984-2019 The MathWorks, Inc.

if nargout>0
    [isaxarr,hh]=matlab.graphics.chart.internal.objArrayDispatch(@subtitle,varargin{:});
else
    isaxarr=matlab.graphics.chart.internal.objArrayDispatch(@subtitle,varargin{:});
end
if isaxarr
    return
end

[ax,args,nargs] = labelcheck('Subtitle',varargin);

if nargs == 0  || (nargs > 1 && (rem(nargs-1,2) ~= 0))
  error(message('MATLAB:title:InvalidNumberOfInputs'))
end

if isempty(ax)
    ax = gca;
    % Chart subclass support
    % Invoke title method with same number of outputs to defer output arg
    % error handling to the method.
    if isa(ax,'matlab.graphics.chart.Chart')
        if(nargout == 1)
            hh = subtitle(ax,args{:});
        else
            subtitle(ax,args{:});
        end
        return
    end
end

titlestr = args{1};
if isempty(titlestr), titlestr=''; end
pvpairs = args(2:end);

% get-set does not support strings as of now
pvpairs = matlab.graphics.internal.convertStringToCharArgs(pvpairs);

matlab.graphics.internal.markFigure(ax);
h = get(ax,'Subtitle');
set(h, 'String', titlestr, pvpairs{:});

if nargout > 0
  hh = h;
end
