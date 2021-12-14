% ParseArgPairs     flexible argument assigning
%
% [ arg1name arg2name ... argNname ] = ParseArgPairs( ArgNames, DefaultValues, varargin )
% 
% where varargin is a list of parameter name/parameter value pairs, i.e. 
%   argName1, argVal1, ..., argNameM, argValM
% 
% for instance, the call
% [ shanknums, Fs, binsize ] = ParseArgPairs( { 'shanknums', 'Fs', 'binsize' }, { 1 : 6, 20000, 20 }, varargin{ : } );
% where varargin = { 'shanknums', 1 : 4, 'binsize', 10 }
%
% will result in assigning the following values:
% shanknums = 1 : 4;
% Fs = 20000;
% binsize = 10;
%
% note that: 
% -the number of elements in ArgNames, DefaultValues, and the number of
% output arguments must be the same
% -a call with mismatching ArgNames and argName will result in outName 
% being assigned the default value (but argName is NOT case sensitive)

% 04-mar-13 ES

function varargout = ParseArgPairs( ArgNames, DefArgs, varargin )

% check input arguments
if isempty( ArgNames )
    ArgNames = { [] };
end
if ~iscell( DefArgs )
    DefArgs = { DefArgs };
end
nDefArgs = length( DefArgs );
if length( ArgNames ) ~= nDefArgs
    error( '%s: mismatch: ArgNames and DefArgs must be the same length', upper( mfilename ) )
end
if nargout ~= nDefArgs
    error( '%s: ArgNames length and number of output arguments must be the same', upper( mfilename ));
end
for i = 1 : nDefArgs
    if ~ischar( ArgNames{ i } )
        error( '%s: ArgNames{ %d } must be a string', upper( mfilename ), i )
    else
        ArgNames{ i } = lower( ArgNames{ i } );
    end
end

% handle the argument pairs
params = cell( 1 );
if length( varargin ) == 1
    if length( varargin{ 1 } ) >= 2
        params = varargin{ : };
    end
elseif length( varargin ) >= 2
    params = varargin;
end
if ~isempty( params ) && length( params ) >= 2
    paramNames = params( 1 : 2 : end );
    paramValues = params( 2 : 2 : end );
    nParams = length( paramNames );
    paramValues = paramValues( 1 : nParams );
else
    nParams = 0;
    paramNames = {};
    paramValues = {};
end
for i = 1 : nParams
    if ~ischar( paramNames{ i } )
        error( '%s: paramNames{ %d } must be a string', upper( mfilename ), i )
    else
        paramNames{ i } = lower( paramNames{ i } );
    end
end
missing = paramNames( ~ismember( paramNames, ArgNames ) );
for i = 1 : length( missing )
    fprintf( '%s: parameter ''%s'' not specified in ArgNames; ignored!\n', upper( mfilename ), missing{ i } )
end

% assign values to output arguments
for i = 1 : nDefArgs
    argName = ArgNames{ i };
    idx = find( ismember( paramNames, argName ) );
    if ~isempty( idx ) && ~isempty( paramValues{ idx } )
        varargout( i ) = { paramValues{ idx } };
    else
        varargout( i ) = { DefArgs{ i } };
    end
end

return

% EOF
