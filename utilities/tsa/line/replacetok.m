% replacetok        replace token repetitively in string
%
% outstr = replacetok( instr, toknew, tokold )
%
% instr         any string array
% toknew        {'#'}; any string 
% tokold        {'.'}; any string
% 
% thus, 
% replacetok( instr );      replaces all '.' with '#'
% replacetok( instr, 'hello' ); replaces all '.' with 'hello'
% replacetok( instr, 'H', 'h' ); replaces all 'h' with 'H'
% replacetok( instr, '\ ', ' ' ); adds a backslash before each space

% 23-dec-12 ES

function outstr = replacetok( instr, toknew, tokold )

outstr = '';
nargs = nargin;
if nargs < 1 || isempty( instr ) || ~ischar( instr )
    return
end
instr = instr( : ).';
if nargs < 2 || isempty( toknew )
    toknew = '#';
end
if nargs < 3 || isempty( tokold )
    tokold = '.';
end

idx = strfind( instr, tokold ); 
idx = [ idx length( instr ) + 1 ];
outstr = instr( 1 : idx( 1 ) - 1 );
for i = 1 : ( length( idx ) - 1 )
    outstr = sprintf( '%s%s%s', outstr, toknew, instr( ( idx( i ) + 1 ) : ( idx( i + 1 ) - 1 ) ) ); 
end


return