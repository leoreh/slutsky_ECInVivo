% removeline            clean analogue signals in file from line influences
%
% call:                 [ rc msg xlta ] = removeline( filename, nchans, lineChanNum, varargin )
%
% gets:                 filename
%                       nchans              e.g. 41 if 41 channels
%                       lineChanNum         e.g. 41.6 (if last channel is composed of 16 digital channels, and the 6th is the line)
%
% optional arguments:   are given as argument name/value pairs:
%                       chansToClean        [] (i.e. all except the line channel)
%                       verbose             {0}
%                       graphics            {0}
%
% returns:              rc                  return code (-1: issues; 0: no problems)
%                       msg                 description of issues
%                       xlta                cell array with, for each cleaned channel, a matrix of the
%                                           interference (as a function of block size)
%
% calls:                for argument:   ParseArgPairs
%                       for algorithm:  lineRemove
%                       for graphics:   firfilt, replacetok
%
% called by:            massagedatfile
%
% see also:             removedcoffset, removechannels, getchannels

% 03-dec-17 ES

% revisions
% 04-dec-17 block design truly supported
% 05-dec-17 debugging, graphics
% 06-dec-17 (1) acquire line signals block-wise
%           (2) align blocks with line zero-crossing signals
% 07-dec-17 (1) DC computation added
%           (2) edges handled perfectly by weighted average in lineRemove

% to do (not critical):
%   presently, the line data are assumed to be in a single bit of an
%   int16 in the same file as the rest of the data. 
%   should extend to support:
%       -loading line data directly from a digital file
%       -loading line data as an int16 from an analog file
%   this does not involve any modifications to lineRemove or to the block
%   structure, just to the way lineChanNum is handled. also will require
%   additional arguments (digital file, additional analog file, etc)

function [ rc, msg, xlta ] = removeline( filename, nchans, lineChanNum, varargin )

% initialize output
rc              = 0;
msg             = '';
xlta            = [];

% constants
nbytes          = 2;            % [bytes/sample]
blocksize       = 1e6;          % [samples/channel]
lineFactor      = 61 / 20000;   % convert from samples to max possible number of events

% arguments
nargs = nargin;
if nargs < 3 || isempty( filename ) || isempty( nchans ) || isempty( lineChanNum )
    return
end
[ chansToClean ...
    , verbose, graphics  ] = ParseArgPairs(...
    { 'chansToClean' ...
    , 'verbose', 'graphics' } ...
    , {  [] ...
    , 1, 0 } ...
    , varargin{ : } );
if isempty( chansToClean )
    chansToClean = setdiff( 1 : nchans, floor( lineChanNum ) );
end
if ~isa( filename, 'char' ) || ~exist( filename, 'file' ) 
    msg         = sprintf( 'missing source %s', filename );
    rc          = -1;
    if verbose
        fprintf( '%s\n', msg )
    end
    return
end
if lineChanNum < 0 || lineChanNum > ( nchans + 1 )
    msg         = sprintf( 'inadequate digital line channel' );
    rc          = -1;
    if verbose
        fprintf( '%s\n', msg )
    end
    return
end

%------------------------------------------------------------------------
% preparations
%------------------------------------------------------------------------

% partition into blocks
info                = dir( filename );
nsamples            = info.bytes / nbytes / nchans;
if ~isequal( nsamples, round( nsamples ) )
    msg             = sprintf( 'incorrect nchans (%d) for file %s', nchans, fname );
    rc              = -1;
    return
end
nblocks             = ceil( nsamples / blocksize );
blocks              = [ 1 : blocksize : blocksize * nblocks; blocksize : blocksize : blocksize * nblocks ]';
blocks( nblocks, 2 ) = nsamples;

% get the line channel zero-crossing times from an analog file
if verbose
    t0          = clock;
    fprintf( 1, 'Getting line data...' )
end
linechannum         = floor( lineChanNum );
bitnum              = round( ( lineChanNum - linechannum ) * 10 );
lineT               = zeros( ceil( nsamples * lineFactor ), 1 );
j                   = 1;
for i = 1 : nblocks
    boff            = ( blocks( i, 1 ) - 1 ) * nbytes * nchans;
    bsize           = ( diff( blocks( i, : ) ) + 1 );
    m               = memmapfile( filename, 'Format', 'int16', 'Offset', boff, 'Repeat', bsize * nchans );
    d               = m.Data( linechannum : nchans : end );
    lineDataI       = rem( floor( double( d ) * pow2( 1 - bitnum ) ), 2 ) + 1;
    lineDiffI       = [ 0; diff( lineDataI ) ];
    lineTi          = find( lineDiffI > 0.5 );
    tidx            = j + [ 0 : length( lineTi ) - 1 ];
    lineT( tidx )   = lineTi  + blocks( i, 1 ) - 1 ;
    j               = j + length( lineTi );
end
clear m
tidx                = j : length( lineT );
lineT( tidx )       = [];
nmax                = max( diff( lineT ) );
if verbose
    et( 1 )         = etime( clock, t0 );
    fprintf( 1, 'done (%0.3g sec).\n', et( 1 ) )
end

% reorganize blocks to be aligned with line data
blocksNew           = blocks;
i = 1;
for j = 1 : length( lineT )
    if lineT( j ) > blocks( i, 2 )
        blocksNew( i, 2 )   = lineT( j ) - 1;
        i                   = i + 1;
        blocksNew( i, 1 )   = lineT( j );
    end
end
blocks              = blocksNew;

%------------------------------------------------------------------------
% go over blocks and compute DC
%------------------------------------------------------------------------
if verbose
    t0          = clock;
    fprintf( 1, 'Computing DC ' )
end
sumd = zeros( nchans, 1 );
for i = 1 : nblocks
    if verbose
        fprintf( '.' )
    end
    boff        = ( blocks( i, 1 ) - 1 ) * nbytes * nchans;
    bsize       = ( diff( blocks( i, : ) ) + 1 );
    m           = memmapfile( filename, 'Format', 'int16', 'Offset', boff, 'Repeat', bsize * nchans, 'writable', false );
    d           = reshape( m.data, [ nchans bsize ] );
    sumd        = sumd + double( sum( d, 2 ) );
    clear d m
end
dc              = sumd / nsamples;
if verbose
    etd         = etime( clock, t0 );
    fprintf( 1, 'done (%0.3g sec).\n', etd )
end

%------------------------------------------------------------------------
% go over blocks and remove line influences
%------------------------------------------------------------------------
xlta            = cell( 1, length( chansToClean ) );
kidx            = setdiff( 1 : nchans, chansToClean );
for i = 1 : nblocks
    if verbose
        t0      = clock;
        fprintf( 1, '\t\tBlock %d/%d ', i, nblocks )
    end
    % load data of all channels
    boff        = ( blocks( i, 1 ) - 1 ) * nbytes * nchans;
    bsize       = ( diff( blocks( i, : ) ) + 1 );
    m           = memmapfile( filename, 'Format', 'int16', 'Offset', boff, 'Repeat', bsize * nchans, 'writable', true );
    d           = reshape( m.data, [ nchans bsize ] );
    % remove line interference and update xlta
    for j = chansToClean
        if verbose
            fprintf( '.' )
        end
        if i == 1
            xin             = [];
        else
            xin             = xlta{ j }( i - 1, : );
        end
        linet               = lineT( lineT >= blocks( i, 1 ) & lineT <= blocks( i, 2 ) ) - blocks( i, 1 ) + 1;
        din                 = double( d( j, : )' );
        [ dout xout ]       = lineRemove( din, linet, xin, dc( j ), nmax );
        d( j, : )           = int16( dout );
        xlta{ j }( i, : )   = xout;
    end
    % write back to disk
    m.data      = d( : );
    clear d m
    if verbose
        et( i + 1 ) = etime( clock, t0 );
        fprintf( 1, 'done (%0.3g sec).\n', et( i + 1 ) )
    end
end

if verbose
    fprintf( 1, '\t\tAll done (%0.3g sec).\n', sum( et ) )
end

%------------------------------------------------------------------------
% graphics
%------------------------------------------------------------------------
if graphics
    
    Fs          = 20000;        % [Hz]
    nSec        = 10;           % MA [s]
    
    dt          = diff( lineT );
    m           = round( Fs / mean( dt ) * nSec);
    dtf         = firfilt( dt, ones( m, 1 ) / m );
    
    ncols       = 8;
    nrows       = ceil( length( chansToClean ) / ncols ) + 2;
    
    figure
    
    subplot( floor( nrows / 2 ), 1, 1 )
    plot( lineT( 2 : end ) / Fs / 60, dtf / Fs * 1000 )
    xlabel( 'Time [min]' )
    ylabel( 'Cycle duration [ms]' )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    lh = line( xlim, [ 20 20 ] );
    set( lh, 'linestyle', '--', 'color', [ 1 0 0 ] )
    lh = line( xlim, mean( dt ) / Fs * 1000 * [ 1 1 ] );
    set( lh, 'linestyle', '--', 'color', [ 0 0 1 ] )
    [ aa bb cc ]        = fileparts( filename );
    [ aa1 bb1 cc1 ]     = fileparts( aa );
    title( replacetok( bb1, '\_', '_' ) )
    
    for i = 1 : length( chansToClean ), 
        subplot( nrows, ncols, i + 2 * ncols )
        imagesc( xlta{ i } )
        axis off
        title( sprintf( '%d [%0.2g;%0.2g]', i, min( xlta{ i }( : ) ), max( xlta{ i }( : ) ) ) ); 
    end
    
end

return

% EOF


% call example:

filename            = '/Users/eranstark/Documents/slab/data/m589_171130_103429/all_in_one.dat';
nchans              = 41;
lineChanNum         = 41.6;
[ rc msg xlta ] = removeline( filename, nchans, lineChanNum, [], 1 );

