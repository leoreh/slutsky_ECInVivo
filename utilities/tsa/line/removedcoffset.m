% removedcoffset        remove DC offset
%
% removedcoffset( fname, nchans, fromfile )
%
% fname         full path to dat file
% nchans        number of channels
% fromfile      -1/{0}/1:
%                   -1  uses presaved values in *dco of same prefix (if
%                           available), does not rewrite
%                   0   computes anew from dat, saves *dco file
%                   1   uses existing *dco file (if available)
%
% does
% 1. computes each channel's mean (entire file)/loads the values in the *dco file
% 2. removes the mean from each channel
% 3. saves a *dco file (if not existing already)
%
% see also      removechannels, reorderchannels

% 09-aug-13 ES

function [ rc, msg ] = removedcoffset( fname, nchans, fromfile )

% initialize output
rc = 0;
msg = '';

% constants
nbytes = 2;         % [bytes/sample]
blocksize = 1e6;    % [samples/channel]

% arguments
if nargin < 2 || isempty( fname ) || isempty( nchans ) 
    return
end
if ~isa( fname, 'char' ) || ~exist( fname, 'file' ) 
    msg = sprintf( 'missing source %s', fname );
    rc = -1;
    return
end
if ~isa( nchans, 'numeric' ) 
    msg = 'improper format for nchans';
    rc = -1;
    return
end
if nargin < 3 || isempty( fromfile )
    fromfile = 0;
end
[ pathname filename ] = fileparts( fname );
dcofname = [ pathname '/' filename '.dco' ];

% partition into blocks
info = dir( fname );
nsamples = info.bytes / nbytes / nchans;
if ~isequal( nsamples, round( nsamples ) )
    msg = sprintf( 'incorrect nchans (%d) for file %s', nchans, fname );
    rc = -1;
    return
end
nblocks = ceil( nsamples / blocksize );
blocks = [ 1 : blocksize : blocksize * nblocks; blocksize : blocksize : blocksize * nblocks ]';
blocks( nblocks, 2 ) = nsamples;

% get the mean for each channel
rewrite = 1;
if fromfile ~= 0 && exist( dcofname, 'file' )
    meand = load( dcofname );
    if length( meand ) ~= nchans
        clear meand
        msg = sprintf( 'mismatch between dco file %s (%d values) and nchans (%d)'...
            , dcofname, length( meand ), nchans );
    end
    rewrite = 0;
elseif fromfile == -1
    rewrite = 0;
end

if ~exist( 'meand', 'var' )
    sumd = zeros( nchans, 1 );
    for i = 1 : nblocks
        boff = ( blocks( i, 1 ) - 1 ) * nbytes * nchans;
        bsize = ( diff( blocks( i, : ) ) + 1 );
        m = memmapfile( fname, 'Format', 'int16', 'Offset', boff, 'Repeat', bsize * nchans, 'writable', true );
        d = reshape( m.data, [ nchans bsize ] );
        sumd = sumd + double( sum( d, 2 ) );
        clear d m
    end
    meand = int16( sumd / nsamples );
end
meand = double( meand );

% go over blocks and remove mean
for i = 1 : nblocks
    boff = ( blocks( i, 1 ) - 1 ) * nbytes * nchans;
    bsize = ( diff( blocks( i, : ) ) + 1 );
    m = memmapfile( fname, 'Format', 'int16', 'Offset', boff, 'Repeat', bsize * nchans, 'writable', true );
    d = reshape( m.data, [ nchans bsize ] );
    d = d - int16( meand * ones( 1, bsize ) );
    m.data = d( : );
    clear d m
end

% save the dco file
if rewrite
    fid = fopen( dcofname, 'w' );
    if fid == -1
        msg = 'cannot open dco file';
        rc = fid;
        return
    end
    fprintf( fid, repmat( '%d\n', [ 1 nblocks ] ), meand( : ).' );
    rc = fclose( fid );
    if rc == -1
        msg = 'cannot save dco file';
    end
end

return

% EOF
