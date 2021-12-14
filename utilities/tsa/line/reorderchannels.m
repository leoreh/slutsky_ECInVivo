% reorderchannels       reorder channnels
%
% reorderChannels( fname, neworder )
%
% fname         full path to dat file
% neworder      vector (nchans elements): neworder (1-based)
%
% does
% reorders the channels in the file
%
% see also      removedcoffset, removechannels, getchannels
% 
% use example:
% assume recording from Buzsaki32sp, plus 15 external channels
%
% neworder = [ 27 25 29 21 18 19 20 17 32 30 24 28 31 26 23 22 7 5 11 3 13 16 15 14 4 2 6 10 8 1 12 9 33 : 47 ];
% reorderchannels( 'orig.dat', neworder );
%
% see also example at removechannels.m, EOF 

% 09-aug-13 ES

% revisions
% 07-dec-17 do not reorder if neworder is trivial

function [ rc, msg ] = reorderchannels( fname, neworder )

% initialize output
rc = 0;
msg = '';

% constants
nbytes = 2;         % [bytes/sample]
blocksize = 1e6;    % [samples/channel]

% arguments
if nargin < 2 || isempty( fname ) || isempty( neworder )
    return
end
if ~isa( fname, 'char' ) || ~exist( fname, 'file' ) 
    msg = sprintf( 'missing source %s', fname );
    rc = -1;
    return
end
if ~isa( neworder, 'numeric' ) 
    msg = 'improper format for neworder';
    rc = -1;
    return
end
neworder = neworder( : );
nchans = length( neworder );
if isequal( [ 1 : nchans ]', neworder )
    return
end

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

% go over blocks and reorder
for i = 1 : nblocks
    boff = ( blocks( i, 1 ) - 1 ) * nbytes * nchans;
    bsize = ( diff( blocks( i, : ) ) + 1 );
    m = memmapfile( fname, 'Format', 'int16', 'Offset', boff, 'Repeat', bsize * nchans, 'writable', true );
    d = reshape( m.data, [ nchans bsize ] );
    d = d( neworder, : );
    m.data = d( : );
    clear d m
end

return

% EOF
