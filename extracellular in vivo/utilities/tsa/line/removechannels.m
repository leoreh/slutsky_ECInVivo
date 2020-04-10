% removechannels        from dat file
%
% [ rc, msg ] = removechannels( fname, nchans, chans )
% 
% fname         full file name 
% nchans        number of channels
% chans         channels to remove
%
% does
% rewrites the file without the chans
%
% see also      removedcoffset, reorderchannels

% 09-aug-13 ES

function [ rc, msg ] = removechannels( fname, nchans, chans )

% initialize output
rc = 0;
msg = '';

% constants
nbytes = 2;         % [bytes/sample]
blocksize = 1e6;    % [samples/channel]

% arguments
if nargin < 3 || isempty( fname ) || isempty( nchans ) || isempty( chans )
    return
end
if ~isa( fname, 'char' ) || ~exist( fname, 'file' ) 
    msg = sprintf( 'missing source %s', fname );
    rc = -1;
    return
end
tmpfname = [ fname '.tmp' ];
if ~isa( nchans, 'numeric' ) || ~isa( chans, 'numeric' ) ...
        || max( chans ) > nchans || min( chans ) < 1 ...
        || sum( chans ~= round( chans ) )
    msg = 'improper format for chans/nchans';
    rc = -1;
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

% open file for writing
fid = fopen( tmpfname, 'w' );
if fid == -1
    msg = 'cannot open file';
    rc = fid;
    return
end

% go over blocks and write out
for i = 1 : nblocks
    boff = ( blocks( i, 1 ) - 1 ) * nbytes * nchans;
    bsize = ( diff( blocks( i, : ) ) + 1 );
    m = memmapfile( fname, 'Format', 'int16', 'Offset', boff, 'Repeat', bsize * nchans, 'writable', true );
    d = reshape( m.data, [ nchans bsize ] );
    d( chans, : ) = [];
    fwrite( fid, d( : ), 'int16' );
    clear d m
end

% close file
rc = fclose( fid );
if rc == -1
    msg = 'cannot save new file';
    return
end

% remove original file and rename new one
if isunix
    cmd = sprintf( '!rm %s; mv %s %s', fname, tmpfname, fname );
elseif ispc
    cmd = sprintf( '!del %s; rename %s %s', fname, tmpfname, fname );
else
    msg = 'unrecognized system';
    rc = -1;
    return
end
eval( cmd )

return

% EOF

