% mergefiles        two files that have exactly the same duration
%
% call              [ RC MSG ] = MERGEFILES( SOURCEFILE1, SOURCEFILE2, NEWFILE, NCHANS, PRECISION, CASTMODE )
%
% merge the data in SOURCEFILE1 and SOURCEFILE2
% and save it in a NEWFILE
% 
% the two channels must have the same duration, but may have a different
% number of channels: NCHANS( 1 ) and NCHANS( 2 )
%
% e.g. if file1, sample1 is:
%       c(1)t(1) c(2)t(1) .. c(n)t(1)
% and file2, sample1 is:
%       c(n+1)t(1) c(n+2)t(1) .. c(n+m)t(1)
% then fileNew, sample1 will simply be:
%       c(1)t(1) c(2)t(1) .. c(n)t(1) c(n+1)t(1) c(n+2)t(1) .. c(n+m)t(1)
% 
% PRECISION is an input to fread/fwrite
% it is a cell array of 2 strings, default is { 'int16', 'uint16' }
% the NEWFILE will be of PRECISION{ 1 }
%
% CASTMODE is how to do type casting if the two files consist of different
% data types; relevant only for uint16->int16 case. the options are -1/0/{1}:
% -1:   the C/MATLAB default for uint16->int16: clip numbers above 2^15-1 (loses data!!)
% 0     the "information preserving" option: translate (shift) to the -2^15 to 2^15-1 range (changes values!!)
% 1:    the "intuitive" option: compress to the 0-2^15-1 range (lose resolution!!)
%
% see also          CONCATFILES, EXTRACTFILE, PARTITION (work on the time dimension)
%                   SPLITFILE (works on the channel dimension)

% 17-mar-16 ES

% revisions
% 05-may-16 actually written the routine (previously just conceptualized)
% 03-dec-17 (1) initialized output; added msg
%           (2) modified handling of blocksize2

function [ rc msg ] = mergefiles( sourcefile1, sourcefile2, newfile, nchans, precision, castmode )

% initialize output
rc      = 0;
msg     = '';

% constants
nfiles      = 2;
BLOCKSIZE   = 2^20; % number of elements/block (not bytes)

% input arguments
nargs = nargin;
if nargs < 3 || isempty( sourcefile1 ) || isempty( sourcefile2 ) || isempty( newfile )
    error( 'missing input parameters' )
end
if nargs < 4 || isempty( nchans ) || length( nchans ) ~= 2
    nchans = [ 32 4 ];
end
if any( nchans <= 0 ) || any( nchans ~= round( nchans ) )
    error( 'nchans should be non-negative integers' )
end
if nargs < 5 || isempty( precision ) || length( nchans ) ~= 2
    precision = { 'int16', 'uint16' };
end
if ~strcmp( precision{ 1 }, precision{ 2 } ) 
    if ~strcmp( precision{ 1 }, 'int16' ) && ~strcmp( precision{ 2 }, 'uint16' )
        error( 'unsupported' )
    end
end
if nargs < 6 || isempty( castmode )
    castmode = 1;
end
if ~ismember( castmode, [ -1 0 1 ] )
    error( 'unsupported' )
end

% build the type casting strings
for i = 1 : nfiles
    precisionstr{ i } = sprintf( '*%s', precision{ i } );
end

% determine number of bytes/sample/channel
nbytes = zeros( nfiles, 1 );
for i = 1 : nfiles
    a = ones( 1, 1, precision{ i } );
    sourceinfo = whos( 'a' );
    nbytes( i ) = sourceinfo.bytes;
end

% get the number of time samples/channel:
a1              = dir( sourcefile1 );
nsamples1       = a1.bytes / nchans( 1 );
a2              = dir( sourcefile2 );
nsamples2       = a2.bytes / nchans( 2 );
datasize1       = [ nchans( 1 ) nsamples1 ];
datasize2       = [ nchans( 2 ) nsamples1 ];
datasize        = [ sum( nchans ) nsamples1 ] ;
if nsamples1 ~= nsamples2
    msg = 'file duration/nchans mismatch';
    rc = -1;
    return
end

% divide into blocks
nelements1      = prod( datasize1 );
nelements2      = prod( datasize2 );
nblocks1        = ceil( nelements1 / BLOCKSIZE );
nblocks2        = ceil( nelements2 / BLOCKSIZE );
nblocks         = max( [ nblocks1 nblocks2 ] );
% % older version:
% blocksize1      = 2.^ceil( log2( nelements1 / nblocks ) );
% blocksize2      = 2.^ceil( log2( nelements2 / nblocks ) );
blocksize1      = 2.^ceil( log2( nelements1 / nblocks ) );
blocksize1      = floor( blocksize1 / nchans( 1 ) ) * nchans( 1 );
blocksize2      = ceil( blocksize1 / nchans( 1 ) * nchans( 2 ) );
if ceil( blocksize1 / blocksize2 ) ~= nchans( 1 ) / nchans( 2 ) ...
        || blocksize1 ~= round( blocksize1 ) ...
        || blocksize2 ~= round( blocksize2 )
    error( 'debugging error' )
end

% open files for reading and writing
fp0             = fopen( newfile, 'w' );
fp1             = fopen( sourcefile1, 'r' );
fp2             = fopen( sourcefile2, 'r' );
if fp0 == -1, error( 'fopen error' ), end
if fp1 == -1, error( 'fopen error' ), end
if fp2 == -1, error( 'fopen error' ), end

% go over the sourcefile in blocks and write the newfile
for bnum = 1 : nblocks

    % collect the data
    if bnum == nblocks
        toload1 = nelements1 - ( nblocks - 1 ) * blocksize1;
        toload2 = nelements2 - ( nblocks - 1 ) * blocksize2;
    else
        toload1 = blocksize1;
        toload2 = blocksize2;
    end
    data1 = fread( fp1, toload1, precisionstr{ 1 } );
    data2 = fread( fp2, toload2, precisionstr{ 2 } );
    
    % convert uint16 to int16 (specifically for uint16 in sourcefile2)
    if strcmp( precision{ 1 }, 'int16' ) && strcmp( precision{ 2 }, 'uint16' )
        switch castmode
            case -1 % 'clip': 0:2^15-1 remains; 2^15:2^16-1 are clipped to 2^16-1
                data2 = int16( data2 );
            case 0 % 'shift': 0:2^16-1 is mapped 1:1 to -2^15:2^16-1
                data2 = int16( double( data2 ) - 2^15 );
            case 1 %'compress': 0:2^16-1 are mapped to 0:2^15-1
                data2 = int16( floor( double( data2 ) / 2 ) );
        end
    end
    
    % merge
    data1hat = reshape( data1, nchans( 1 ), size( data1, 1 ) / nchans( 1 ) );
    data2hat = reshape( data2, nchans( 2 ), size( data2, 1 ) / nchans( 2 ) );
    data0hat = [ data1hat; data2hat ];
    data0 = data0hat( : );
  
    % write out
    fwrite( fp0, data0, precision{ 1 } );
    
end

% close files
rc0 = fclose( fp0 );
rc1 = fclose( fp1 );
rc2 = fclose( fp2 );
rc = [ rc0 rc1 rc2 ];

return

% EOF


% testing:
fname1 = '/Volumes/Data3/phaser10/dan/uLED_1-160316/uLED_1-160316-01.dat';
fname2 = '/Volumes/Data3/phaser10/dan/uLED_1-160316/uLED_1-160316-01_analogin.dat';
newfile = '/Volumes/Data3/phaser10/dan/uLED_1-160316/uLED_1-160316-01_merged.dat';
nchans = [ 32 1 ];
precision = { 'int16', 'uint16' };
rc = mergefiles( fname1, fname2, newfile, nchans, precision );
