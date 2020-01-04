% massagedatfile        remove channels, reorder channels, remove DC, clean line
%
% call:             rc = massagedatfile( filename, varargin )
%
% does:             takes a *dat file and 
%                   (1) {is  done by default} backs it up
%                   (2) {is  done by default} checks general integrity and removes samples which are non-integer multiples of nchans 
%                   (3) {not done by default} removes unnecessary channels
%                   (4) {not done by default} reorders remaining channels (e.g. to match adjacency matrix) 
%                   (5) {not done by default} removes DC from each channel
%                   (6) {not done by default} removes line-locked artifacts from each channel
%                   (7) {is  done by default} adds extra channels from other files of the same duration (e.g. digital channels)
%
% arguments:        given as argument name/value pairs, e.g.
%                   [ rc msg ] = massagedatfile( filename, 'newOrder', newOrder );
%
% gets:             filename            full path, name, and suffix
%
%                   makeBackup          {0} steps over the original file
%                                        1  before anything, makes a backup called filename.orig.suffix
%
%                   checkIntegrity      -1  check whether binary file consists of interger multiple of number of channels
%                                        0  do not check at all
%                                       {1} check and delete extra samples
%
%                   nChannelsOrig       can be left empty if newOrder is given
%
%                   chansToRemove       {[]}; 1-based, in original dat file order, can be left empty.
%
%                   newOrder            {[]}; 1-based, in original dat file order. can be left empty.
%
%                   dcoFile              0   computes anew from dat, saves *dco file
%                                        1   if available, uses existing *dco file; saves anew if not
%                                       -1   same as 1, but does not rewrite
%                                       {-2} does not remove DC
%
%                   filenameSecond      {'digitalin'}; will be merged head-to-tail
%                   nchansSecond        {1}; number of channels in 2nd file
%                   precisionSecond     {'uint16'}; precision of 2nd file
%
%                   lineChannel         {[]}; 1-based, in original dat file order or in lineFilename
%                                           if an analog channel in dat file (int16) - a 1-based integer
%                                                   NOT SUPPORTED IN THIS VERSION
%                                           if a digital channel in digitalin file (uint16) - a 1-based integer
%                                                   NOT SUPPORTED IN THIS VERSION
%                                           if a digital channel in a dat file (int16) - a decimal number
%                                               (e.g. 41.6 is the 6th bit in the 41st channel)
%                                                   FULLY SUPPORTED IN THIS VERSION
%                   lineFilename        {[]}; full path, name, and suffix. if left empty, assumes same file as filename
%                   lineFileNchans      {16}; how many channels in lineFilename
%                   chansToClean        {[]}; which channnels to clean from line influences (only neuronal)
%
% calls:            ParseArgPairs:      argument/value handling
%                   removesamples:      makes a new file without removed channels, then deletes original; d.t. aberrant recording/data transfer
%                   removechannels:     makes a new file without removed channels, then deletes original
%                   reorderchannels:    reorders the channels (steps over file)
%                   removeline:         removes the line from each channel (steps over file)
%                   removedcoffset:     removes the mean from each channel (steps over file)
%                   mergefiles:         concatenate head-to-tail (same durations)
%
%
% NOTE:             does not concatenate/split files
%                   that can be done using the following routines:
%
%                   mergefiles:         different channel count, same durations (e.g. recorded simultaneously on two files/systems)
%                   concatfiles:        same channel count, different durations (e.g. recorded sequentially on a given system)
%                   partition:          the reverse action of concatfiles (split a single big file into multiple shorter ones)
%                   extractfile:        the workhorse of partition
%
% xml files:        any modification of the channel order or number of
%                   channels automatically renders the xml file irrelevant.
%                   the new xml must contain:
%                       Number of channels:         nchans + nchansSecond
%                       Sampling rate:              same as original (e.g. 20000)
%                       Resolution (bits):          same as original (16)
%                       Initial offset:             same as original (0)
%                       Voltage range:              same as original (e.g. 2.45, or 49)
%                       Amplification:              same as original (e.g. 192, or 3840)
%
% NewOrder:
%   BUZ32 recorded by Intan RHD2132 with 8 analog input channels:
%   newOrder        = [ 21 20 22 18 31 17 30 16 24 25 28 26 23 27 19 29 8 7 6 5 3 4 2 12 11 10 13 9 14 0 15 1 32 : 39 ] + 1;
%   chansToClean    = 1 : 32; 
%
%   Edge32 recorded by Intan RHD2132:
%

% 12-aug-13 ES

% revisions
% 14-nov-17 (1) changed dcoFile default to -2
%           (2) added an option to back up the file
%           (3) modified argument handling
%           (4) added line cleaning option
%           (5) added integrity check and extra sample removal
% 03-dec-17 (1) finished writing removesamples.m and its call
%           (2) organized workflow, tested with files
% 05-dec-17 (1) integrated removeline, tested
%           (2) back up xml file

% to do:
% (1) update xml file
% (2) add two-*dat merging

function [ rc msg ] = massagedatfile( filename, varargin )

% constants
verbose         = 1;
nbytes          = 2;        % files are binary, int16
nsteps          = 7;

% initialize
rc              = zeros( nsteps, 1 );
et              = zeros( nsteps, 1 );
msg             = cell(  nsteps, 1 );
t0all           = clock;

% arguments
nargs = nargin;
if nargs < 1, error( 'filename required' ), end
[ makeBackup, nChannelsOrig, checkIntegrity ...
    , newOrder, chansToRemove ...
    , lineChannel, lineFilename, lineFileNchans, chansToClean ...
    , filenameSecond, nchansSecond, precisionSecond ...
    , dcofilemode ] = ParseArgPairs(...
    { 'makeBackup', 'nChannelsOrig', 'checkIntegrity' ...
    , 'newOrder', 'chansToRemove' ...
    , 'lineChannel', 'lineFilename', 'lineFileNchans', 'chansToClean' ...
    , 'filenameSecond', 'nchansSecond', 'precisionSecond' ...
    , 'dcofilemode' } ...
    , {  0, [], 1 ...
    , [], [] ...
    , [], [], 16, [] ...
    , 'digitalin', 1, 'uint16' ...
    , -2 } ...
    , varargin{ : } );
if isempty( nChannelsOrig )
    if ~isempty( newOrder )
        nChannelsOrig = length( newOrder );
    else
        error( 'nChannelsOrig required' )
    end
end

%------------------------------------------------------------------------
% preparations
%------------------------------------------------------------------------

% file management
if verbose
    fileinfo = dir( filename );
    fprintf( 1, '%s: processing %s (%3.2fGB): \n', upper( mfilename ), filename, fileinfo.bytes / 2^30 )
end
if ispc
    sep = '\';
    cpCmd = 'copy';
    rmCmd = 'del';
    mvCmd = 'rename';
else
    sep = '/';
    cpCmd = 'cp';
    rmCmd = 'rm';
    mvCmd = 'mv';
end
[ path filebase suffix ] = fileparts( filename );
if ~isempty( filenameSecond )
    [ path2 filebase2 suffix2 ] = fileparts( filenameSecond );
    if isempty( path2 )
        if isempty( suffix2 )
            suf = suffix;
        else
            suf = suffix2;
        end
        filenameSecond = sprintf( '%s%s%s%s', path, sep, filebase2, suf );
    end
    if exist( filenameSecond, 'file' ) ~= 2
        filenameSecond = '';
        msg{ 7 } = 'cannot find second file';
        fprintf( '%s\n', msg{ 7 } )
    end
end


% channel management (if any channels are to be removed, remap newOrder)
if ~isempty( newOrder )
    rmv = ismember( newOrder, chansToRemove );
    neworder = newOrder( ~rmv );
    nchans = length( neworder );
    if sum( neworder > nchans )
        mat = [ unique( neworder )' ( 1 : nchans )' ];
        for i = 1 : nchans
            neworder( i ) = mat( mat( :, 1 ) == neworder( i ), 2 );
        end
    end
else
    nchans = nChannelsOrig;
end

%------------------------------------------------------------------------
% STEP 1: make a backup
%         write a new file and work on (overwrite) the old one
%------------------------------------------------------------------------
if makeBackup
    bkp = sprintf( '%s%s%s.orig%s', path, sep, filebase, suffix );
    if exist( bkp, 'file' ) == 2
        if verbose
            fprintf( 1, '\tStep 1: \tNOT \tstepping on existing backup file %s!!\n', bkp )
        end
    else
        if verbose
            fprintf( 1, '\tStep 1: Backing up...' )
            t0 = clock;
        end
        cmd = sprintf( '!%s %s %s', cpCmd, filename, bkp );
        eval( cmd )
    end
    xmlfilename = [ path '/' filebase '.xml' ];
    xmlbkp      = [ path '/' filebase '.orig.xml' ];
    if exist( xmlfilename, 'file' ) == 2
        cmd = sprintf( '!%s %s %s', cpCmd, xmlfilename, xmlbkp );
        if exist( xmlbkp, 'file' ) == 2
            if verbose
                fprintf( 1, '\tStep 1: \tNOT \tstepping on existing backup xml file %s!!\n', xmlbkp )
            end
        else
            eval( cmd )
        end
    end
    if verbose
        et( 1 ) = etime( clock, t0 );
        fprintf( 1, 'done (%0.3g sec).\n', et( 1 ) )
    end
else
    if verbose
        fprintf( 1, '\tStep 1: \tNo backup made.\n' )
    end
end

%------------------------------------------------------------------------
% STEP 2: check file integrity (integer multiples)
%         write a new file and erase old one
%         removesamples.m
%------------------------------------------------------------------------
if checkIntegrity
    info = dir( filename );
    nsamples = info.bytes / nbytes / nChannelsOrig;
    if isequal( nsamples, round( nsamples ) )
        rc( 2 ) = 0; 
        msg{ 2 } = 'file integrity intact';
        if verbose
            fprintf( 1, '\tStep 2: Checked for aberrant samples... %s.\n', msg{ 2 } )
        end
    else
        if checkIntegrity == -1
            msg{ 2 } = sprintf( 'incorrect nchans (%d) for file %s', nchans, filename );
            return
        end
        % remove the extra data
        sampsToRemove = ( floor( nsamples ) * nChannelsOrig + 1 ) : nsamples * nChannelsOrig;
        if sampsToRemove > 0
            if verbose
                t0 = clock;
                fprintf( 1, '\tStep 2: Removing %d aberrant samples...', length( sampsToRemove ) )
            end
            [ rc( 2 ) msg{ 2 } ] = removesamples( filename, nChannelsOrig, sampsToRemove );
            if ~isempty( msg{ 2 } )
                fprintf( '%s\n', msg{ 2 } )
            end
            if verbose
                et( 2 ) = etime( clock, t0 );
                fprintf( 1, 'done (%0.3g sec).\n', et( 2 ) )
            end
        end
    end
else
    if verbose
        fprintf( 1, '\tStep 2: \tIntegrity not checked.\n' )
    end
end

%------------------------------------------------------------------------
% STEP 3: remove extra channels
%         write a new file and erase old one
%         removechannels.m
%------------------------------------------------------------------------
if isempty( chansToRemove )
    if verbose
        fprintf( 1, '\tStep 3: \tNo channels to remove.\n' )
    end
else
    if verbose
        t0 = clock;
        fprintf( 1, '\tStep 3: Removing %d channels... ', length( chansToRemove) )
    end
    [ rc( 3 ), msg{ 3 } ] = removechannels( filename, nChannelsOrig, chansToRemove );
    if ~isempty( msg{ 3 } )
        fprintf( '%s\n', msg{ 3 } )
    end
    if verbose
        et( 3 ) = etime( clock, t0 );
        fprintf( 1, 'done (%0.3g sec).\n', et( 3 ) )
    end
end

%------------------------------------------------------------------------
% STEP 4: reorder the remaining channels (e.g. according to geometry)
%         steps on file
%         reorderchannels.m
%------------------------------------------------------------------------
if isempty( newOrder )
    if verbose
        fprintf( 1, '\tStep 4: \tNo channels to reorder.\n' )
    end
else
    if verbose
        t0 = clock;
        fprintf( 1, '\tStep 4: Reordering channels... ' )
    end
    [ rc( 4 ), msg{ 4 } ] = reorderchannels( filename, neworder );
    if ~isempty( msg{ 4 } )
        fprintf( '%s\n', msg{ 4 } )
    end
    if verbose
        et( 4 ) = etime( clock, t0 );
        fprintf( 1, 'done (%0.3g sec).\n', et( 4 ) )
    end
end

%------------------------------------------------------------------------
% STEP 5: remove DC offset from those channels (compute for entire file)
%         steps on file
%         removedcoffset.m
%------------------------------------------------------------------------
if dcofilemode == -2
    if verbose
        fprintf( 1, '\tStep 5: \tNot removing DC.\n' )
    end
else
    if verbose
        t0 = clock;
        fprintf( 1, '\tStep 5: Removing DC... ' )
    end
    [ rc( 5 ), msg{ 5 } ] = removedcoffset( filename, nchans, dcofilemode );
    if ~isempty( msg{ 5 } )
        fprintf( '%s\n', msg{ 5 } )
    end
    if verbose
        et( 5 ) = etime( clock, t0 );
        fprintf( 1, 'done (%0.3g sec).\n', et( 5 ) )
    end
end

%------------------------------------------------------------------------
% STEP 6: merge data from other files
%         create a new file composed of two other files
%         mergefiles.m
%------------------------------------------------------------------------
if isempty( filenameSecond )
    if verbose
        fprintf( 1, '\tStep 6: \tNot merging second file.\n' )
    end
else
    % preps
    if verbose
        t0 = clock;
        fprintf( 1, '\tStep 6: Merging files... ' )
    end
    sourcefile1 = filename;
    sourcefile2 = filenameSecond;
    newfile     = sprintf( '%s%s%s_new%s', path, sep, filebase, suffix );
    nchansFirst = nchans;
    nchansMerge = [ nchansFirst nchansSecond ];
    precision   = { 'int16', precisionSecond };
    castmode    = 0;
    % actual call
    [ rc6 msg{ 6 } ] = mergefiles( sourcefile1, sourcefile2, newfile, nchansMerge, precision, castmode );
    % organize
    rc( 6 ) = mean( rc6 );
    if ~isempty( msg{ 6 } )
        fprintf( '%s\n', msg{ 6 } )
    end
    if rc( 6 ) == 0
        nchans  = sum( nchansMerge );
        cmd = sprintf( '!%s %s', rmCmd, filename );
        eval( cmd )
        cmd = sprintf( '!%s %s %s', mvCmd, newfile, filename );
        eval( cmd )
        if verbose
            et( 6 ) = etime( clock, t0 );
            fprintf( 1, 'done (%0.3g sec).\n', et( 6 ) )
        end
    end
end


%------------------------------------------------------------------------
% STEP 7: remove line influences from all requested channels
%         steps on file
%         removeline.m
%------------------------------------------------------------------------
if isempty( lineChannel ) || isempty( filenameSecond ) % temporary 
    if verbose
        fprintf( 1, '\tStep 7: \tNot removing mains influences.\n' )
    end
else
    if verbose
        t0 = clock;
        fprintf( 1, '\tStep 7: Removing mains influences... ' )
    end
    if floor( lineChannel ) > nchans
        msg{ 7 }    = 'file integrity intact';
        rc( 7 )     = -1;
    else
        lineGraphics = 1;
        [ rc( 7 ) msg{ 7 } ] = removeline( filename, nchans, lineChannel, 'chansToClean', chansToClean, 'graphics', lineGraphics );
    end
    if ~isempty( msg{ 7 } )
        fprintf( '%s\n', msg{ 7 } )
    end
    if verbose
        et( 7 ) = etime( clock, t0 );
        fprintf( 1, 'done (%0.3g sec).\n', et( 7 ) )
    end
end

if verbose
    fileinfo = dir( filename );
    etAll = etime( clock, t0all );
    fprintf( 1, 'All done (%3.2fGB; %0.3g sec).\n', fileinfo.bytes / 2^30, etAll )
end

return

% EOF

% call examples:

% BUZ32 w/ 4 diodes:
filename            = '/Users/eranstark/Documents/slab/data/m589_171130_103429/all_in_one.dat';
filename            = '/Users/eranstark/Documents/slab/data/m589_171130_194716/all_in_one.dat';
makeBackup          = 1;            % at least for the testing phase
% 32 neuronal, 4 movement (X, Y, theta, AMZ), 4 diodes
mat                 = [ 21 20 22 18 31 17 30 16; 24 25 28 26 23 27 19 29; 8 7 6 5 3 4 2 12; 11 10 13 9 14 0 15 1; 32 : 39 ]';
newOrder            = mat( : )' + 1;
lineChannel         = 41.6;
[ rc msg ] = massagedatfile( filename, 'makeBackup', makeBackup, 'newOrder', newOrder, 'lineChannel', lineChannel );

% 3 tetrodes w/ 1 diode:
filename            = '/Users/eranstark/Documents/slab/data/m746_171130_172237/all_in_one.dat';
filename            = '/Users/eranstark/Documents/slab/data/m746_171128_132518/all_in_one.dat';
makeBackup          = 1;
% 12 neuronal, 4 movement (X, Y, theta, AMZ), 1 diode
newOrder            = 1 : 17; 
%nChannelsOrig       = 17;
lineChannel         = 18.6;
[ rc msg ] = massagedatfile( filename, 'makeBackup', makeBackup, 'newOrder', newOrder, 'lineChannel', lineChannel );
%[ rc msg ] = massagedatfile( filename, 'makeBackup', makeBackup, 'nChannelsOrig', nChannelsOrig, 'lineChannel', lineChannel );

