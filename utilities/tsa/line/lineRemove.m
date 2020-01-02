% lineRemove            clean analogue signal from line influence
%
% call:                 [ X xlta wvec ] = lineRemove( X, LINE, xlta, DC, MA, nmax, w )
%
% gets:                 X           input signal (vector)
%                       LINE        line (zero-crossing) times (same sampling rate as X)
%                       xlta        {[]}    IC for the line-triggered average
%                                           may be overloaded (two rows - xlta, wvec)
%                       nmax        {[]}    maximal possible length for xlta
%                       DC          {[]}    value to add after removing xlta
%                                           defaults to <X>
%                       MA           1      uses a running average
%                                   {0}     uses a weighted average
%                       w           {0.999} update weight
%
% returns:              X           signal vector, without the piece-wise
%                                   linear effect of the interference
%                       xlta        the effect of the interference
%                       wvec        weight vector
%
% calls:                nothing
%
% algorithmic notes:
%
% (1) choice of update method:
%   to track variations in the noise (the interference), the MA method is
%   optimal. However, it also tracks variations in the signal. For instance,
%   if there is a MA of 60 seconds and there is a short deflection (e.g.
%   pulse) of 60 ms, the contamination factor is 1:1000. But if there is a
%   pulse train with 50% duty cycle, the contamination will be 50%. In the
%   case of a zero-baseline signal with a pulse train of amplitude V, the
%   mean of the xlta will thus be V/2 (instead of close to zero in the
%   weighted average method), which will induce a huge shift in the signal.
%
% (2) signal mean:
%   in case of a zero-mean raw signal, the w-method will work perfectly even in
%   the presence of temporary deflections. However, in a DC-coupled signal
%   with a mean M, the mean of the xlta will also be M. This will result in
%   the line-removed waveform having a zero-mean. If the raw signal has
%   zero-baseline but occasional deflections of amplitude N (e.g. the mean
%   is 0.5 and the max value is 1, i.e. a duty cycle d of 50%), then
%   the mean of the xlta will be 0.5. Subtracting will yield a zero-mean
%   signal, with a baseline of -0.5 and a max value of 0.5.    
%   this problem is solved by adding the DC to the xlta. Thus, the formula
%   is
%
%       y(t) = x(t) - n(t) + myu(x)
%
%   in the above example, adding the mean of 0.5 will shift the signal back
%   to have zero baseline and max of 1. 
% 
%   the situation is less than ideal when the duty cycle is small and the
%   method is MA. For instance, assume dc is 1%, the signal baseline M is 
%   e.g. -1, and the max deflection N is e.g. 2. Then, the
%   signal mean is 
%           M * ( 1 - dc ) + ( M + N ) * dc = 
%           M + N * dc = 
%           -0.98
%   if the instantaneous xlta is for a period in which the mean is M, then
%   the error will be 
%
%           -M + [ M + N * dc ] = N * dc
%
%   which, in our example, 2*0.01 = 0.02 positive shift. Of course, this is
%   still much better than not removing the mean
%
% (3) the I.C.should be computed for a period without any major deflections. 
%
% (4) rarely-sampled points:
%   are treated by a weighted average between the sample, the last 
%   fully-sampled point, and the cyclically next fully-sampled point. note
%   that this is relevant only for the w-average mode
%
% see also:             lineDetect, lineClean
%                       removeline (wrapper)

% 04-apr-12 ES

% revisions
% 03-dec-17 (1) optional MA version (about 7000% slower)
%           (2) optional external ICs to support full block design
% 04-dec-17 (1) speeded up MA so it is only 10% slower; made this the default
%           (2) removed the Fs argument (IC determined by w and cycle duration only)
% 05-dec-17 (1) added nmax (external argument)
% 06-dec-17 (1) first cycle subtracted backwards
%           (2) w the default method
%           (3) added algorithmic notes
%           (4) added DC addition
% 07-dec-17 (1) modified method to verify X within range of LINE
%           (2) first clean, then update 
%           (3) weighted mean for rarely-sampled points
%           (4) added DC to first cycle

function [ X xlta ] = lineRemove( X, LINE, xlta, DC, nmax, MA, w )

% arguments
nargs = nargin;
if nargs < 3 || isempty( xlta )
    xlta = [];
end
if nargs < 4 || isempty( DC )
    DC = [];
end
if nargs < 5 || isempty( nmax )
    nmax = [];
end
if nargs < 6 || isempty( MA )
    MA = 0;
end
if nargs < 7 || isempty( w )
    w = 0.999; 
end

% preps
n       = ceil( 1 / ( 1 - w ) );
LINE    = LINE( : );
i       = find( LINE > 0, 1 );
if i > 1
    LINE( 1 : ( i - 1 ) ) = [];
end
X       = X( : ).';
nX      = length( X );
dL      = diff( LINE );
T       = max( dL );
WIN     = [ 0 T * n ];
WIN( 2 ) = min( max( LINE ) - T, T * n);

% initialization
if isempty( DC )
    DC = mean( X );
end
if isempty( xlta )
    lineIC  = LINE( find( LINE >= WIN( 1 ), 1 ) : find( LINE <= ( WIN( 2 ) - T ), 1, 'last' ) );
    rows    = length( lineIC );
    idx     = repmat( lineIC, 1, T ) + repmat( 0 : T - 1, rows, 1 );
    xmat    = X( idx );
    xlta    = mean( xmat, 1 );
    idx     = [];
    wvec    = size( xmat, 1 ) * ones( 1, size( xmat, 2 ) );
else
    if size( xlta, 2 ) == 2
        wvec = xlta( 2, : );
        xlta = xlta( 1, : );
    else
        wvec = 1/w * ones( 1, size( xlta, 2 ) );
    end
    xlta    = [ xlta zeros( 1, T - length( xlta ) ) ];
	xmat    = repmat( xlta, [ n 1 ] );
    wvec    = [ wvec zeros( 1, T - length( xlta ) ) ];
end
if size( xmat, 2 ) < nmax
    nz      = nmax - size( xmat, 2 );
    xlta    = [ xlta zeros( 1, nz ) ];
    xmat    = [ xmat zeros( size( xmat, 1 ), nz ) ];
end
[ n m ]     = size( xmat );
w           = 1 - 1/n;

% clean first cycle
if LINE( 1 ) > 1
    t1              = 1;
    t2              = LINE( 1 ) - 1;
    nx              = t2 - t1 + 1;
    idx             = length( xlta ) + [ ( 1 - t2 ) : 0 ];
    x               = X( t1 : t2 );
    X( t1 : t2 )    = x - xlta( idx ) + DC;
end
% adaptive subtraction
for i = 1 : length( LINE ) - 1
    t1 = LINE( i );
    t2 = LINE( i + 1 ) - 1;
    if t2 > nX
        t2 = nX; 
    end
    if t2 - t1 + 1 > T
        t2 = t1 + T - 1;
    end
    % clean
    x               = X( t1 : t2 );
    nx              = t2 - t1 + 1;
    idx             = 1 : nx;
    X( t1 : t2 )    = x - xlta( idx ) + DC;
    % update average
    if MA
        ihat            = mod( i - 1, n ) + 1;
        xz              = [ x zeros( 1, m - nx ) ];
        xlta            = xlta - ( 1 - w ) * xmat( ihat, : ) + ( 1 - w ) * xz;
        xmat( ihat, : ) = xz;
    else
        xlta( idx )     = w * xlta( idx ) + ( 1 - w ) * x;
    end
    % handle rarely-sampled points
    fidx            = ( nx + 1 ) : m;
    bidx            = [ 1 nx ];
    wvec( fidx )    = max( wvec( fidx ) - 1, 0 );
    num             = sum( xlta( bidx ) .* wvec( bidx ) ) + xlta( fidx ) .* wvec( fidx );
    den             = sum( wvec( bidx ) ) + wvec( fidx );
    xlta( fidx )    = num ./ den;
    if t2 >= nX
        break
    end
end

return

% EOF
