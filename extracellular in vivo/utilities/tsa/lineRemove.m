% lineRemove            clean analogue signal from line influence
%
% call:                 [ X tsa wvec xmat] = lineRemove( X, LINE, tsa, DC, MA, TW, nmax, w )
%
% gets:                 X           input signal (vector)
%                       LINE        line (zero-crossing) times (same sampling rate as X)
%                       tsa         {[]}    IC for the line-triggered average
%                                           may be overloaded (two rows - xlta, wvec)
%                       nmax        {[]}    maximal possible length for xlta
%                       DC          {[]}    value to add after removing xlta
%                                           defaults to <X>
%                       MA           1      uses a running average
%                                   {0}     uses a weighted average
%                       TW           1      uses time warping to detect and remove PLI
%                                   {0}     uses clipping to detect and remove PLI
%                       w           {0.999} update weight
%                                   {1000}  window size for running average
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
% 09-jul-18 (1) desloping added
%           (2) debugmode added
%           (3.1) changed intialization of wvec according to sample CDF
%           (3.2) TEMPORARILY cancelled special update of rarely-sampled
%           (3.3) unclear if the deslope removes the need for the special handling
% 19-jul-18 (1) added time wraping by linear interpolation
%           (2) changed xlta to tsa

function [ X tsa wvec xmat ] = lineRemove( X, LINE, tsa, DC, MA, TW, nmax, w )

% arguments
nargs = nargin;
if nargs < 3 || isempty( tsa )
    tsa = [];
end
if nargs < 4 || isempty( DC )
    DC = mean( X );
end
if nargs < 5 || isempty( MA )
    MA = 0;
end
if nargs < 6 || isempty( TW )
    TW = 0;
end
if nargs < 7 || isempty( nmax )
    nmax = [];
end
if nargs < 8 || isempty( w )
    w = 0.999;
end

% preps
n       = ceil( 1 / ( 1 - w ) );
LINE    = LINE( : );
i       = find( LINE > 0, 1 );
if i > 1
    LINE( 1 : ( i - 1 ) ) = [];
end
X         = X( : ).';
nX        = length( X );
dL        = diff( LINE );
mT        = round(mean(dL));
T         = max( dL );
WIN       = [ 0 T * n ];
WIN( 2 )  = min( max( LINE ) - T, T * n);

%% initialization %%
if isempty( tsa )
    if TW
        xmat    = zeros(n, mT);
        wvec    = ones(mT, 1) * n;
    else
        xmat    = zeros(n, T);
        % handle rarely-sampled points by calculating weighted mean
        lineIC              = LINE( find( LINE >= WIN( 1 ), 1 ) : find( LINE <= ( WIN( 2 ) - T ), 1, 'last' ) );
        dLIC                = diff( lineIC );
        edges               = 0.5 : 1 : max( dLIC ) + 0.5;
        h                   = histc( dLIC, edges );
        h( end )            = [];
        bins                = ( edges( 1 : end - 1 ) + edges( 2 : end ) ) / 2;
        Fx                  = ( sum( h ) - cumsum( h ) );
        wvec                = [n; Fx];
        wvec(end)           = [];
    end
    for i = 1 : WIN(2) / T
        t1              = LINE(i);
        t2              = LINE(i + 1) - 1;
        if t2 > nX
            t2 = nX;
        end
        xmat(i,:)       = prepseg(X(t1 : t2), size(xmat, 2), TW);
    end
    tsa                 = sum(xmat, 1) ./ wvec';
else
    if size( tsa, 2 ) == 2
        wvec    = tsa( 2, : );
        tsa     = tsa( 1, : );
    else
        wvec    = 1/w * ones( 1, size( tsa, 2 ) );
    end
    tsa         = [ tsa zeros( 1, T - length( tsa ) ) ];
    xmat        = repmat( tsa, [ n 1 ] );
    wvec        = [ wvec zeros( 1, T - length( tsa ) ) ];
end
if size( xmat, 2 ) < nmax
    nz          = nmax - size( xmat, 2 );
    tsa         = [ tsa zeros( 1, nz ) ];
    xmat        = [ xmat zeros( size( xmat, 1 ), nz ) ];
end
[n m]         = size(xmat);

%% clean first cycle %%
% if LINE(1) > 1
%     t1              = 1;
%     t2              = LINE(1) - 1;
%     idx             = length(tsa) + [(1 - t2) : 0];
%     X(t1 : t2)      = X(t1 : t2) - tsa(idx);
% end

%% adaptive subtraction %%
for i = 1 : length( LINE ) - 1
    t1              = LINE( i );
    t2              = LINE( i + 1 ) - 1;
    nx              = t2 - t1 + 1;
    if t2 > nX
        t2 = nX;
        X(t1 : t2)  = X(t1 : t2) - tsa(1 : length(t2));
        return
    end
    if TW
        idx             = 1 : (mT - 1) / (nx - 1) : mT;
        tsax            = tsa(round(idx));
        X(t1 : t2)      = X(t1 : t2) - tsax;
    else
        idx             = 1 : nx;
        X(t1 : t2)      = X(t1 : t2) - tsa(idx);
    end
    %% update average %%
    if i > WIN(2) / T
        x                       = prepseg(X(t1 : t2), m, TW);
        if MA
            ihat                = mod(i - 1, n) + 1;
            tsa                 = tsa - (1 - w) * xmat(ihat, :) + (1 - w) * x;
            xmat(ihat, :)       = x;
            if ~TW
                wvec            = sum(xmat ~= 0)';
            end
        else
            if TW
                tsa             = w * tsa + (1 - w) * x;
            else
                tsa(idx)        = w * tsa(idx) + (1 - w) * x(idx);
                wvec(idx)       = wvec(idx) + ones(nx, 1);
            end
        end
    end
    if t2 >= nX
        break
    end
    %             if 0                %nx > mean( dL )
    %             fidx            = ( nx + 1 ) : m;
    %             bidx            = [ 1 nx ];
    %             wvec( fidx )    = max( wvec( fidx ) - 1, 0 );
    %             num             = sum( tsa( bidx ) .* wvec( bidx ) ) + tsa( fidx ) .* wvec( fidx );
    %             den             = sum( wvec( bidx ) ) + wvec( fidx );
    %             tsa( fidx )     = num ./ den;
end

%-------------------------------------------------------------------------%
    function xo = prepseg(x, m, TW)
        % deslope, ddc and warp
        nx          = length(x);
        dx          = (x(nx) - x(1)) / (nx - 1);
        slp         = fliplr(cumsum([x(1) ones(1, nx - 1) * dx]));
        xo          = x + slp;
        xo          = xo - xo(1);
        if TW
            idx     = 1 : (nx - 1) / (m - 1) : nx;
            xo      = xo(round(idx));
        else
            xo      = [xo zeros(1, m - nx)];
        end
    end

end
% EOF
