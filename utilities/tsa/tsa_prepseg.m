function xo = tsa_prepseg(x, m, tw)

x = x(:)';

% deslope, ddc and warp
nseg = length(x);
dx  = (x(nseg) - x(1)) / (nseg - 1);
slp = fliplr(cumsum([x(1) ones(1, nseg - 1) * dx]));
xo = x + slp;
xo = xo - xo(1);
if tw
    segIdx = 1 : (nseg - 1) / (m - 1) : nseg;
    xo = xo(round(segIdx));
else
    xo = [xo, zeros(1, m - nseg)];
end
end