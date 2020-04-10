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
            xo      = [x zeros(1, m - nx)];
        end
    end