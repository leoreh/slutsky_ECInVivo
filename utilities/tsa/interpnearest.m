    function xo = interpnearest(idx1, x, idx2)
        xo = zeros(length(idx2), 1);
        for j = 1 : length(idx2)
            [v y] = min(abs(idx2(j) - idx1));
            xo(j) = x(idx1(y));
        end
        xo = xo(:).'
    return