function [cpair, cidx] = cidx2cpair(nunits, cidx, cpair)

if size(cpair, 1) == 2
    cpair = cpair';
end

% create simulated mat of indices
mat_cidx = reshape(1 : numel(ones(nunits, nunits)), nunits, nunits);
mat_cidx(diag(mat_cidx)) = nan;
for iunit = 1 : nunits
    mat_cidx(:, iunit) = mat_cidx(:, iunit) - 1 * iunit;
end
mat_cidx(triu(true(nunits))) = mat_cidx(triu(true(nunits))) + 1;

% find cidx from cpairs
if isempty(cidx)
    cpair(cpair(:, 1) == cpair(:, 2), :) = [];
    for ipair = 1 : size(cpair, 1)
        cidx(ipair) = mat_cidx(cpair(ipair, 1), cpair(ipair, 2));
    end
end

% find cpairs from cidx
if isempty(cpair)
    cidx(isnan(cidx)) = [];
    cpair = nan(length(cidx), 2);
    for ipair = 1 : length(cidx)
        [cpair(ipair, 1), cpair(ipair, 2)] = find(mat_cidx == cidx(ipair));
    end
end

% 
% if isempty(cidx)
%     cidx = (cpair(:, 2) - 1) .* (nunits - 1) + cpair(:, 1);
% end
% 
% 
% nunits = 10;
% cpair = [1, 2];
% cidx12 = cpair(:, 1) + cpair(:, 2) .* (nunits - 1)
% cpair = [10, 1];
% cidx34 = cpair(:, 1) + cpair(:, 2) .* (nunits - 1)
% 
% cpair = [1, 2];
% temp = cpair(:, 2) * (nunits - 1);
% cidx12 = temp + cpair(:, 1)
% cpair = [10, 1];
% temp = cpair(:, 2) * (nunits - 1);
% cidx34 = temp + cpair(:, 1)
% 
% nunits = 10;
% bs = nunits - 1;
% u1 = 1;
% u2 = 2;
% cidx12 = ( u2 - 1 ) * bs + u1
% u3 = 10;
% u4 = 1;
% cidx34 = ( u4 - 1 ) * bs + u3
% 
% 
% for cidx = 1 : 11
% cpair = [];
% if isempty(cpair)
%     cpair(:, 2) = ceil(cidx ./ (nunits - 1));
%     cpair(:, 1) = cidx - (cpair(:, 2) - 1) .* (nunits)
% end
% end
% 
% 
% 
% if isempty(cpair)
%     cpair(:, 2) = ceil(cidx ./ (nunits - 1));
%     cpair(:, 1) = cidx - (cpair(:, 2) - 1) .* (nunits - 1);
% end
% 
% 
%  mat_cidx = reshape(1 : numel(ones(nunits, nunits)), nunits, nunits);
%  mat_cidx(diag(mat_cidx)) = nan;
%  for iunit = 1 : nunits
%      mat_cidx(:, iunit) = mat_cidx(:, iunit) - 1 * iunit;
%  end
%  mat_cidx(triu(true(nunits))) = mat_cidx(triu(true(nunits))) + 1;

end