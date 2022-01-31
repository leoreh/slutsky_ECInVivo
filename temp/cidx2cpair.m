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

% % another approach - could not make it work
% % find cidx from cpairs
% if isempty(cidx)
%         tmp = sub2ind([nunits nunits], cpair(:, 1), cpair(:, 2));
%         cidx = tmp - round(tmp ./ nunits);
% end
% 
% % find cpairs from cidx
% if isempty(cpair)
%     [u1, u2] = ind2sub([nunits, nunits], cidx);
%     cpair = [u1 + round(cidx ./ nunits), u2];
% end

end