function [dist_mat] = wv2dist(wv,dist_fun, assume_symmetry)
arguments
    wv (:,:) double {mustBeNumeric}
    dist_fun (1,1) function_handle = @(x,y) 1 - max(abs(xcorr(x,y,"normalized")));
    assume_symmetry (1,1) logical {mustBeNumericOrLogical} = true;
end

if assume_symmetry
    % collect matrix lower trig
    for iWv1 = size(wv,1):-1:1
        for iWv2 = iWv1:-1:1
            dist_mat(iWv1,iWv2) = dist_fun(wv(iWv1,:),wv(iWv2,:));
        end
    end
    
    % copy the matrix across the diag
    dist_mat = dist_mat' + dist_mat;
else
    % collect full matrix
    for iWv1 = size(wv,1):-1:1
        for iWv2 = size(wv,1):-1:1
            dist_mat(iWv1,iWv2) = dist_fun(wv(iWv1,:),wv(iWv2,:));
        end
    end
end
% EOF