
data = [1 2; 3 4];
centers = [1.5; 3.5];
scales = [0.5; 0.5];

disp('--- Normalize with Vector Center/Scale ---')
try
    % normalize(A, dim, method, ...)
    % dim=2 (rows)
    % center: vector of length size(A,1) ?
    N = normalize(data, 2, 'center', centers, 'scale', scales);
    disp(N)
    % Expected: (1-1.5)/0.5 = -1, (2-1.5)/0.5 = 1
    %           (3-3.5)/0.5 = -1, (4-3.5)/0.5 = 1
catch ME
    disp(ME.message)
end
