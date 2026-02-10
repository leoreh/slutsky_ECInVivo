
% Test script to identify the warning ID for "perfect fit"
try
    % Create data that might trigger "perfect fit"?
    % Small variance, very accurate model
    x = (1:100)';
    y = x + randn(100,1)*0.0000001; % almost perfect
    tbl = table(y,x);
    
    warning('on', 'all');
    lastwarn(''); % clear
    
    % Trying to force the warning about "Initial estimate of residual noise standard deviation is almost zero"
    % This usually happens in mixed models when variance component is near zero.
    % Let's try fitting a GLME with a distribution that might struggle or be perfect.
    % The user saw this with Inverse Gaussian.
    
    % If we can't reproduce exactly, we can try to guess the ID or search common ones.
    % But let's try to fit something.
    try
        fitglme(tbl, 'y ~ x + (1|x)', 'Distribution', 'InverseGaussian', 'Link', 'Log');
    catch
    end
    
    [msg, id] = lastwarn;
    fprintf('Last Warning ID: %s\n', id);
    fprintf('Last Warning Msg: %s\n', msg);

catch ME
    fprintf('Error: %s\n', ME.message);
end
