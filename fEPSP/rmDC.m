function mat = rmDC(mat, art, fs)

% removes DC from baseline

for i = 1 : size(mat, 2)  
    bl = mean(mat(1 : art - 0.002 * fs, i));
    mat(:, i) = mat(:, i) - bl;    
end

end