function mat = vec2nanmat(v, rowsz)

% 19 jun 22 LH  

v = v(:);

cnt = 1;
mat = cell(1, length(rowsz));
for irow = 1 : length(rowsz)
    mat{irow} = v(cnt : cnt + rowsz(irow) - 1, :);
    cnt = cnt + rowsz(irow);
end
mat = cell2nanmat(mat, 2);

% EOF


