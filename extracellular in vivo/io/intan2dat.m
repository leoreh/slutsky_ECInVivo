% recordings were made with the RHD2132 by Intan. the pcb omentics (10x2)
% was connected to the left side of the headstage omnetics (20x2) such that
% channels [16:31] were used. see
% http://intantech.com/RHD2132_RHD2216_amp_board.html

% params
fs = 30000;
ch = 16 : 31; % intan uses zero based mapping

% load
j = 1;
for i = ch
    if i < 10
        prefix = 'amp-A-00';
    else
        prefix = 'amp-A-0';
    end
    filename = [prefix num2str(i) '.dat'];

fin = fopen(filename, 'r');
data(j, :) = fread(fin, Inf, 'int16');
fclose(fin);

j = j + 1;
end

% inspect
idx = 5 * fs : 6 * fs;
for i = 1 : length(ch) 
    figure
    plot(idx / fs, data(i, idx))
    axis tight
end

length(data) / fs / 60;

% write and save
newname = 'bruce_140519.dat';
fout = fopen(newname, 'w');
fwrite(fout, data(:), 'int16');         
fclose(fout);

