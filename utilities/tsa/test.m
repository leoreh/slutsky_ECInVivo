% % %%%%% load digital input %%%%%
% % filename = 'D:\Google Drive\PhD\Equipment\PLI\mc4_180612_222728\digitalin.dat';
% % numDig = 16;
% % fid = fopen(filename, 'r');
% % fileinfo = dir(filename);
% % num_samples = fileinfo.bytes/2; % uint16 = 2 bytes
% % digitalin = fread(fid, num_samples, 'uint16');
% % fclose(fid);
% % dig = fliplr(rem(floor(digitalin*pow2(1-numDig:0)) ,2));
% % clear digitalin
% % clear fileinfo
% 
% clear all
% close all
% 
% nchans = 17;
% Fs = 20000;
% 
% filename = 'D:\Google Drive\PhD\Equipment\PLI\mc4_180612_222728\all_in_one.dat';
% m = memmapfile( filename, 'Format', 'int16');
% d = reshape( m.data, [ nchans length(m.Data) / nchans ] );
% 
% mode = 'fir';
% LINE = lineDetect( d(17,:), Fs, [], mode, [], [] );
% 
% MA = 0;
% TW = 1;
% X = double(d(3,:)) * 0.195;
% 
% tic
% [Y1 tsa1 wvec1] = lineRemove(X, LINE, [], [], MA, TW, [], []);
% toc
%   
% TW = 0;
% tic
% [Y2 tsa2 wvec2] = lineRemove(X, LINE, [], [], MA, TW, [], []);
% toc
% 
% figure
% idx = Fs*20 : Fs*20 + Fs / 10;
% subplot(2,1,1)
% plot(X(idx))
% hold on
% plot(Y1(idx))
% plot(Y2(idx))
% % line([1 length(idx)], [mean(int16(Y1)), mean(int16(Y1))], 'lineWidth', 2, 'Color', 'k')
% % line([1 length(idx)], [mean(Y2), mean(Y2)], 'lineWidth', 2, 'Color', 'k')
% legend('raw', 'w/ TA', 'wo/ TA')
% title('PLI removal')
% axis tight
% subplot(2,1,2)
% plot(tsa1)
% hold on
% plot(tsa2, 'lineWidth', 2, 'lineStyle', '--')
% % line([1 length(tsa1)], [mean(tsa1), mean(tsa1)], 'lineWidth', 2, 'Color', 'k')
% % line([1 length(tsa2)], [mean(tsa2), mean(tsa2)], 'lineWidth', 2, 'Color', 'k')
% legend('w/ TA', 'wo/ TA')
% title('TSA')
% axis tight
% 
% % i               = 87;
% % t1              = LINE( i );
% % t2              = LINE( i + 1 ) - 1;
% % nx              = t2 - t1 + 1;
% % x               = X( t1 : t2 );
% % idx1            = 1 : nx;
% % idx2            = 1 : (nx - 1) / (mT - 1) : nx;
% % xl              = interp1(idx1, x, idx2);
% % 
% % i               = 6;
% % t1              = LINE( i );
% % t2              = LINE( i + 1 ) - 1;
% % nx              = t2 - t1 + 1;
% % x2              = X( t1 : t2 );
% % idx1            = 1 : nx;
% % idx2            = 1 : (nx - 1) / (300 - 1) : nx;
% % xl2              = interp1(idx1, x2, idx2);
% % xl3              = interp1(idx1, x2, idx2, 'nearest');
% % 
% % figure
% % subplot(2,1,1)
% % plot(x)
% % hold on
% % plot(xl)
% % axis tight
% % subplot(2,1,2)
% % plot(x2)
% % hold on
% % plot(xl2)
% % plot(xl3, 'k')
% % axis tight
% % 
% % % dt = 0.00005;
% % % T = 0.02;
% % % t = ( dt : dt : T )';
% % % f = 50;
% % % x = sin( 2 * pi * f * t );
% % % m = length( x );
% % % idx0 = ( 1 : m )'; 
% % % 
% % % nx = 450; 
% % % idx1 = ( 1 : ( m - 1 ) / ( nx - 1 ) : m )';
% % % x1 = interp1(idx0,x,idx1);
% % % 
% % % nx = 350; idx2 = ( 1 : ( m - 1 ) / ( nx - 1 ) : m )';
% % % x2 = interp1(idx0,x,idx2);
% % % 
% % % figure
% % % subplot( 3, 1, 1 ), plot( x ); axis tight
% % % subplot( 3, 1, 2 ), plot( x1 ); axis tight
% % % subplot( 3, 1, 3 ), plot( x2 ); axis tight
% % 
% 
% % mT              = 400;
% % i               = 6;
% % t1              = LINE( i );
% % t2              = LINE( i + 1 ) - 1;
% % nx              = t2 - t1 + 1;
% % x               = X( t1 : t2 );
% % 
% % idx1            = 1 : nx;
% % idx2            = 1 : (nx - 1) / (mT - 1) : nx;
% % xt              = interp1(idx1, x, idx2, interpm);
% % 
% % 
% % out = zeros(length(idx2), 1);
% % for j = 1 : length(idx2)
% %     [v i] = min(abs(idx2(j) - idx1))
% %     out(j) = x(idx1(i));
% % end
% 

%%
clear all
close all

Fs = 20000; % sampling rate
T = 20; % length of the windon in seconds
t = 0:(1/Fs):T;
F1 = 40;
F2 = 50;
randnoise = randn(1, length(t)) / 10;
randnoise = zeros(1, length(t));
x(1, :) = 2*cos(2*pi*F1*t);                   % signal of interest
x(2, :) = 2*sin(2*pi*F2*t);                   % PLI
x(3, :) = x(1, :) + x(2, :) + randnoise;    % raw data

% Line = find(x(2, :) > -0.00001 & x(2, :) < 0.00001 & [0 diff(x(2, :))] > 0);
Line = [1 : 400 : length(t)];
MA = 0;
TW = 1;
[x(4, :) tsa wvec] = lineRemove(x(3, :), Line, [], [], MA, TW, [], []);


 

figure
plot(x(1, :))
hold on
plot(x(2, :))
plot(x(3, :))
plot(x(4, :), 'k', 'lineWidth', 1)
xlim([19 20] * Fs)

figure, plot(x(2,1:400)), hold on, plot(tsa, 'k'), axis tight

% x = x + abs(min(x(3, :)));

figure
plot((x(3, :) - x(4, :)), 'lineWidth', 1);
hold on
plot((x(3, :) - x(1, :)), '*', 'lineWidth', 0.1);
plot(x(2, :), 'k');
xlim([0 1000])
% axis tight


% cleancontamination = sum(abs(x(4, :) - x(3, :)));
% rawcontamination = sum(abs(x(2, :) - x(3, :)));
% accuracy = cleancontamination * 100 / rawcontamination;

rPLI = (raw - clean) / max(PLI);
PLI = PLI / max(PLI);
accuracy = (abs(abs(rPLI) - abs(oPLI))) * 100;



rPLI = (x(3, :) - x(4, :) - randnoise) / max(x(2, :));
nPLI = x(2, :) / max(x(2, :));
accuracy = abs(rPLI) - abs(oPLI);

figure, plot(accuracy)
hold on
% plot(x(3, :) - x(4, :) - randnoise)
% plot(x(2, :))
plot(rPLI)
plot(nPLI)
xlim([0 1000])

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x(1, :) = cos(2*pi*F1*t);                   % signal of interest
% x(2, :) = sin(2*pi*F2*t);                   % PLI
% x(3, :) = x(1, :) + x(2, :) + randnoise;    % raw data
% x(4, :) = linRemove(x(3, :));               % clean data
% 
% % If our algorithm is perfect, then x(3, :) - x(4, :) - randnoise == x(2, :)
% % Any difference between the two sides of this equation indicated that the signal is still infected with PLI.
% % One way to normalize this measurement between 0 and 1 is:
% 
% residualPLI = x(3, :) - x(4, :) - randnoise;
% abs(x(2, :) - residualPLI) / max(abs(x(2, :) - residualPLI));
% 
% % However, the difference 