clear all
close all

Fs = 10000;
T = 20;
t = 0:(1/Fs):T;

df = 20;
% df = 0;
ps = 0;

cf = zeros(1, length(df));
for i = 1:length(df)

f1 = 50;
f2 = f1 + df(i);

x1 = sin(2*pi*f1*t)';
x2 = sin(2*pi*f2*t + ps)';

x = x1 + x2;

% extract zero crossings from signal
phase0 = -pi/2 - 1/Fs;
phs = angle( hilbert( x1 ) );
LINE = find( [ 0; diff( phs > phase0 ) ] == 1 );

% bandpass = [10 100];
% mode = [];
% graphics = 0;
% LINE = lineDetect( x, Fs, bandpass, mode, phase0, graphics );


% remove PLI
xlta = [];
DC = [];
nmax = [];
MA = [];
w = [];
[ X, xlta ] = lineRemove( x, LINE, xlta, DC, nmax, MA, w );
X = X';
xlta = xlta';

% contamination factor
idxPLI = [1:length(xlta)];
cf(i) = sum(abs(xlta - x1(idxPLI)));

end

% graphics
figure
subplot(2,2,[1 2])
plot(t, x1)
hold on
plot(t, x2)
plot(t, x)
plot(t, X)
separators( LINE / Fs, [], [ 0 0.7 0 ] );
txtf1 = sprintf('f1 = %.0fHz', f1);
txtf2 = sprintf('f2 = %.0fHz', f2);
legend(txtf1, txtf2, 'f1 + f2', 'f1 + f2 - f1');
xlabel('time [s]')
ylabel('amplitude')
ylim([-2 2])
xlim([5 5.05])

subplot(2,2,3)
plot(idxPLI, x1(idxPLI))
hold on
plot(idxPLI, xlta)
xlim([0 1/f1*Fs])
ylim([-2 2])
txt = sprintf('contamination factor = %.2f', cf(i));
text(10, 1, txt);
legend('original PLI', 'extracted PLI')
xlabel('sample')
ylabel('amplitude')

subplot(2,2,4)
semilogx(df, cf)
xlabel('log(\Delta frequency)')
ylabel('Contamination factor')