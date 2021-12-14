clear all
close all

%%%%% test interpolation methods %%%%%

Fs = 1000;
t = 0:1/Fs:1-1/Fs;
nx = length(t);
f0 = 1;
f1 = 100;
t1 = t(end);
x = chirp(t,f0,t1,f1,'linear', 180);
% x = randn(nx, 1);

nx2 = nx * 2;
nx3 = nx / 2;
idx1 = 1 : nx;
idx2 = 1 : (nx - 1) / (nx2 - 1) : nx;
idx3 = 1 : (nx - 1) / (nx3 - 1) : nx;

interpm = {'linear', 'nearest', 'spline'};

for i = 1 : length(interpm)
    x2(i, :) = interp1(idx1, x, idx2, interpm{i});
    x3(i, :) = interp1(idx1, x, idx3, interpm{i});
    
    r2(i, :) = fft(x2(i, :))';
    r3(i, :) = fft(x3(i, :))';
    
    p2 = abs(r2(i, :) / nx2);
    p1 = p2(1 : nx2 / 2 + 1);
    p1(2:end-1) = 2*p1(2:end-1);
    z2(i, :) = p1;
    
    p2 = abs(r3(i, :) / nx3);
    p1 = p2(1 : nx3 / 2 + 1);
    p1(2:end-1) = 2*p1(2:end-1);
    z3(i, :) = p1;
end

f = Fs * (0 : (nx / 2)) / nx;
f2 = Fs * (0 : (nx2 / 2)) / nx2;
f3 = Fs * (0 : (nx3 / 2)) / nx3;

r(i, :) = fft(x)';
p2 = abs(r(i, :) / nx);
p1 = p2(1 : nx / 2 + 1);
p1(2:end-1) = 2*p1(2:end-1);

figure
subplot(3,2,1)
plot(x, 'k')
xlabel('t [sec]')
axis tight

subplot(3,2,2)
plot(f, p1, 'k')
ylabel('|p1(F)|')
xlabel('f [Hz]')
axis tight

subplot(3,2,3)
plot(x2(1, :))
hold on
plot(x2(2, :))
plot(x2(3, :))
xlabel('t [sec]')
axis tight

subplot(3,2,4)
plot(f2, z2(1, :))
hold on
plot(f2, z2(2, :))
plot(f2, z2(3, :))
ylabel('|p1(F)|')
xlabel('f [Hz]')
legend('linear', 'nearest', 'spline')
axis tight

subplot(3,2,5)
plot(x3(1, :))
plot(x3(2, :))
plot(x3(3, :))
xlabel('t [sec]')
axis tight

subplot(3,2,6)
plot(f3, z3(1, :))
hold on
plot(f3, z3(2, :))
plot(f3, z3(3, :))
ylabel('|p1(F)|')
xlabel('f [Hz]')
axis tight


%----------------------------------------------------%
Fs = 1000;
t = 0:1/Fs:1-1/Fs;
nx = length(t);
f0 = 1;
f1 = 100;
t1 = t(end);
x = chirp(t,f0,t1,f1,'linear', 180);
% x = randn(nx, 1);
xs = x(1:2:end);

interpm = {'linear', 'nearest', 'spline'};

idx1 = 1 : length(xs);
idx2 = 1 : (length(xs) - 1) / (length(x) - 1) : length(xs);

for i = 1 : length(interpm)
    xl(i, :) = interp1(idx1, xs, idx2, interpm{i});
    
    r5(i, :) = fft(xl(i, :))';
    
    p2 = abs(r5(i, :) / length(x));
    p1 = p2(1 : length(x) / 2 + 1);
    p1(2:end-1) = 2*p1(2:end-1);
    z5(i, :) = p1;
    
end

f = Fs * (0 : (length(x) / 2)) / length(x);

r(i, :) = fft(x)';
p2 = abs(r(i, :) / length(x));
p1 = p2(1 : length(x) / 2 + 1);
p1(2:end-1) = 2*p1(2:end-1);

figure
subplot(2,1,1)
plot(x, 'k')
hold on
plot(xl(1, :))
plot(xl(2, :))
plot(xl(3, :))
xlabel('t [sec]')
legend('raw', 'linear', 'nearest', 'spline')
axis tight

subplot(2,1,2)
plot(f, p1, 'k')
hold on
plot(f, z5(1, :))
plot(f, z5(2, :))
plot(f, z5(3, :))
ylabel('|p1(F)|')
xlabel('f [Hz]')
axis tight
