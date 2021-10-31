% signal; chirp amp-modulated by a sine wave
fs = 1250;
f = 20;
t = 1 : 1 / fs : 10;
mod = 5 * sin(2 * pi * f / 10 * t);
sig = mod .* chirp(t, f, t(end), f * 10);
    
h = hilbert(sig);
% instantaneous phase and amplitude
phase = angle(h);
amplitude = abs(h);
% instantaneous frequency
instfrq = fs / (2 * pi) * diff(medfilt1(unwrap(phase), 1));
% This is from bz_RippleStats. DID NOT UNDERSTAND bz_Diff
% unwrapped = unwrap(phase);
% frequency = bz_Diff(medfilt1(unwrapped,12),timestamps,'smooth',0);
% frequency = frequency/(2*pi);

figure
subplot(2, 1, 1)
plot(t, sig, 'k')
hold on
plot(t, real(h), 'k')
plot(t, amplitude)
plot(t, phase)
xlim([1 1.2])
subplot(2, 1, 2)
plot(t, sig, 'k')
hold on
yyaxis right
plot(t(1 : end - 1), instfrq)
xlim([1 10])

