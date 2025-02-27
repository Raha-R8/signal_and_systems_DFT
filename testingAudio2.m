
Fs = 8000; 
t = 0:1/Fs:1-1/Fs; 

f1 = 500; 
f2 = 1000; 
f3 = 1500; 
noise_amp = 0.2;

signal = sin(2*pi*f1*t) + 0.5*sin(2*pi*f2*t) + 0.3*sin(2*pi*f3*t);

normalAudio = signal/max(abs(signal));
plot(t, signal);

audiowrite('testingImp.wav', normalAudio, Fs);