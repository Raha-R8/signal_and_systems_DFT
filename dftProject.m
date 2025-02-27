% Specify the path to the audio file

% filePath = '/Users/raha/Desktop/test2.mp3';
filePath = '/Users/raha/Desktop/synthetic_audio_1sec.wav';

% Read the audio file
[audioData, sampleRate] = audioread(filePath);
sound(audioData,sampleRate)

%___________________________________________________________________
%these lines can be used to test the correctness of code
% sound(audioData, sampleRate);

% t = (0:length(audioData)-1) / sampleRate; % Time vector
% plot(t, audioData);
% xlabel('Time (s)');
% ylabel('Amplitude');
% title('Audio Signal');
%___________________________________________________________________

% downsample is only needed for audios that are two large to be handled 
% using myDft function, therefore we will reduce the samples to be able to
% handle it using myDft

 % Downsample the audio data (e.g., take every 10th sample)
 % downsampleFactor = 10;
 % audioData = audioData(1:downsampleFactor:end);
%___________________________________________________________________



 function X = myDFT(x)
     N = length(x);
     X = zeros(1, N); 
     for k = 0:N-1
         for n = 0:N-1
             X(k+1) = X(k+1) + x(n+1) * exp(-1j * 2 * pi * k * n / N);
         end
     end
 end

 function x = myIDFT(signalDFT)
    N = length(signalDFT);
    x = zeros(1, N);
    for n = 0:N-1
        for k = 0:N-1
            x(n+1) = x(n+1) + (1/N) * signalDFT(k+1) * exp(1j * 2 * pi * k * n / N);
        end
    end
end



function X = myHPF(signal, frequencies, cuttingFreq)
    %finding the index of cutoff frequency 
    cutoffIndex = find(frequencies >= cuttingFreq, 1);
    signal(1:cutoffIndex) = 0; 
    X = signal;
end


 function X = myLPF(signal,frequencies,cuttingFreq)
    cutoffIndex = find(frequencies >= cuttingFreq, 1);
    signal(cutoffIndex+1:end) = 0; 
    X = signal;
 end

 function X = myBPF(signal,frequencies,bandStartFreq,bandEndFreq)
    cutoffIndex1 = find(frequencies >= bandStartFreq, 1);
    cutoffIndex2 = find(frequencies >= bandEndFreq,1);
    signal(1:cutoffIndex1) = 0;
    signal(cutoffIndex2+1:end) = 0;
    X = signal;
 end


% applying DFT
X = myDFT(audioData);
test = X;

%using this to check the result 
% answer = fft(audioData);

% Computing the frequency axis
N = length(audioData);
frequencies = (0:N-1) * (sampleRate / N);


%applying filter
%X = myBPF(X,frequencies,400,1600);
X = myHPF(X,frequencies,6000);
%X = myLPF(X,frequencies,1000);


%applying inverse DFT
newAudioData = myIDFT(X);



% making sure data is real and absolute so that it can be played using
% sound function
newAudioData = real(newAudioData);
newAudioData = newAudioData / max(abs(newAudioData));



%playing the result
% if the audio is downsampled the new sample rate must be given to the
% sound function:
% sound(newAudioData, sampleRate / downsampleFactor);

sound(newAudioData, sampleRate);



% Plotting the magnitude of the frequency spectrum and time spectrum to
% check if the filters are working correctly
subplot(4,1,1)
plot(frequencies, abs(test));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum of the Audio Signal');


subplot(4,1,2)
plot(frequencies,abs(X));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum of filtered Signal');



subplot(4,1,3)
 t = (0:length(audioData)-1) / sampleRate; % Time vector
 plot(t, audioData);
 xlabel('Time (s)');
 ylabel('Amplitude');
 title('Audio Signal');


 subplot(4,1,4)
 t = (0:length(newAudioData)-1) / sampleRate; % Time vector
 plot(t, newAudioData);
 xlabel('Time (s)');
 ylabel('Amplitude');
 title('New Audio Signal');