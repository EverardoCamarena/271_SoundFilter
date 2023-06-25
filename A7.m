% Fast Fourier Transform
NFFT = 4096;    
OVERLAP = 0.5;

% dB scale
dB_scale = 80;  
weightOption = 1;

% provided sound file in
[data,Fs] = audioread('mysteriousSound.wav');
temp = 1;
signal = data(:,temp);

% erasing portions
wait = 1;
if wait>1
    signal = decimate(signal,wait);
    Fs = Fs/wait;
end
samples = length(signal);

%graph to show audio file
% averaged FFT spectrum
[Frequency, sensor] = FFTPK(signal,Fs,NFFT,OVERLAP);
% converting to dB scale (ref = 1)
spectrumdB = 20*log10(sensor);
% applying weight if needed
if weightOption == 1
    AdB = Afunc(Frequency);
    spectrumdB = spectrumdB+AdB;
    labels = ('Amplitude dB');
else
    labels = ('Amplitude dB y');
end
figure(1),plot(Frequency,spectrumdB,'b');grid
title(['Average FFT = ' num2str(Frequency(2)-Frequency(1)) ' Hz ']);
xlabel('Frequency (Hz)');ylabel(labels);

% weight
function A_dB = Afunc(f)
	n = ((12200^2*f.^4)./((f.^2+20.6^2).*(f.^2+12200^2).*sqrt(f.^2+107.7^2).*sqrt(f.^2+737.9^2)));
	r = ((12200^2*1000.^4)./((1000.^2+20.6^2).*(1000.^2+12200^2).*sqrt(1000.^2+107.7^2).*sqrt(1000.^2+737.9^2))) * ones(size(f));
	pondA = n./r;
	A_dB = 20*log10(pondA(:));
end

% noise removal
filter = data;
min = 0;
max = 14.99;
index = 1+fix(min*Fs:max*Fs-1);
flow = 500;
fhigh = 2500;
N = 4;

% signal filtered 
[b,a] = butter(N,2/Fs*[flow fhigh]);
filter(index,:) = filtfilt(b,a,data(index,:));

% filtered audio
audiowrite('original.wav',filter,Fs);

function  [Freq,FFTspecs] = FFTPK(signal, Fs, nfft, Overlap)
[samples,channels] = size(signal);
if samples<nfft
    temporary = zeros(nfft,channels);
    temporary((1:samples),:) = signal;
    signal = temporary;
    samples = nfft;
end
window = hanning(nfft);
window = window(:);
%Fast Fourier Transform and overlap 
 offset = fix((1-Overlap)*nfft);
 spectnum = 1+ fix((samples-nfft)/offset); 
    %Hanning and averaging scaling
    FFTspecs = 0;
    for i=1:spectnum
        start = (i-1)*offset;
        sw = signal((1+start):(start+nfft),:).*(window*ones(1,channels));
        FFTspecs = FFTspecs + (abs(fft(sw))*4/nfft);
    end
   
    FFTspecs = FFTspecs/spectnum; 
    if rem(nfft,2) 
        select = (1:(nfft+1)/2)';
    else
        select = (1:nfft/2+1)';
    end
FFTspecs = FFTspecs(select,:);
Freq = (select - 1)*Fs/nfft;
end