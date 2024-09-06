function [convres,frex,trimmedT] = MyCWT(epoched_signal,num_trials,timevec,trim_time,srate)
    arguments
        epoched_signal
        num_trials
        timevec
        trim_time = 1500;
        srate = 1000;
    end
    
% Initialize parameters of signal
signal_length = length(epoched_signal);
num_elecs = size(epoched_signal,3);

% Create wavelet family
[wavelets,frex,nFrex,pnts] = MakeWaveletsWaveNumber(srate); 

% Ensure linear convolution
nConv = signal_length + pnts - 1;
halfK = floor(pnts/2);

% FFT of signal
sigX = fft(epoched_signal,nConv);
   
% Initialize result of convolution matrix
convres = zeros(nFrex,nConv,num_trials,num_elecs);

% Convolution per frequency
for fi = 1:nFrex
    
    % FFTs of wavelets
    waveX = fft(wavelets(fi,:),nConv).';
    
    % Element-wise multiplication and inverse FFT to get back into time
    % domain 
    convres(fi,:,:,:) = ifft( sigX.*waveX);
end

% Remove edges of convolution
convres = convres(:,halfK+1:end-halfK,:,:);

% Trimming buffers off of trials (was always 1500 ms)

convres = convres(:,1+trim_time:end-trim_time,:,:);

trimmedT = timevec(1+trim_time:end-trim_time); % New time vector
end