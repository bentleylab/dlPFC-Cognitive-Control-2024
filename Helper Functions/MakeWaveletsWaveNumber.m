function [wavelets,frex,nFrex,pnts] = MakeWaveletsWaveNumber(srate,length_wavelet,width)
    arguments
        srate = 1000;
        length_wavelet = 2;
        width = 6
        end

% Parameters

frex = 2.^([10:76]/10); % frequencies to compute CWT, log-spaced, 10 scales per octave, 2-194 Hz
nFrex = length(frex); % number of scales
wtime = -length_wavelet:1/srate:length_wavelet; % time vector for wavelet kernel centered 0 s ( -2 to 2 s, step size 1 ms)
pnts = length(wtime); % number of samples in kernel
wavelets = NaN(nFrex,pnts); % initialize wavelets


% Create complex Morlet wavelet family
for wi = 1:nFrex
    
    % Define parameters for time-frequency trade-off
    sf = frex(wi)/width; % spread factor or bandwidth in frequency domain for given frequency
    s = 1/(2*pi*sf); % time domain standard deviation of the Gaussian
    
    % Create complex wavelet
    A = 1/sqrt(s*sqrt(pi)); % normalization factor to ensure all scales have the same energy (unity)
    wavelets(wi,:) = A*exp(2*pi*1i*frex(wi).*wtime) .* exp(-wtime.^2 ./ (2*s^2)); % taper complex sine wave with gaussian to construct complex Morlet wavelet and normalize its energy
end
end