%Authors: Wendy Mendoza & Shahrear Khan Faisal 
%Due Date: October 18, 2021 
%Mid-term Project PHYS 5387 
addpath 'https://drive.matlab.com/files/'
%(i):Loading Data for LIGO S5 strain 
load S5_878486500_878486600_DARM.mat
% data1 = load('dat1.mat');
data1 = dat1;
fs = 16384; %Sampling frequency 

tdata = (0:length(data1)-1/fs); %Total length of the data time axis in seconds  
size(tdata); %total length of the data in seconds is 1 x 3145728
timeVec = 0:1/fs:(tdata-1)/fs; %time vector 


%(ii):Estimation of the power spectrum of dat1 using pwelch methods. 
%Power spectrum density measure the signal power versus frequency(width).
%PDS shows the strong and weak frequencies variations(energy). For the PSD
%is use the terms of Welch method to approach a spectral density estimation. 
% Using PWELCH methods to calculate the spectral density.   
[pxx, f] = pwelch(data1,[],[],[],fs); 

% (iii):Plotting the time series and the power of spectrum. 
figure(1)
subplot(2,1,1);
plot(tdata,dat1); 
xlabel('Time(sec)');
ylabel('Amplitude');
title('Time Series');
subplot(2,1,2);
plot(f, 10*log10(pxx)); 
xlabel('Frequency(Hz)');
ylabel('PSD (dB/Hz)');
title('Power Spectral Density(PWELCH)')
 
%(iv):Applying the lowpass fiter to data1 
Wn = 2048/fs/2;  %Cutoff frequency of 0.025
[b,a] = butter(6, Wn, 'low'); %Bandpass of 6th order Butterworth filter 
%Filtering the data 
dat_low_pass = filter(b,a,data1); 
[pxx1, f1] = pwelch(dat_low_pass,[],[],[],fs);
figure(2)
subplot(3,1,1);
plot(tdata,dat_low_pass); 
xlabel('Time(seconds)');
ylabel('Amplitude');
title('Filtered Time Series'); 
subplot(3,1,2);
plot(f1, 10*log10(pxx1))
xlabel('Frequency (Hz)');
ylabel('Magnitude(dB)');
title('Magnitude Response (dB)'); 



%The low pass filter retains frequencies below a given cut-off which means
% it eliminate the higher frequencies to allow the lower frequencies to pass
%through. 

%(v):Resample the data to a lower sampling frequency
p = 1; %Resampling factors
q = 4; %Resampling factors 
dat_low_pass1=resample(dat_low_pass,1,4); %change sampling rate
%new sample frequency will be 4096Hz 
fs1 = fs * p/q; 
tdata1 = (0:length(data1)-1/fs1); %duration in seconds 
timeVec1 = 0:1/fs:(tdata1-1)/fs; %time vector 

%(vi)pwelch the new sampling rate 
[pxx2, f2] = pwelch(dat_low_pass1,[],[],[],fs1); 
figure(3);   
plot(f2, 10*log10(pxx2)); 
xlabel('Frequency(Hz)');
ylabel('PSD (dB/Hz)');
title('Power Spectral Density with Resample Data');
% Comparing the power spectrum of the data before and after low pass
% filering, we can see that before the data had frequency components above
% 8 kHz, but after applying the filter the frequency components above 2048
% Hz have been discarded, which is expected, as out cut of frequency of the
% butterworth filter was 2048 Hz. 

%(vii):Whiten the data
pxx_median_est = rngmed2(pxx,256);  
%Median reduces the uncertainty while ignore the biased estimation if there
%is a strong signal in addition to the noise.  

%Inverse the median 
freq=0:0.125/(fs*4):1;
%This function design filters with frequency inverse and the magnitude characterisitics of
%the white noise which is the containing vectors desires in the freq.  
N = 500;
bfilt=fir2(500,freq',1./sqrt(pxx_median_est)); %N is number of order, frequency, & magnitude
%fir2 is digital filter order specified as an integer scalar. For
%configuration with passband at the Nyquist. This remove of certain
%frequency band in the data. 
dat_whitened=fftfilt(bfilt,dat_low_pass); 
[pxx3, f3] = pwelch(dat_whitened,[],[],[],fs1); 

figure(4); 
subplot(4,1,1)
plot(tdata1,dat_whitened); 
xlabel('Time(sec)');
ylabel('Amplitude');
title('Time Series of Whitened Data');
subplot(4,1,2);
plot(f3, 10*log10(pxx3)); 
xlabel('Frequency(Hz)');
ylabel('PSD (dB/Hz)');
title('Magnitude Responds of Whiten Data')

%(viii) Applying the Bandpass filter to the whitened data 
[b,a]=butter(6,[100 1024]/2048);  
dat_whitened1=filtfilt(b,a,dat_whitened); 
[pxx4, f4] = pwelch(dat_whitened1,[],[],[],fs1); 

figure(5); 
subplot(5,1,1)
plot(tdata1,dat_whitened1); 
xlabel('Time(sec)');
ylabel('Amplitude');
title('Time Series of Bandpass Whiten Data');

subplot(5,1,2);
plot(f4, 10*log10(pxx4));
xlabel('Frequency(Hz)');
ylabel('PSD (dB/Hz)');
title('Magnitude Responds of Bandpass Whiten Data');
%The bandpassfilter retains frequencies between a given lower cut-off and a
%higher cut-off. % Why we bandpassed the whitened data - Because LIGO is most sensitive
% between 100 Hz and 1000 Hz. Outside this frequency range signals other
% than our desired signals (noise) prevail, which makes it nearly
% impossible to detect gravitational wave in those frequency ranges. 

%(ix):Resample the whitened data to 2018 HZ 
dat_low_pass1=resample(dat_low_pass,1,4);
nfft = 2018 * 4; 
[pxx5, f5] = pwelch(dat_low_pass1,[],[],nfft); 
figure(6);
plot(f5, 10*log10(pxx5)); 
xlabel('Frequency(Hz)');
ylabel('PSD (dB/Hz)');
title('Magnitude Responds of Resample Data to 2018 Hz'); 

y = f5;
x = pxx5;
[pks,locs,widths,proms] = findpeaks(y,x); 






