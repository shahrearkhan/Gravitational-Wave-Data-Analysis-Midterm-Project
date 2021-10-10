%Authors: Wendy Mendoza & Shahrear Khan Faisal 
%Due Date: October 18, 2021 
%Mid-term Project PHYS 5387 

%(i):Loading Data for LIGO S5 strain 
data1 = load('dat1');
fs = 16384; %Sampling frequency 

tdata = (0:length(data1)-1/fs); %Total length of the data in seconds  
size(tdata) %not sure if this is right I got 1 1 answer for the data in seconds 

%(ii):Estimation of the power spectrum of dat1 using pwelch methods. 
%Power spectrum density measure the signal power versus frequency we will
%be using PWELCH methods to calculate the spectral density.   
[pxx, f] = pwelch(dat1,[],[],fs/4,fs);  

% (iii):Plotting the time series and the power of spectrum. 
figure;
subplot(2,1,1);
plot(f,10*log10(pxx));
xlabel('Time(sec)');
ylabel('Amplitude');
title('Time Series');
figure; 
subplot(2,1,2); 
[pxx, f] = pwelch(dat1,[],[],fs/4,fs); 
 
%(iv):Applying the lowpass fiter to data1 
Wn = 2048/fs/2;  %Cutoff frequency of 0.25
[b,a] = butter(6, Wn, 'low'); %Bandpass of 6th order Butterworth filter 

%Filter the data 
dat_low_pass = filter(b,a,dat1); 
figure; 
plot(dat_low_pass); 
title('Time series after Butterworth bandpass filter')
%The low pass filtered allow us to lower the level of overshoot noise 



