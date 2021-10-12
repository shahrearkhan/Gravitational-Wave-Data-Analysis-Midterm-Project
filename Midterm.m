%Authors: Wendy Mendoza & Shahrear Khan Faisal 
%Due Date: October 18, 2021 
%Mid-term Project PHYS 5387 

%(i):Loading Data for LIGO S5 strain 
data1 = load('dat1.mat');
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
[pxx, f] = pwelch(dat1,[],[],fs/4,fs); %Recheck pwelch 
 
%(iv):Applying the lowpass fiter to data1 
Wn = 2048/fs/2;  %Cutoff frequency of 0.025
[b,a] = butter(6, Wn, 'low'); %Bandpass of 6th order Butterworth filter 

%Filter the data 
dat_low_pass = filter(b,a,dat1); 
figure; 
plot(tdata,dat_low_pass); 
title('Data stream dat.mat')
xlabel('Frequency');
ylabel('Amplitude');
%Retains frequencies below a given cut-off and
%eliminate the higher frequencies 

%(v):Resample the data to a lower sampling frequency
p = 1; %Resampling factors
q = 4; %Resampling factors 
dat_low_pass1=resample(dat_low_pass,1,4); %change sampling rate
tdata1=(0:(length(dat_low_pass1)-1))*4/(1*fs); 
%new sample frequency 

%(vi)pwelch the new sampling rate 
[pxx1,f1]=pwelch(dat_low_pass1,[],[],[],1024); 
figure; 
plot(f1,10*log10(pxx1)); 

%(vii):Whiten the data
pxx_median_est = rngmed2(pxx1,256);  

%Median reduces the uncertainty while ignore the biased estimation if there
%is a strong sigal in addition to the noise. 
%Inverse the median 
freq = 0:0.25/fs:1;
%This function desugn filters with frequency inverse and the magnitude characterisitics of
%the white noise which is the containing vectors desires in the freq.  
bfilt=fir2(500,freq',1./sqrt(pxx_median_est) ); 
%fir2 is the filter order specidfied as an integer scalar. For confifuration with passband at the Nyqiost
%frequency fir2 always uses even order. 
 







