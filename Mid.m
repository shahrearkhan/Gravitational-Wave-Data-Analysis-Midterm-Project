%% Midterm Project
addpath 'D:\COURSE\Fall_2021\Gravitational Wave\Midterm_Project'

%% Step I
load S5_878486500_878486600_DARM.mat
% Sampling Frequency
fs = 16384; % Hz
% Length of the data in seconds
t = length(dat1)/fs; % seconds

%% Step II
% Apply pwelch to get the power spectrum
% Window length
w = 512;
novrlp = w/2;
[pxx, f] = pwelch(dat1, w, novrlp, [], fs);


%% Step III
% Plot data and power spectrum
figure;
subplot(2,1,1);
plot(dat1);
xlabel('time(sec)');
ylabel('amplitude');
title('Time Series');
subplot(2,1,2);
pwelch(dat1, w, novrlp, [], fs);


%% Step IV
% Apply low pass butterworth filter
[b,a] = butter(6, 2048/fs/2, 'low');
% Produce filtered data
dat_low_pass = filter(b,a,dat1);
% Plot the filtered time series and it's power spectrum
figure;
subplot(2,1,1);
plot(dat_low_pass);
title('Filtered Time Series');
xlabel('time(sec)');
ylabel('Amplitude');
subplot(2,1,2);
pwelch(dat_low_pass, w, novrlp, [], fs);


%% Step V
% Resample filtered data
dat_low_pass = resample(dat_low_pass,1,4);
figure;
plot(dat_low_pass);
% New sampling frequency
fs = fs*(1/4);


%% Step VI
% Calculate the power spectrum of dat_low_pass
[pxx, f] = pwelch(dat_low_pass, w, novrlp, [], fs);
figure;
pwelch(dat_low_pass, w, novrlp, [], fs);


%% Step VII
% Whiten the data using rngmed2
pxx_median_est = rngmed2(pxx, 256);

