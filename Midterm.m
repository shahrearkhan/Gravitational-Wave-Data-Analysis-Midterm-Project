%% Midterm Project
addpath 'D:\COURSE\Fall_2021\Gravitational Wave\Midterm_Project'


%% Step I
load S5_878486500_878486600_DARM.mat
% Sampling Frequency
fs = 16384; % Hz
% Length of the data in seconds
t = length(dat1)/fs; % duration in seconds
% Length of the data
len = t*fs; % duration*sampling frequency
% Time vector
timeVec =  0:1/fs:(len-1)/fs;


%% Step II
% Power Spectral Density is the distribution of power in different
% frequency components of the signal

% Apply pwelch to get the power spectrum
% Window length
w = 512;
novrlp = w/2;
% [pxx, f] = pwelch(dat1, w, novrlp, [], fs);
% [pxx, f] = pwelch(dat1, 4096, 2048, [], fs);
% plot(f, 10*log10(pxx));
[pxx,f] = pwelch(dat1, [], [], [], fs);
% pwelch(dat1, [], [], [], fs);



%% Step III
% Plot data and power spectrum
figure;
subplot(2,1,1);
plot(timeVec, dat1);
xlabel('time(sec)');
ylabel('amplitude');
title('Time Series');
subplot(2,1,2);
% pwelch(dat1, w, novrlp, [], fs);
pwelch(dat1, [], [], [], fs);


%% Step IV
% Apply low pass butterworth filter
[b,a] = butter(6, 2048/fs/2, 'low'); % The cutoff frequency here is 0.25 (Is it actually 2048 Hz?)
% Produce filtered data
dat_low_pass = filter(b,a,dat1); % Have to write down the filter transfer function
% Plot the filtered time series and it's power spectrum
figure;
subplot(2,1,1);
plot(timeVec,dat_low_pass);
title('Filtered Time Series');
xlabel('time(sec)');
ylabel('Amplitude');
subplot(2,1,2);
% pwelch(dat_low_pass, w, novrlp, [], fs);
pwelch(dat_low_pass, [], [], [], fs);

% Why did we low pass the data? - Because above a certain frequency (above
% 1 kHz) there is shot noise, which we want to discard


%% Step V
% Resample filtered data
dat_low_pass = resample(dat_low_pass,1,4);
% New sampling frequency
fs = fs*(1/4);
t1 = length(dat_low_pass)/fs; % duration in seconds
% Length of the data
len1 = t1*fs; % duration*sampling frequency
% Time vector
timeVec1 =  0:1/fs:(len1-1)/fs;
figure;
plot(timeVec1, dat_low_pass);


%% Step VI
% Calculate the power spectrum of dat_low_pass
% pwelch(dat_low_pass, w, novrlp, [], fs);
[pxx, f] = pwelch(dat_low_pass, [], [], [], fs);
% [pxx, f] = pwelch(dat_low_pass, 4096, 2048, [], fs);
figure;
plot(f,10*log10(pxx));

% Comparing the power spectrum of the data before and after low pass
% filering, we can see that before the data had frequency components above
% 8 kHz, but after applying the filter the frequency components above 2048
% Hz have been discarded, which is expected, as out cut of frequency of the
% butterworth filter was 2048 Hz.



%% Step VII
% Whiten the data using rngmed2
pxx_median_est = rngmed2(pxx, 256);
freq = 0:0.125/(fs*4):1;
bfilt = fir2(500,freq',1./sqrt(pxx_median_est));
dat_whitened = fftfilt(bfilt,dat_low_pass);
% Plot the time series and power spectral density
figure;
subplot(2,1,1);
plot(timeVec1, dat_whitened);
title('Time series of whitened data');
subplot(2,1,2);
pwelch(dat_whitened,4096,2048,[],fs);
% pwelch(dat_whitened,[],[],[],fs);
% fir2 is a finite impulse response filter which is a linear digital filter 
% Characterize dat_whitened


%% Step VIII
% Apply bandpass filter to the whitened data
[b,a]=butter(6,[100 1024]/2048);
dat_whitened=filtfilt(b,a,dat_whitened);
% Plot the bandpassed time series and power spectral density
figure;
subplot(2,1,1);
plot(timeVec1, dat_whitened);
title('Bandpassed whitened time series');
subplot(2,1,2);
pwelch(dat_whitened,[],[],[],fs);

% Why we bandpassed the whitened data - Because LIGO is most sensitive
% between 100 Hz and 1000 Hz. Outside this frequency range signals other
% than our desired signals (noise) prevail, which makes it nearly
% impossible to detect gravitational wave in those frequency ranges


%% Step IX
dat_whitened_resampled = resample(dat_low_pass,1,2);
% New sampling frequency
fs_new = 2048; % fs/2
% Create a power spectrum plot with the newly resampled data
[pxx, f] = pwelch(dat_whitened_resampled,[],[],4*fs_new,fs_new);
figure;
plot(f, 10*log10(pxx));
title('Welch Power Spectral density estimate after the new resampling');


%% Step X
% Trying to be smart 
% Find peaks in the power specrtal density plot
[pks, locs, w, p] = findpeaks(10*log10(pxx));
p_sorted = sort(p, 'descend');
for i = 1:11
    for j = 1:length(locs)
        if p_sorted(i) == pks(j)
            indx(i) = j;
        end
    end
end
freq_list = f(locs(indx));


% Not so smart
% List of peak frequencies from the PSD plot
% f_list = [393, 347.5, 344.75, 343.5, 60, 46.75, 16.75, 13.75, 12, 10.25, 2]; % kHz
% list = [freq_list pks];
% We can see 11 narrowband noise in the final power spectrum. The one at 60
% Hz is the noise due to power lines. The ones that are below 50 Hz are due
% to seismic noise. 


