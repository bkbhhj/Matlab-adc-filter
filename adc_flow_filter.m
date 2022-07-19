% Try to figure out if we can see the doppler offset with an adaptive
% filter that cancels the TX signal and leaves only the RX signal.

% Doppler equations
% freq_observer = (vel_sound+vel_observer)/(vel_sound-vel_source)*freq_actual
% vel_source = vel_sound*(freq_observed-freq_actual)/freq_observed
% 1Hz offset => ((640001 - 640000)/640001)*1484

close all;
clear all;
min_vel = 0.03;         % Desired resolution of the velocity (m/sec)
N = 2*ceil(3/min_vel);  % Works out to be about 200+ bins.

Fs = 6.4e6;             % Sampling rate.
TX_Freq = 639999.98569488525390625;

% Parameters for demodulation: TX_Freq_demod will move the spectrum from
% around TX to around 5kHz, new sampling rate will be 40kHz (160x
% compression)
Fs_t = 4.0e4;
TX_Freq_demod = TX_Freq - 5e3;

% Filter design for the bandpass that provides fairly close notches to
% minimize the residual of some of the spurious clocks in the data -
% demodulated.
%later:f = [0 1e3/(Fs_t/2) 3e3/(Fs_t/2) 7e3/(Fs_t/2) 9e3/(Fs_t/2) 1];
%later:a = [0 0.0 1.0 1.0 0.0 0];
%later:b = firpm(101,f,a)';


%create empty figure for the final subplot
figure(1)

% configuration
% Stepsize configuration
stepsize1 = 7;
stepsize2 = 11;
stepsize_threshold = 4096;

% Threshold for convergence
conv_threshold = 4096;

% Other config
exclude_saturated = 1;    % Exclude saturated blocks from test
use_cached = 1;           % Use cached results from disk to speed up test
make_threshold_plots = 0; % Create threshold plots for fading
display_stats = 0;        % Display fading stats

%create empty figure for the final subplot
figure(1)
%01_16_2019_13_24_19.bin - some results

fileID = fopen('03_25_2020_20_44_27.bin');%02_22_2019_13_15_20.bin'); %01_17_2019_12_12_21.bin');%01_17_2019_12_03_17.bin');%01_17_2019_11_10_13.bin');%01_16_2019_13_24_34.bin'); %01_16_2019_13_24_15.bin'); %01_16_2019_13_24_15.bin');%01_16_2019_13_24_48.bin'); %01_16_2019_13_25_11
ecg_sygnal = fread(fileID, [2,Inf], 'uint16'); %, 2

%fileID2 = fopen('a.bin','w+');
%fwrite(fileID2,ecg_sygnal);
%fclose(fileID2);

Fd = 50000000; %Частота дискретизации
dt = 1/Fd; % период дискретизации
%t = 0:dt:(length(ecg_sygnal)-1)*dt; % массив времени

%first_sygnal = ecg_sygnal(2,1:end); % is Ok for 1st channel
first_sygnal = ecg_sygnal(2,:);%1:end);
second_sygnal = ecg_sygnal(1,:);%1:end); % is Ok for 1st channel
%first_sygnal = ecg_sygnal(1,1:end);% is  for 2nd channel

%ecg_sygnal(:,1);
%

N=length(first_sygnal); %ecg_sygnal);
t = [0:N-1]/Fd;

w = 50/(250/2);
bw = w;
[num, den] = iirnotch(w,bw);
ecg_notch = filter(num,den,first_sygnal);%ecg_sygnal);
ecg_notch2 = filter(num,den,second_sygnal);%ecg_sygnal);

N1=length(ecg_notch);
t1 = [0:N1-1]/Fd;


%plot(t(1:100), ecg_sygnal(1:100)) % первые 100 точек
%xlabel('t, c')
%fclose all
figure 
%plot(t, ecg_sygnal);
hold on

%a= fir1(100, [0.045 0.056],'stop');
%fir_ecg = filter(a,1,ecg_sygnal);


grid on
grid minor 

windowSize = 1;
b = (1/windowSize ) * ones(1, windowSize);
avg = filter(b, 1, ecg_notch);
avg2 = filter(b, 1, ecg_notch2);

%fir_ecg = filter(a,1,avg);
%plot(t1, ecg_notch);%fir_ecg);

plot(t, avg, t, avg2);

hold off

%plot(A);
