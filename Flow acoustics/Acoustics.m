%%
close all
clear all
clc
%% specifiy the Directory
workDir = pwd;
filename = 'Labview.txt';
LabviewRec = fullfile(workDir, filename);                           %LabviewRec contains the full Directory to the textfile

%% Import Data
delimiterIn = ',';                                                          %Column separator character
headerlinesIn = 23;                                                         %reading numeric data starting from line headerlinesIn+1
meas = importdata(LabviewRec, delimiterIn, headerlinesIn);                  

%% Seperate Columns
time = meas.data (:, 1);                                                    
piezo = meas.data(:, 2);                                                    
clock = meas.data(:, 3);                                                    
photodiode = meas.data(:, 4);                                               

%% sampling Frequency
f_s = 1/(time(3)-time(2));                                                  %sampling Frequency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
index=zeros(length(clock),1);
for i=1:length(clock)
index(i)=((clock(i))>2);
end

index1 = find(index);
c1 = index1(1); 

for i=2:length(index1)
if (index1(i)-index1(i-1))>1
    c2 = index1(i);
    break
end
end

f_cs = 1/(time(c2)-time(c1));                                               %sampling Frequency of camera

%%

WINDOW = ones(1, length(piezo));
NFFT = length(piezo);

[psd, f] = pwelch(piezo-mean(piezo), WINDOW,0,NFFT ,f_s);


figure(1);
semilogy(f,psd)
xlim([0 500])
ylim([10^-10 10^5])
title('Piezo Sensor Frequency')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
set(gca,'FontSize',18);

%%
WINDOW = ones(1, length(photodiode));
NFFT = length(photodiode);

[psd, f] = pwelch(photodiode-mean(photodiode), WINDOW,0,NFFT ,f_s);


figure(2);
semilogy(f,psd)
xlim([0 500])
ylim([10^-15 10])
title('Photodiode Frequency')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
set(gca,'FontSize',18);
%%
WINDOW = ones(1, length(clock));
NFFT = length(clock);

[psd, f] = pwelch(clock-mean(clock), WINDOW,0,NFFT ,f_s);

figure(3);
semilogy(f,psd)
xlim([0 500])
ylim([10^-10 10^5])
title('Camera frequency')
xlabel('Time [s]')
ylabel('Frequency [Hz]')

set(gca,'FontSize',18);

%% Plot of the Power Spectral Density

nfft=100;
fs=31.26;
WINDOW=ones(1,350);

figure();

for i=1:5
    i
    [psd,f]=pwelch((eta(i,:)-mean(eta(i,:))),WINDOW,0,nfft,fs);
    plot(f,abs(psd)+(i-1)*100)
    hold on;
end

title('PSD modes 1 to 5')
xlabel('Frequency [Hz]')
ylabel('Power Amplitude')
legend('mode 1', 'mode 2', 'mode 3', 'mode 4', 'mode 5')

figure(1)
semilogy(f, psd);
xlim ([0 500])
xlabel ('Frequency [Hz]')
ylabel ('Power Amplitude')
title ('Power Spectral Density')

[a,b]= max(psd)
f_b = f(b);                                                                 %frequency of bottle

set(gca,'FontSize',18);
%% Plot piezo PSD

WINDOW = ones(1, length(piezo));
NFFT = length(piezo);

[psd, f] = pwelch(piezo-mean(piezo), WINDOW,0,NFFT ,f_s);

figure(1)
semilogy(f, psd);
xlim ([0 500])
xlabel ('Frequency [Hz]','FontSize',16)
ylabel ('Power Amplitude','FontSize',16)
title ('Power Spectral Density of Piezo 2','FontSize',16)

[a,b]= max(psd);
f_b = f(b);                                                                 %frequency of bottle

%% Plot photodiode PSD

WINDOW = ones(1, length(photodiode));
NFFT = length(photodiode);

[psd, f] = pwelch(photodiode-mean(photodiode), WINDOW,0,NFFT ,f_s);

figure(2)
semilogy(f, psd);
xlim ([0 500])
xlabel ('Frequency [Hz]','FontSize',16)
ylabel ('Power Amplitude','FontSize',16)
title ('Power Spectral Density of Photodiode 2','FontSize',16)
