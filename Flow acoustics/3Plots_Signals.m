%% Werte
nfft=300000;
fs=5128;
WINDOW=ones(1,300000);
%% Channel 1 - Piezo
X1=sampling1.VarName2;
[psd4,f4]=pwelch(X1,WINDOW,0,nfft,fs)
plot(f4,psd4)
title('Piezo')
%% Channel 2 - Camera
X2=sampling1.VarName3;
pwelch(X2,WINDOW,0,nfft,fs)
title('Camera')
%% Channel 3 - Diode
X3=sampling1.VarName4;
pwelch(X3,WINDOW,0,nfft,fs)
title('Diode')