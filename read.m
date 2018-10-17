%2000 us de periodo de muestreo
%115200 puerto

clear,clc
datalb= dlmread('14_F_LB.log');
data=dlmread('14_F_S.log');
%From 35000 --> Anger     150 seg
% From 45000 --> Amusement 165 seg

%Hoja 1 
%slb=datalb(end-30000:end,1).*5/1024;%Linea base conductividad
%plb=datalb(end-30000:end,2).*5/1024;%Linea base pulso

slb=datalb(1:45000,1).*5/1024;%Linea base conductividad
plb=datalb(1:45000,2).*5/1024;%Linea base pulso
%slb=datalb(:,1).*5/1024;%Linea base conductividad
%plb=datalb(:,2).*5/1024;%Linea base pulso
 
%Sad

%Audition
s=data(:,1).*5/1024; %Conductividad
p=data(:,2).*5/1024; %Pulso
%Anger
%s=data(40001:75000,1).*5/1024; %Conductividad
%p=data(40001:75000,2).*5/1024; %Pulso
% % Amusement
% s=data(47501:82500,1).*5/1024; %Conductividad
% p=data(47501:82500,2).*5/1024; %Pulso


% %6 Primeras pruebas
% s=data(72501:107500,1).*5/1024; %Conductividad
% p=data(72501:107500,2).*5/1024; %Pulso
Fs=500; %F. sampling
T=1/Fs; %T. sampling
%%Prueba eliminacion 500 muestras cada 6 segundos
%%
Tlb=length(slb);    %Numero de muestras de linea base
Tamp=length(p);     %Numero de muestras de tratamiento

%Vector de tiempo  tratamiento
t=(0:Tamp-1)*T; 
t1=t(1001:end); %Vector de tiempo para conductividad
t2=(0:Tamp-1)*T;%Vector de tiempo para la señal de Temperatura
t=t(51:end);    %Vector de tiempo para ritmo

%Vector de tiempo  linea base
Tbase=(0:Tlb-1)*T;
tglb=Tbase(1001:end);%Vector de tiempo para conductividad
tplb=Tbase(51:end); %Vector de tiempo para ritmo

%Definicion ventana hamming 
 
win=hamming(101);   %Tamaño de ventana para Ritmo: 101
win2=hamming(1001); %Tamaño de ventana para Conductividad: 1001

%Filtro pasabanda tipo FIR - Ritmo cardiaco
f1=0.5;
f2=10;

w1=2*f1/Fs;
w2=2*f2/Fs;
B=fir1(100,[w1 w2],'bandpass',win); %Filtro pasabanda tipo FIR
y=filter(B,1,p);
y=y(51:end);        %Eliminacion de primeras 50 muestras
pulsoaf=y;   

%Filtro pasabanda tipo FIR - Ritmo cardiaco linea base
ylb=filter(B,1,plb);
ylb=ylb(51:end);        %Eliminacion de primeras 50 muestras
pulsolb=ylb;   


%Filtrado pasabajos GSR
f3=1;
w3=2*f3/Fs;
B2=fir1(1000,w3,'low',win2);
y1=filter(B2,1,s);
y1=y1(1001:end);      %Eliminacion de primeras 1000 muestras
gsr = y1;
%Filtro de media movil
gsrmm=movmean(gsr,500); %Mejor resultado tomando el valor medio cada 500 muestras (1 segundo)


%Filtrado linea base GSR
y2=filter(B2,1,slb);
y2=y2(1001:end);      %Eliminacion de primeras 1000 muestras
gsrlb = y2;
%Filtro de media movil
gsrmmlb=movmean(gsrlb,500); %Mejor resultado tomando el valor medio cada 500 muestras (1 segundo)

% gsrlb=gsrlb-mean(gsrlb);

%Segmentacion en ventanas de 5 segundos
 
 bx = 2500; %Division del total de los clips en ventanas de 5 segundos  
 na = numel(pulsoaf); %Numero de elementos de la señal de Ritmo cardiaco
 c = mat2cell(pulsoaf,diff([0:bx:na-1,na]));%Division de señal en ventanas de 5000 muestras
 c=c(1:end-1);
 n_iter = length(c);
 
 
%Segmentacion en ventanas de 5 segundos (Linea base)

 blb= 2500; %Division del total de los clips en ventanas de 5 segundos  
 nlb = numel(pulsolb); %Numero de elementos de la señal de Ritmo cardiaco
 clb = mat2cell(pulsolb,diff([0:blb:nlb-1,nlb]));%Division de señal en ventanas de 5000 muestras
 clb=clb(1:end-1);
 n_iterlb = length(clb);
 
%% 
%Definicion ventana hamming para FFT
winh=hamming(2500);

%Calculo de la FFT para la linea base del ritmo cardiaco
 for i=1:n_iterlb
 Xlb=fft(clb{i}.*winh,15000);
 bpmlb=abs(Xlb);
 bpmlb(1:30)=0;                 
 bpmlb(54:end)=0;
 [ylb,in]=max(bpmlb);
 hrlb(i)=(in-1)*60*Fs/(length(Xlb)); % Vector de ritmo cardiaco de linea base en BPM
 end
 
plm=mean(hrlb); %Media de linea base de ritmo cardiaco
pls=std(hrlb);  %Desv standar de linea base de ritmo cardiaco

%Calculo de la FFT para la señal de Ritmo
for i=1:n_iter
 X=fft(c{i}.*winh,15000);
 bpm=abs(X);
 bpm(1:30)=0;            
 bpm(54:end)=0;
[yvalue,k]=max(bpm);
 hr(i)=(k-1)*60*Fs/(length(X)); % Vector de ritmo cardiaco en BPM
end 
%%
%Time domain features for BVP
%Emotion
[amplitud,muestra,ancho,prominence]=findpeaks(pulsoaf,'MinPeakDistance',250,'MinPeakProminence',0.1,'MaxPeakWidth',250);
NN=diff(muestra)/500;
interbeat=60./NN;

%Baseline
[amplitud,muestra,ancho,prominence]=findpeaks(pulsolb,'MinPeakDistance',250,'MinPeakProminence',0.1,'MaxPeakWidth',250);
NNlb=diff(muestra)/500;
interbeatlb=60./NNlb;
%%
%Estimacion de conductividad promedio cada 5 segundos
nk=250;
bk=arrayfun(@(i) mean(gsr(i:i+nk-1)),1:nk:length(gsr)-nk+1);
% siem=1024+2*s+10000/512-s; %Valor en Komhs
%%
%GSR features

diff_gsrmm_lb=diff(gsrmmlb); %Derivate of the baseline GSR
avg_diff_gsrmm_lb=mean(diff_gsrmm_lb) %Average of the GSR's baseline derivate
amountneg_diff_gsrmm_lb=length(diff_gsrmm_lb(diff_gsrmm_lb<0))%Amount of negative samples
ratio_amountneg_lb=amountneg_diff_gsrmm_lb/length(diff_gsrmm_lb)%Proportion of negative values vs the whole lb values
%------------------------------------------------------------------------------------------
diff_gsrmm   =diff(gsrmm);   %Derivate of the GSR during the stimulus
avg_dif_gsrmm = mean(diff_gsrmm) %%Average of the GSR response during the stimulus
amountneg_diff_gsrmm=length(diff_gsrmm(diff_gsrmm<0))%Amount of negative samples
ratio_amountneg=amountneg_diff_gsrmm/length(diff_gsrmm)%Proportion of negative values vs the whole stimulus values
%------------------------------------------------------------------------------------------
%Convolucion GSR con Bartlett 20 puntos
scrb=downsample(gsr,25);
scr=scrb-mean(scrb);
comd=diff(scr);
co=conv(bartlett(20),comd);
%------------------------------------------------------------------------------------------

%Convolucion GSR linea base con Bartlett 20 puntos
scrlb=downsample(gsrlb,25);
scrlb=scrlb-mean(scrlb);
comlb=diff(scrlb);
colb=conv(bartlett(20),comlb);

%%     
%Las señales rojas son las correspondientes a la prueba
figure(1)
subplot(2,2,1)
plot(t1,gsrmm,'r')

subplot(2,2,2)
plot(tglb,gsrmmlb,'b')
 
subplot(2,2,3);
plot(t,pulsoaf,'r')

subplot(2,2,4); 
plot(tplb,pulsolb,'b')

figure(2)
subplot(2,1,1);
plot(gsrmm,'r')
subplot(2,1,2);
plot(diff(gsrmm),'b')

%%

disp('Features extraction')

disp('Baseline data')

%Baseline

%GSR
fprintf('SCR Neutral  Mean %8f .\n',mean(gsrmmlb));
fprintf('SCR Neutral  Standard deviation %8f .\n',std(gsrmmlb));
fprintf('SCR Neutral  Dynamic range %8f .\n',max(gsrmmlb)-min(gsrmmlb));
fprintf('SCR Average Derivate %8f .\n',avg_diff_gsrmm_lb);
fprintf('SCR Amount of negative samples %8f .\n',amountneg_diff_gsrmm_lb);
fprintf('SCR Proportion of negative values vs the whole stimulus ones %8f .\n',ratio_amountneg_lb);

%BVP
fprintf('Heart Rate Neutral Mean %8f .\n',mean(hrlb));
fprintf('Heart Rate Neutral Standard deviation %8f .\n',std(hrlb));
fprintf('Heart Rate Neutral Dynamic range %8f .\n',(max(hrlb)-min(hrlb)));
fprintf('Heart Rate Neutral Mode (Frequency) %8f .\n',mode(hrlb));
fprintf('Heart Rate Mean (Time) %8f .\n',mean(interbeatlb));
fprintf('Heart Rate SDNN (Time) %8f .\n',std(diff(NNlb)));
fprintf('Heart Rate RMSSD (Time) %8f .\n',rms(diff(NNlb)));

%Emotion

disp('Emotion')
fprintf('SCR Mean %8f .\n',mean(gsrmm));
fprintf('SCR Standard deviation %8f .\n',std(gsrmm));
fprintf('SCR Neutral  Dynamic range %8f .\n',(max(gsrmm)-min(gsrmm)));
fprintf('SCR Average Derivate %8f .\n',avg_dif_gsrmm);
fprintf('SCR Amount of negative samples %8f .\n',amountneg_diff_gsrmm);
fprintf('SCR Proportion of negative values vs the whole stimulus values %8f .\n',ratio_amountneg);


fprintf('Heart Rate Mean (FFT) %8f .\n',mean(hr));
fprintf('Heart Rate Standard deviation %8f .\n',std(hr));
fprintf('Heart Rate Dynamic range %8f .\n',(max(hr)-min(hr)));
fprintf('Heart Rate Mode (Frequency) %8f .\n',mode(hr));
fprintf('Heart Rate Mean (Time) %8f .\n',mean(interbeat));
fprintf('Heart Rate SDNN (Time) %8f .\n',std(diff(NN)));
fprintf('Heart Rate RMSSD (Time) %8f .\n',rms(diff(NN)));


%Array info

infolb=[mean(gsrmmlb),std(gsrmmlb),max(gsrmmlb)-min(gsrmmlb),avg_diff_gsrmm_lb,amountneg_diff_gsrmm_lb,ratio_amountneg_lb,mean(hrlb),std(hrlb),max(hrlb)-min(hrlb),mode(hrlb),mean(interbeatlb),std(diff(NNlb)),rms(diff(NNlb))];
infoem=[mean(gsrmm),std(gsrmm),max(gsrmm)-min(gsrmm),avg_dif_gsrmm,amountneg_diff_gsrmm,ratio_amountneg,mean(hr),std(hr),max(hr)-min(hr),mode(hr),mean(interbeat),std(diff(NN)),rms(diff(NN))];

fprintf('%8f,',infolb)
disp('.')
fprintf('%8f,',infoem)


