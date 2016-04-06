%# Audrey LUSTIG
%# Last update: 28_08_2012
%# 
%# Transformée Fourier rapide d'une série temporelle
%# Determine et affiche le spectre de Fourier
%# Calcul des différents indices à ajouter dans
clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# Lecture des donnees
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%# Ouverture du fichier
fileName = 'D:\BVH\Haiti_data_streamflows\Final\Data\Raw_data\21_LA GORGE_BLANCHE_170.08.txt'
fid = fopen(fileName,'r');  %# Open the file
data = textscan(fid,'%s %f','treatAsEmpty',{'NA','na'}); %# read data as a list
fclose(fid); %# close the file

%# Lecture donnees
Debits = cell2mat(data(2)); %# transform the list in matric   
S=transpose(Debits); %# Debits est en colonne, S en ligne
N=length(S); %# Longueur de la serie

%# Creation d'un vecteur temps (recherche premier jour, mois, annee)
[TOKEN1,REMAIN1]=strtok(data{1}(1),'-');
[TOKEN2,REMAIN2]=strtok(REMAIN1,'-');
if str2double(TOKEN2)==1
    day1=str2double(strtok(REMAIN2,'-'));
end
if str2double(TOKEN2)==2
    day1=str2double(strtok(REMAIN2,'-'))+31;
end
if str2double(TOKEN2)==3
    day1=str2double(strtok(REMAIN2,'-'))+59;
end
if str2double(TOKEN2)==4
    day1=str2double(strtok(REMAIN2,'-'))+90;
end
if str2double(TOKEN2)==5
    day1=str2double(strtok(REMAIN2,'-'))+120;
end
if str2double(TOKEN2)==6
    day1=str2double(strtok(REMAIN2,'-'))+151;
end
if str2double(TOKEN2)==7
    day1=str2double(strtok(REMAIN2,'-'))+181;
end
if str2double(TOKEN2)==8
    day1=str2double(strtok(REMAIN2,'-'))+212;
end
if str2double(TOKEN2)==9
    day1=str2double(strtok(REMAIN2,'-'))+243;
end
if str2double(TOKEN2)==10
    day1=str2double(strtok(REMAIN2,'-'))+273;
end
if str2double(TOKEN2)==11
    day1str2double(strtok(REMAIN2,'-'))+304;
end
if str2double(TOKEN2)==12
    day1=str2double(strtok(REMAIN2,'-'))+334;
end
time = [day1:day1+length(Debits)-1]*1/365 + str2double(TOKEN1) ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# Spectre de Fourier (FFT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[W,E]=ezfft(time,S);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# Intervalle de confiance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# Genere nb_simul bruit rouge de même autocorrelation de premier ordre que le signal etudie
[ACF,Lags,Bounds]=autocorr(S,3);
foac=ACF(2);
nb_simul=5;
for i=1:nb_simul
	a=10*rand(1,1);                         
	randn('seed',floor(sum(10*clock)/a));
	G(1,1) = randn(1,1);
    for k=1:N-1
    	G(1,k+1) = (foac*G(1,k)) + ((1-foac)*randn(1,1));
    end
    G = G - repmat(mean(G),1,1);
    G = G/std(G);
    for j=1:N-1
    	G(1,j) = G(1,j) .* sqrt(nanstd(S)^2);
    end
    S_tab(i,1:N)=G;
    clear G a clock j k
end

%# Calcul le spectre de Fourier pour chacun de ces bruits rouges
cf_fft=[];
for i=1:nb_simul
	[Wi,Ei] = ezfft(time,S_tab(i,:));
	cf_fft=[cf_fft;Ei];
    clear Ei Wi
end

%# Determine l'intervalle de confiance a 99%
fft_up=quantile(cf_fft,0.99);
fft_down=quantile(cf_fft,0.01);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# Figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# Affiche le spectre de Fourier et son intervalle de confiance
figure(1);
clf;
subplot(2,1,1)
plot(time,log(S))
set(gca,'Xlim',[time(1),time(end)])
set(gca,'Ylim',[min(log(S)),max(log(S))])
Xlabel('Jours')
Ylabel('log(Debits) (m^3/j)')
title('Chroniques')

subplot(2,1,2)
ezfft(time,S)
hold on
loglog(W/(2*pi),2*pi*fft_up,'r-')
loglog(W/(2*pi),2*pi*fft_down,'r-')
hold off
set(gca,'Xlim',[W(2)/(2*pi),W(end)/(2*pi)])
title('Spectre de Fourier')
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# Statistiques quantitatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# recherche valeur au dessus et en dessous de l'intervalle de confiance
positive=[];
negative=[];
for i=1:length(E)
	if E(i)>fft_up(i)
    	positive=[positive, (E(i)-fft_up(i))]; %# valeurs significatives positives
    end
    if E(i)<fft_down(i)
        negative=[negative, -(E(i)-fft_down(i))];%# valeurs significatives negatives
    end
end
%# Ensemble des valeurs significatives
siginifcatif=[negative,positive];


%# Calcul des indices caracteristiques pour le spectre de Fourier    
 if length(siginifcatif)>=1
    E_siginifcatif_moy=mean(siginifcatif); 
    if length(siginifcatif)>=2
    	E_siginifcatif_std=std(siginifcatif);
    else
    	E_siginifcatif_std=0;
    end
    E_siginifcatif_length=length(siginifcatif)/length(E);
else
    E_siginifcatif_moy=0; 
    E_siginifcatif_std=0;
    E_siginifcatif_length=0;
 end 
disp('Moyenne du spectre significatif')
disp(E_siginifcatif_moy)
disp('Ecart-type du spectre significatif')
disp(E_siginifcatif_std)
disp('Longueur du spectre significatif')
disp(E_siginifcatif_length) 


    
