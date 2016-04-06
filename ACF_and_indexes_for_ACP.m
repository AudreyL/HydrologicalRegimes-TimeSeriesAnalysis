
%# Audrey LUSTIG
%# Last update: 28_08_2012
%# 
%# calcule et trace l'ACF d'une série temporelle
%# Calcul des différents indices à ajouter dans l'ACP

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# Calcul de l'ACF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%# Calcule l'ACF du lag 1 ou lag length(S)-1
%# En theorie on peut s'arreter a 2/3 de la longueur
nLag=floor(length(S)*2/3);
[ACF,Lags,Bounds]=autocorr(S,nLag);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# Intervalle de confiance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# Etape 1
%# Genere nb_simul bruit rouge avec AR1 (ici foac) equivalent a la serie etudiee
foac=ACF(2);
nb_simul=100;
for i=1:nb_simul
    %# Change la seed
    a=10*rand(1,1);                         
    randn('seed',floor(sum(10*clock)/a));
    %# genere le bruit rouge
    G(1,1) = randn(1,1);
    for k=1:N-1
        G(1,k+1) = (foac*G(1,k)) + ((1-foac)*randn(1,1));
    end
    G = G - repmat(mean(G),1,1);
    G = G/std(G);
    for j=1:N-1
        G(1,j) = G(1,j) .* sqrt(nanstd(S)^2);
    end
    %# stockage des bruits rouges
    S_tab(i,1:N)=G;
    clear G a clock k j
end

%# calcul l'ACF pour chacun des bruits rouges du lag 1 ou lag length(S)-1
ACF_ic=[];
for i=1:nb_simul
    ACF_red_noise=autocorr(S_tab(i,:),nLag);
    ACF_ic=[ACF_ic;ACF_red_noise];
    clear ACF_red_noise
end

%# calcul l'intervalle de confiance à 99%
ACF_up=quantile(ACF_ic,0.99);
ACF_down=quantile(ACF_ic,0.01);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# Figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# Plot log de la serie et ACF avec intervalle de confiance
%# WARNING: pour les anomalies ne pas tracer les débits en log


%# Nom de la figure

 
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
plot(Lags(2:end),ACF(2:end))
set(gca,'Xlim',[Lags(2),Lags(end)])
Xlabel('Jours')
Ylabel('ACF')
hold on
plot(Lags(2:end),ACF_up(2:end),'r-')
plot(Lags(2:end),ACF_down(2:end),'r-')
hold off
set(gca,'Ylim',[min(min(ACF(2:end)),min(ACF_down(2:end))),max(max(ACF(2:end)),max(ACF_up(2:end)))])

%# Pour sauvegarder l'image dans un dossier
% name_f=['D:\Inde\BVH\Results\ACF\ACF_Jeu_Continu\',char(name{1}(1)),'_',char(name{1}(2)),'_', char(name{1}(3)),'.jpg'];
% saveas(gcf,name_f)
%     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# Statistiques quantitatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# recherche valeur au dessus et en dessous de l'intervalle de confiance
positive=[];
negative=[];
for i=1:length(ACF)
	if ACF(i)>ACF_up(i)
    	positive=[positive, (ACF(i)-ACF_up(i))]; %# valeurs significatives positives
    end
    if ACF(i)<ACF_down(i)
        negative=[negative, -(ACF(i)-ACF_down(i))];%# valeurs significatives negatives
    end
end
%# Ensemble des valeurs significatives
siginifcatif=[negative,positive];

% Calcul pour le spectre ACF
if length(siginifcatif)>=1
    E_siginifcatif_moy=mean(siginifcatif);
    if length(siginifcatif)>=2 	
    	E_siginifcatif_std=std(siginifcatif);
    else
   	E_siginifcatif_std=0;
    end
    E_siginifcatif_length=length(siginifcatif)/length(ACF);
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













