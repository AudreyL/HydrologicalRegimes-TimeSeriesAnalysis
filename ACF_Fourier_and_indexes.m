%# Afficher transformée de Fourier et ACF 
%# Calcul les indices caracteristiques pour l'ACP
%# Des données brutes et des anoamlies
%# Audrey LUSTIG
%# Last update : 29-08-2012

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# Parametre utilisateur
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Repertoire= 'D:\BVH\Haiti_data_streamflows\Final\Data\Raw_data'; % repertoire des données brutes
Repertoire_ano= 'D:\BVH\Haiti_data_streamflows\Final\Data\Anomalies'; % repertoire des anomalies

dirName=Repertoire; %# directory name
dirData = dir(dirName); %# Get the data for the current directory
dirIndex = [dirData.isdir]; %# Find the index for directories
fileList = {dirData(~dirIndex).name}'; %# Get a list of the files

dirName_ano=Repertoire_ano; %# directory name
dirData_ano = dir(dirName_ano); %# Get the data for the current directory
dirIndex_ano = [dirData_ano.isdir]; %# Find the index for directories
fileList_ano = {dirData_ano(~dirIndex_ano).name}'; %# Get a list of the files

% mat={'id','station','riviere','ACF_moy','ACF_std','ACF_l','fft_moy','fft_std','fft_l'};
% mat_ano={'id','station','riviere','ACF_moy','ACF_std','ACF_l','fft_moy','fft_std','fft_l'};
for f=1:length(fileList)
    %# Lecture des debits
    fileName = strcat(dirName,'\',char(fileList(f))); %# get directory
    fid = fopen(fileName,'r');  %# Open the file
    data = textscan(fid,'%s %f','treatAsEmpty',{'NA','na'}); %# read data as a list
    fclose(fid); %# close the file
    Debits = cell2mat(data(2)); %# transform the list in matric   
    S=transpose(Debits); %# Debits est en colonne, S en ligne
    N=length(S); %# Longueur de la serie
    clear fid

    %# Lecture des anoamlies
    fileName_ano = strcat(dirName_ano,'\',char(fileList_ano(f))); %# get directory
    fid_ano = fopen(fileName_ano,'r');  %# Open the file
    data_ano = textscan(fid_ano,'%s %f','treatAsEmpty',{'NA','na'}); %# read data as a list
    fclose(fid_ano); %# close the file
    Ano = cell2mat(data_ano(2)); %# transform the list in matric   
    S_ano=transpose(Ano); %# Debits est en colonne, S en ligne
    clear fid_ano
    
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
    time = [day1:day1+N-1]*1/365 + str2double(TOKEN1) ;
    clear TOKEN1 REMAIN1 TOKEN2 REMAIN2 day1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %# Calcul de l'ACF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %# Calcule l'ACF du lag 1 ou lag length(S)*2/3
    nLag=floor(N*2/3);
    [ACF,Lags]=autocorr(S,nLag); % donnees brutes
    [ACF_ano]=autocorr(S_ano,nLag); % anoamlie
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %# Spectre de Fourier (FFT)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [W,E]=ezfft(time,S); % donnees brutes
    [W_ano,E_ano]=ezfft(time,S_ano); % anomalies
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %# Intervalle de confiance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %# Genere nb_simul bruit rouge de même autocorrelation de premier ordre que le signal etudie
	nb_simul=100;
    foac=ACF(2);
    S_tab=[];
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
        S_tab=[S_tab;G];
        clear G a k j clock
    end
    
    % bruit rouge pour les anomalies
    foac_ano=ACF_ano(2);
    S_tab_ano=[];
    for i=1:nb_simul
        a=10*rand(1,1);                         
        randn('seed',floor(sum(10*clock)/a));
        G(1,1) = randn(1,1);
        for k=1:N-1
            G(1,k+1) = (foac_ano*G(1,k)) + ((1-foac_ano)*randn(1,1));
        end
        G = G - repmat(mean(G),1,1);
        G = G/std(G);
        for j=1:N-1
            G(1,j) = G(1,j) .* sqrt(nanstd(S_ano)^2);
        end
        S_tab_ano=[S_tab_ano;G];
        clear G a k j clock
    end
    
    % IC ACF donnees brutes
    %# calcul l'ACF pour chacun des bruits rouges du lag 1 ou lag length(S)-1
    ACF_ic=[];
    for i=1:nb_simul
        ACF_red_noise=autocorr(S_tab(i,:),nLag);
        ACF_ic=[ACF_ic;ACF_red_noise];
        clear ACF_red_noise
    end
    ACF_up=quantile(ACF_ic,0.99);
    ACF_down=quantile(ACF_ic,0.01);
    
    % IC ACF anoamlies
    ACF_ic_ano=[];
    for i=1:nb_simul
        ACF_red_noise=autocorr(S_tab_ano(i,:),nLag);
        ACF_ic_ano=[ACF_ic_ano;ACF_red_noise];
        clear ACF_red_noise
    end
    ACF_up_ano=quantile(ACF_ic_ano,0.99);
    ACF_down_ano=quantile(ACF_ic_ano,0.01);
    
    % IC Fourier donnees brutes
    cf_fft=[];
    for i=1:nb_simul
        [Wi,Ei] = ezfft(time,S_tab(i,:));
        cf_fft=[cf_fft;Ei];
        clear Ei Wi
    end
    fft_up=quantile(cf_fft,0.99);
    fft_down=quantile(cf_fft,0.01);
    
    % IC Fourier anoamlies
    
    cf_fft_ano=[];
    for i=1:nb_simul
        [Wi_ano,Ei_ano] = ezfft(time,S_tab_ano(i,:));
        cf_fft_ano=[cf_fft_ano;Ei_ano];
        clear Ei_ano Wi_ano
    end
    fft_up_ano=quantile(cf_fft_ano,0.99);
    fft_down_ano=quantile(cf_fft_ano,0.01);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %# Figure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % A decommenter pour affichier
    % Une colonne de droite avec données brutes, ACF et Fourier associés
    % Une colonne de gauche avec anoamlies ACF et Fourier associés
    name=fileList(f)
    name=regexp(name,'_','split');
    name_plot=['Station: ',char(name{1}(2)),', riviere: ', char(name{1}(3)),', id: ',char(name{1}(1))];
    figure(1)
    clf
    subplot(3,2,1)
    plot(time,S)
    set(gca,'Xlim',[time(1),time(end)])
    set(gca,'Ylim',[min(S),max(S)])
    Xlabel('Jours')
    Ylabel('Debits (m^3/j)')
    title(name_plot)
    
    subplot(3,2,2)
    plot(time,S_ano)
    set(gca,'Xlim',[time(1),time(end)])
    set(gca,'Ylim',[min(S_ano),max(S_ano)])
    Xlabel('Jours')
    Ylabel('Anoamlies (m^3/j)')
    title('Anomalies')
    
    subplot(3,2,3)
    plot(Lags(2:end),ACF(2:end))
    set(gca,'Xlim',[Lags(2),Lags(end)])
    Xlabel('Jours')
    Ylabel('ACF')
    title('Autocorrelation')
    hold on
    plot(Lags(2:end),ACF_up(2:end),'r-')
    plot(Lags(2:end),ACF_down(2:end),'r-')
    hold off
    set(gca,'Ylim',[min(min(ACF(2:end)),min(ACF_down(2:end))),max(max(ACF(2:end)),max(ACF_up(2:end)))])

    subplot(3,2,4)
    plot(Lags(2:end),ACF_ano(2:end))
    set(gca,'Xlim',[Lags(2),Lags(end)])
    Xlabel('Jours')
    Ylabel('ACF anomalies')
    title('Autocorrelation')
    hold on
    plot(Lags(2:end),ACF_up_ano(2:end),'r-')
    plot(Lags(2:end),ACF_down_ano(2:end),'r-')
    hold off
    set(gca,'Ylim',[min(min(ACF_ano(2:end)),min(ACF_down_ano(2:end))),max(max(ACF_ano(2:end)),max(ACF_up_ano(2:end)))])

    subplot(3,2,5)
    ezfft(time,S,'freq')
    title('Transformee de Fourier Rapide')
    hold on
    loglog(W/(2*pi),2*pi*fft_up,'r-')
    loglog(W/(2*pi),2*pi*fft_down,'r-')
    hold off
    set(gca,'Xlim',[W(2)/(2*pi),W(end)/(2*pi)])
    
    subplot(3,2,6)
    ezfft(time,S_ano,'freq')
    title('Transformee de Fourier Rapide')
    hold on
    loglog(W_ano/(2*pi),2*pi*fft_up_ano,'r-')
    loglog(W_ano/(2*pi),2*pi*fft_down_ano,'r-')
    hold off
    set(gca,'Xlim',[W_ano(2)/(2*pi),W_ano(end)/(2*pi)])
    
    %# Pour sauvegarder l'image dans un dossier
    name_f=['D:\Inde\BVH\Results\ACF\Jeu Continu\',char(name{1}(1)),'_',char(name{1}(2)),'_', char(name{1}(3)),'.jpg'];
	saveas(gcf,name_f)    
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %# Statistiques quantitatives
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Recherche des élements significatifs pour l'ACF
%     ACFp=[];
%     ACFn=[];
%     for i=1:length(ACF)
%         if ACF(i)>ACF_up(i)
%             ACFp=[ACFp, (ACF(i)-ACF_up(i))]; %# valeurs significatives positives
%         end
%         if ACF(i)<ACF_down(i)
%             ACFn=[ACFn, -(ACF(i)-ACF_down(i))];%# valeurs significatives negatives
%         end
%     end
%     ACFs=[ACFp,ACFn];
%     
%     % Recherche des élements significatifs pour l'ACF ano 
%     ACF_ano_p=[];
%     ACF_ano_n=[];
%     for i=1:length(ACF_ano)
%         if ACF_ano(i)>ACF_up_ano(i)
%             ACF_ano_p=[ACF_ano_p, (ACF_ano(i)-ACF_up_ano(i))]; %# valeurs significatives positives
%         end
%         if ACF_ano(i)<ACF_down_ano(i)
%             ACF_ano_n=[ACF_ano_n, -(ACF_ano(i)-ACF_down_ano(i))];%# valeurs significatives negatives
%         end
%     end
%     ACF_ano_s=[ACF_ano_p,ACF_ano_n];
%     
%     % Recherche des élements significatifs pour la fft
%     Ep=[];
%     En=[];
%     for i=1:length(E)
%         if E(i)>fft_up(i)
%             Ep=[Ep, (E(i)-fft_up(i))]; %# valeurs significatives positives
%         end
%         if E(i)<fft_down(i)
%             En=[En, -(E(i)-fft_down(i))];%# valeurs significatives negatives
%         end
%     end
%     Es=[Ep,En];
%     
%     % Recherche des élements significatifs pour la fft ano
%     E_ano_p=[];
%     E_ano_n=[];
%     for i=1:length(E_ano)
%         if E_ano(i)>fft_up_ano(i)
%             E_ano_p=[E_ano_p, (E_ano(i)-fft_up_ano(i))]; %# valeurs significatives positives
%         end
%         if E_ano(i)<fft_down_ano(i)
%             E_ano_n=[E_ano_n, -(E_ano(i)-fft_down_ano(i))];%# valeurs significatives negatives
%         end
%     end
%     E_ano_s=[E_ano_p,E_ano_n];
%     
%     % Statistiques quantitatives de l'ACF
%     if length(ACFs)>=1
%         ACF_moy=mean(ACFs); 
%         if length(ACFs)>=2
%             ACF_std=std(ACFs);
%         else
%             ACF_std=0;
%         end
%         ACF_length=length(ACFs)/length(ACF);
%     else
%         ACF_moy=0; 
%         ACF_std=0;
%         ACF_length=0;
%     end 
%     
%     % Statistiques quantitatives de l'ACF ano
%     if length(ACF_ano_s)>=1
%         ACF_moy_ano=mean(ACF_ano_s); 
%         if length(ACF_ano_s)>=2
%             ACF_std_ano=std(ACF_ano_s);
%         else
%             ACF_std_ano=0;
%         end
%             ACF_length_ano=length(ACF_ano_s)/length(ACF_ano);
%    	else
%         ACF_moy_ano=0; 
%         ACF_std_ano=0;
%         ACF_length_ano=0;
%     end
%     
%     % Statistiques quantitatives de la fft
%     if length(Es)>=1
%         E_moy=mean(Es); 
%         if length(Es)>=2
%             E_std=std(Es);
%         else
%             E_std=0;
%         end
%         E_length=length(Es)/length(E);
%     else
%         E_moy=0; 
%         E_std=0;
%         E_length=0;
%     end 
%     
%     % Statistiques quantitatives de la fft continu
%     if length(E_ano_s)>=1
%         E_moy_ano=mean(E_ano_s); 
%         if length(E_ano_s)>=2
%             E_std_ano=std(E_ano_s);
%         else
%             E_std_ano=0;
%         end
%         E_length_ano=length(E_ano_s)/length(E_ano);
%     else
%         E_moy_ano=0; 
%         E_std_ano=0;
%         E_length_ano=0;
%     end 
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %# Stockage des informations quantitiaves 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     m={char(name{1}(1)),char(name{1}(2)),char(name{1}(3)),num2str(ACF_moy),num2str(ACF_std),num2str(ACF_length),num2str(E_moy),num2str(E_std),num2str(E_length)};
%     m_ano={char(name{1}(1)),char(name{1}(2)),char(name{1}(3)),num2str(ACF_moy_ano),num2str(ACF_std_ano),num2str(ACF_length_ano),num2str(E_moy_ano),num2str(E_std_ano),num2str(E_length_ano)};
%     
%     mat=[mat;m];
%     mat_ano=[mat_ano;m_ano];
    
    clear fileName data Debits S N 
    clear fileName_ano data_ano Ano S_ano
    clear nLag ACG Lags ACF_ano  
    clear W E W_ano E_ano
    clear foac foac_ano nb_simul S_tab S_tab_ano
    clear ACF_ic ACF_up ACF_down ACF_ic_ano ACF_up_ano ACF_down_ano
    clear cf_fft fft_up fft_down cf_fft_ano fft_up_ano fft_down_ano
    clear name name_plot
    clear ACFp ACFn ACFs ACF_ano_p ACF_ano_n ACF_ano_s
    clear Ep En Es E_ano_p E_ano_n E_ano_s
    clear ACF_moy ACF_std ACF_length ACF_moy_ano ACF_std_ano ACF_length_ano
    clear E_moy E_std E_length E_moy_ano E_std_ano E_length_ano
    clear m m_ano
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# Ecriture dans un fichier
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fid=fopen('D:\Inde\BVH\Results\ACF\Jeu Commun 2529\ACF_Fourier_quantitatif.txt','wt'); %# Nom du fichier
% [rows,cols]=size(mat);
% for i=1:rows
%     fprintf(fid,'%s\t',mat{i,1:end-1});
%     fprintf(fid,'%s\n',mat{i,end});
% end
% fclose(fid);
% 
% fid=fopen('D:\Inde\BVH\Results\ACF\Jeu Commun 2529\ACF_Fourier_ano_quantitatif.txt','wt'); %# Nom du fichier
% [rows,cols]=size(mat_ano);
% for i=1:rows
%     fprintf(fid,'%s\t',mat_ano{i,1:end-1});
%     fprintf(fid,'%s\n',mat_ano{i,end});
% end
% fclose(fid);

