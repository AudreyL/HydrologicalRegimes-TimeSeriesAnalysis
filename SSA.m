% SSA Singular Spectrum Analysis of two hydorlogic Time series (Original and Anomalies)
%  Perform SSA (Allen and Smith, 1996) and develop confidence intervals
%  for SSA selection using Monte Carlo simulations
%
% Inputs:
%  repertoire1 and  = file name with two colonns. First on are date DD/MM/YYYY,
%  Second one is the studied streamflow
%  repertoire2 = file name with two colonns. First on are date DD/MM/YYYY,
%  Second one is the studied streamflow anoamlies
%  M = embedding dimension (1 x 1) (lag covariance)
%  nb_simul = number of MC simulation (1 x 1)
%
% Outputs:
%  A = Principal Components for the original series and anomalies
%  evec = eigenvectors for the original series and anomalies
%  eval = eigenvalues for the original series and anomalies
%  varexp = fraction variance explained for each eigenvalue for the original series and anomalies
%  figure = structure containing the confidence intervals from 
%          white and/or red noise Monte Carlo simulation for the original series and anomalies  
%
%
% References:
% Allen M., Smith L.A., 1996: Monte Carlo SSA: Detecting irregular 
% oscillations in the presence of coloured noise,  J. Clim.,  9 , 3373-3404. 
%
% Ghil M., R. M. Allen, M. D. Dettinger, K. Ide, D. Kondrashov, M. E. Mann, 
% A. Robertson, A. Saunders, Y. Tian, F. Varadi, and P. Yiou, 2002: 
% Advanced spectral methods for climatic time series, Rev. Geophys.
% 40(1), pp. 3.1-3.41, 10.1029/2000GR000092.
%
% Tsonis, A.A., and J.B. Elsner, 1992: Oscillating global temperature, Nature, 356, 751.
%
% Written by Audrey LUSTIG (audreylustig@gmail.com)
% Inspired from Kevin Anchukaitis (kja@ldeo.columbia.edu) 
% Last updated: 30-08-2012


clear all
close all
clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# Lecture des données
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Repertoire= 'D:\BVH\Haiti_data_streamflows\Final\Data\Raw_data'; % repertoire des données brutes
Repertoire_ano= 'D:\BVH\Haiti_data_streamflows\Final\Data\Anomalies'; % repertoire des anomalies
M=800;
nb_simul=200;


dirName=Repertoire; %# directory name
dirData = dir(dirName); %# Get the data for the current directory
dirIndex = [dirData.isdir]; %# Find the index for directories
fileList = {dirData(~dirIndex).name}'; %# Get a list of the files

dirName_ano=Repertoire_ano; %# directory name
dirData_ano = dir(dirName_ano); %# Get the data for the current directory
dirIndex_ano = [dirData_ano.isdir]; %# Find the index for directories
fileList_ano = {dirData_ano(~dirIndex_ano).name}'; %# Get a list of the files

mat={'id','station','riviere','ACF_l','ACF_moy','ACF_std','ACF_var'};
mat_ano={'id','station','riviere','ACF_l','ACF_moy','ACF_std','ACF_var'};
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
    
    X=(S-mean(S))/std(S);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ALGO SSA pour données brutes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%% Etape 1: reorganisation des données 
    if M > N; error(['Window length greater than series length! M should be less than N']); return; end;
    D = NaN .* ones(M,N-M+1);
    for i=1:M
        D(i,:)=X(1,i:N-M+i); 
    end 

    %%% Etape 2 : gestion des valeurs manquantes (cela devrait être bon)
    nnan = ~isnan(D); 
    xnan = ind2sub(size(nnan),nnan); 
    xsize = ((xnan*xnan')); 
    D(find(isnan(D))) = 0;
    
    %%% Etape 3: Calcul de la matrice de variance covariance (C)
    C = (D*D').*(1./xsize); 

    %%% Etape 4 : decomposition en valeur singulière de la matrice C
    [U,S,V] = svd(C);  
    
    %%% Etape 5: Classifification des valeurs manquantes en ordre décroissant
    evec = U; % vecteurs propres
    eval = diag(S);  % valeurs propres
    varexp = diag(S)/trace(S); % valeurs propres noramlisées
    [eval_sort,ind]=sort(varexp,'descend' ); % Valeurs propres triées en ordre
                                         % décroissant
                                         
    
    %  Monte Carlo (MC) based significance test for red noise spectrums
    temp=X(~isnan(X)); % eliminer les valeurs manquantes pour le calcul de l'autocorrelation
    [ACF,lags,bounds] = autocorr(temp,1);
    foac=ACF(2);
    lambda_red=[];
    for i=1:nb_simul
        a=randn(1,1)*10;
        randn('seed',floor(sum(10*clock)/a)); % pour changer la seed
        G(1,1) = randn(1,1);
        for k=1:N
            G(1,k+1) = (foac*G(1,k)) + ((1-foac)*randn(1,1));
        end
        G = G - repmat(mean(G),1,1);
        G = G/std(G);
        for j=1:N
            G(1,j) = G(1,j) .* sqrt(nanstd(X)^2);
        end
        for j=1:M
            Dr(j,:)=G(:,j:N-M+j); 
        end
        Cr = (Dr*Dr').*(1./(N-M+1));
        varexpr = diag(U'*Cr*U)'/trace(U'*Cr*U)';
        [evalr,ind]=sort(varexpr,'descend');
        lambda_red = [lambda_red;evalr] ;
        clear evalr Cr Dr evalr varexpr ind G a clock
    end

    inter_up=quantile(lambda_red,0.99);
    signif=[];
    for ind=1:M
        if eval_sort(ind)>inter_up(ind)
            signif=[signif,(eval_sort(ind)-inter_up(ind))];
        end
    end
    
    % Graphique des valeurs propres avec l'intervalle de confiance à 99% pour
    % du bruit blanc.
    name=fileList(f)
    name=regexp(name,'_','split');
    name_plot=['Station: ',char(name{1}(2)),', riviere: ', char(name{1}(3)),', id: ',char(name{1}(1))];
    figure(1);
    clf;
    subplot(1,2,1)
    plot(log(eval_sort),'x-');
    set(gca,'XLim',[1 length(eval_sort)]);
    hold on
    plot(log(quantile(lambda_red,0.99)),'r-');
    hold off
    xlabel('Rang');
    ylabel('log(Variance)');
    title(name_plot)
    grid;

    
    % Graphe des vecteurs propres en pair
    figure (2);
    clf;
    subplot(3,2,1);
    plot(evec(:,1:2), '-');
    legend('1', '2');
    title('Vecteurs Propres 1 & 2')
    subplot(3,2,3);
    plot(evec(:,3:4), '-');
    legend('3', '4');
    title('Vecteurs Propres 3 & 4')
    subplot(3,2,5);
    plot(evec(:,5:6), '-');
    legend('5', '6');
    title('Vecteurs Propres 5 & 6')
    
    
    clear X temp
    clear D nnan xnan xsize C U S V evec eval varexp eval_sort ind RC A
    clear ACF lags bounds foac lambda_red inter_up
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ALGO SSA pour anomalies
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
     X=(S_ano-mean(S_ano))/std(S_ano);
    %%% Etape 1: reorganisation des données 
    if M > N; error(['Window length greater than series length! M should be less than N']); return; end;
    D = NaN .* ones(M,N-M+1);
    for i=1:M
        D(i,:)=X(1,i:N-M+i); 
    end 
        
    %%% Etape 2 : gestion des valeurs manquantes (cela devrait être bon)
    nnan = ~isnan(D); 
    xnan = ind2sub(size(nnan),nnan); 
    xsize = ((xnan*xnan')); 
    D(find(isnan(D))) = 0;
    
    %%% Etape 3: Calcul de la matrice de variance covariance (C)
    C = (D*D').*(1./xsize); 
    
    %%% Etape 4 : decomposition en valeur singulière de la matrice C
    [U,S,V] = svd(C);  
       
    %%% Etape 5: Classifification des valeurs manquantes en ordre décroissant
    evec = U; % vecteurs propres
    eval = diag(S);  % valeurs propres
    varexp = diag(S)/trace(S); % valeurs propres noramlisées
    [eval_sort,ind]=sort(varexp,'descend' ); % Valeurs propres triées en ordre
                                         % décroissant
                                         
    
    %  Monte Carlo (MC) based significance test for red noise spectrums
    temp=X(~isnan(X)); % eliminer les valeurs manquantes pour le calcul de l'autocorrelation
    [ACF,lags,bounds] = autocorr(temp,1);
    foac=ACF(2);
    lambda_red=[];
    for i=1:nb_simul
        a=randn(1,1)*10;
        randn('seed',floor(sum(10*clock)/a)); % pour changer la seed
        G(1,1) = randn(1,1);
        for k=1:N
            G(1,k+1) = (foac*G(1,k)) + ((1-foac)*randn(1,1));
        end
        G = G - repmat(mean(G),1,1);
        G = G/std(G);
        for j=1:N
            G(1,j) = G(1,j) .* sqrt(nanstd(X)^2);
        end
        for j=1:M
            Dr(j,:)=G(:,j:N-M+j); 
        end
        Cr = (Dr*Dr').*(1./(N-M+1));
        varexpr = diag(U'*Cr*U)'/trace(U'*Cr*U)';
        [evalr,ind]=sort(varexpr,'descend');
        lambda_red = [lambda_red;evalr] ;
        clear evalr Cr Dr evalr varexpr ind G a clock
    end

    inter_up=quantile(lambda_red,0.99);
    signif_ano=[];
    for ind=1:M
        if eval_sort(ind)>inter_up(ind)
            signif_ano=[signif_ano,(eval_sort(ind)-inter_up(ind))];
        end
    end
    
    % Graphique des valeurs propres avec l'intervalle de confiance à 99% pour
    % du bruit blanc.
    figure(1)
    subplot(1,2,2)
    plot(log(eval_sort),'x-');
    set(gca,'XLim',[1 length(eval_sort)]);
    hold on
    plot(log(inter_up),'r-');
    hold off
    xlabel('Rang');
    ylabel('log(Variance)');
    title('Anomalies')
    grid;

    % Pour sauvegarder les images dans un repertoire
%     name_f=['D:\Inde\BVH\Results\2.SSA\Jeu Commun 2529\',char(name{1}(1)),'_',char(name{1}(2)),'_', char(name{1}(3)),'.jpg'];
% 	saveas(gcf,name_f)    
    
    
    figure (2);
    subplot(3,2,2);
    plot(evec(:,1:2), '-');
    legend('1', '2');
    title('Vecteurs Propres 1 & 2')
    subplot(3,2,4);
    plot(evec(:,3:4), '-');
    legend('3', '4');
    title('Vecteurs Propres 3 & 4')
    subplot(3,2,6);
    plot(evec(:,5:6), '-');
    legend('5', '6');
    title('Vecteurs Propres 5 & 6')
    
    % Pour sauvegarder les images dans un repertoire
%     name_f=['D:\Inde\BVH\Results\2.SSA\Jeu Commun 2529\',char(name{1}(1)),'_',char(name{1}(2)),'_', char(name{1}(3)),'_vec.jpg'];
% 	saveas(gcf,name_f) 
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Statistiques quantitatives
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    SSA_length=length(signif);
    SSA_moy=mean(signif);
    SSA_std=std(signif);
    SSA_var=var(signif);
    
    
    SSA_length_ano=length(signif_ano);
    SSA_moy_ano=mean(signif_ano);
    SSA_std_ano=std(signif_ano);
    SSA_var_ano=var(signif_ano);
    
    
    
    m={char(name{1}(1)),char(name{1}(2)),char(name{1}(3)),num2str(SSA_length),num2str(SSA_moy),num2str(SSA_std),num2str(SSA_var)};
    m_ano={char(name{1}(1)),char(name{1}(2)),char(name{1}(3)),num2str(SSA_length_ano),num2str(SSA_moy_ano),num2str(SSA_std_ano),num2str(SSA_var_ano)};
    
    mat=[mat;m];
    mat_ano=[mat_ano;m_ano];
    
    
    
    
    clear X time Debits data S_ano data_ano fileName_ano fileName
    clear N D nnan xnan xsize C U S V evec eval varexp eval_sort ind RC A
    clear ACF lags bounds foac lambda_red name temp nb_simul inter_up
    clear inter_up signif signif_ano
    clear SSA_length SSA_moy SSA_std SSA_var SSA_length_ano SSA_moy_ano SSA_std_ano SSA_var_ano m m_ano

end 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# Ecriture dans un fichier des statistiques 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fid=fopen('D:\Inde\BVH\Results\2.SSA\Jeu Commun 2529\SSA_quantitatif.txt','wt'); %# Nom du fichier
% [rows,cols]=size(mat);
% for i=1:rows
%     fprintf(fid,'%s\t',mat{i,1:end-1});
%     fprintf(fid,'%s\n',mat{i,end});
% end
% fclose(fid);
% clear fid
% 
% fid=fopen('D:\Inde\BVH\Results\2.SSA\Jeu Commun 2529\SSA_ano_quantitatif.txt','wt'); %# Nom du fichier
% [rows,cols]=size(mat_ano);
% for i=1:rows
%     fprintf(fid,'%s\t',mat_ano{i,1:end-1});
%     fprintf(fid,'%s\n',mat_ano{i,end});
% end
% fclose(fid);
% clear fid
%     
%     
