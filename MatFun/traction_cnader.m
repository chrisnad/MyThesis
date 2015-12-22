function traction_cnader
%% Depouillement des resultats des calculs
% Cette fonction depouille les resultats FIDES. Elle permet de tracer les 
% courbes effort-deplacement ainsi que de calculer les moyenne et les ecart-type 
% de la resistance et la traction et de l'energie. 
%clc
clear all

%Donnees modifiables

%calcul de la population
petit=1.e-06;
calc='tirant_HA12_1.4.2_allong';          %  nom du fichier de calcul
sdir='ElementdeBase/tirant_HA_1.4.2_allong/';            %  nom du sous repertoire dans FIDES_resultats

idepl = 7;                 %  numero de col du fichier resultats pour depl en y
iforc = 4;                 %  numero de col du fichier resultats pour forc en y
%__________________ Fin donnees modifiables ________________

dirOS='/home/cnader/Documents/FIDES/FIDES_Resultats_calculs/';

%dirc=[dirOS,'/Doctorat/Projet/Modelisation/FIDES/FIDES_Resultats_calculs/']; % path du stockage des resultats
dirc=dirOS;

%Calcul du nombre de fichier .reac dans le repertoire actuel
cd([dirc,sdir])
d = dir('*.reac');
ncalc=numel(d);

W=zeros(1,ncalc);          % matrice des energies
WfromFit=zeros(1,ncalc);          % matrice des energies

xmin=0;
xmax=0;

%%CALCUL NUMERIQUE
%Etude prealable du nombre de pts necessaires a la rediscretisation des courbes
for iicalc=1:ncalc,
    nom=[dirc, sdir, calc,'-',num2str(iicalc),'.reac'];
    %lecture du fichier texte
    s=importdata(nom);
    dp=s.data(:,idepl);
    fo=s.data(:,iforc)/1000;
    
    % Pour tracer force vs deplacement
    dp=full(dp);
    fo=full(fo);
    depl(1:size(dp,1),iicalc)=dp;
    forc(1:size(fo,1),iicalc)=fo;
    [dp1,I,J]=unique(dp);
    fo1=fo(I);
    
    %----- Recuperation de la valeur maximale de la force et du deplacement
    %----- correspondant
    [Fmax,imax]=max(fo1);
    Dini(iicalc)=dp1(1);
    Dmax(iicalc)=dp1(imax);
    Dfin(iicalc)=dp1(numel(dp1));
end

    disp(' ');

%On verifie que les calculs ont bien ete realises jusqu'au bout ...
% (ils doivent tous avoir le meme Dfin)
[Dfi,idfi,jdfi] = unique(Dfin);
if (size(Dfi,2)~=1)
    disp('Attention : probleme dans les calculs');
    disp('Certains calculs n''ont pas abouti');
end

Dmax_val = mean(Dmax);
Dini_val = mean(Dini);
Dfin_val = Dfi(numel(Dfi));
rapD = floor((Dfin_val - Dmax_val)/(Dmax_val - Dini_val));
if (rapD==0), rapD=1;end
npts1 = 100;
npts2 = rapD*npts1;

ical = 0;
icalnc = 0;
for iicalc=1:ncalc,
    
    %Recuperation des resultats des calculs
    dp=depl(:,iicalc);
    fo=forc(:,iicalc);
    [dp1,I,J]=unique(dp);
    fo1=fo(I);
    %Rediscretisation
    Di=Dini(iicalc);
    Dm=Dmax(iicalc);
    Df=Dfin(iicalc);
    
    %Comptage des calculs invalides
    if (Df==Dm)
        icalnc=icalnc+1;
    end
    
    %Comptage des calculs valides et recevables pr courbe moyenne
    if (abs((Df-Dfin_val)/Dfin_val)<petit),
        ical=ical+1;
        lpts1 = (0:1:npts1)/npts1;
        lpts2 = (0:1:npts2)/npts2;
        dpmod1 = Di + lpts1 .* (Dm-Di);
        dpmod2 = Dm + lpts2 .* (Df-Dm);
        dpmod(:,ical) = [dpmod1';dpmod2'];
        fomod1=interp1(dp1,fo1,dpmod1');
        fomod2=interp1(dp1,fo1,dpmod2');
        fomod(:,ical)=[fomod1;fomod2];

    else
        disp(['=> calcul ', num2str(iicalc),' non recevable pour courbe moyenne']);
    end
end
ncal = ical;



%% Courbe barycentrique
%Calcul de la courbe barycentrique en force
dpmbar=mean(dpmod,2);
fombar=mean(fomod,2);

%% Graphique courbes numeriques
figNum=figure('name','Load-displacement num');                                    %Figure 1
plot(depl,forc,'Color','black');hold on
xlabel('Displacement (mm)')
ylabel('Load (KN)')
%xlim([0 limitex])
%ylim([0 5.5])
grid on;
saveas(figNum,'Courbes numeriques','jpg');
saveas(figNum,'Courbes numeriques','fig');
hold off;

%% Graphique courbes numeriques et barycentrique
%Calcul de la courbe barycentrique en contraintes
figNumetBary=figure('name','Average load-displacement');                          % Figure 2
p3=plot(depl,forc);
set(p3,'Color',[0.4,0.4,0.4]);
hold on
p4=plot(dpmbar,fombar,'--');
leg2 = legend(p4,'Numerical Mean Curve');
set(leg2,'Location','NorthEast');
%xlim([0 limitex])
%ylim([0 5.5])
set(p4,'LineWidth',1.5);
set(p4,'Color','black');
%xlabel('Displacement (mm)')
%ylabel('Stress (MPa)')

grid on;
saveas(figNumetBary,'Courbe numeriques et bary','jpg');
saveas(figNumetBary,'Courbe numeriques et bary','fig');
hold off;




return



