
function PTS = scancourb
%
ouvrir
%
Text = strvcat('Récupération de l''échelle de la figure ',...
  'Pointez sur l''origine, sur un point de l''axe des abscisses',...
  'et sur un point de l''axe des ordonnées (O,MX,MY)');
h = msgbox(Text,'Message');
uiwait(h);
ISTOP=1;
while ISTOP,
  [X,Y]=ginput(3);
  O  = [X(1),Y(1)];
  X1 = [X(2),Y(2)];
  Y1 = [X(3),Y(3)];
  OX1 = [X1(1)-O(1);X1(2)-O(2)];
  OY1 = [Y1(1)-O(1);Y1(2)-O(2)];
  NOX1 = sqrt(OX1'*OX1);
  NOY1 = sqrt(OY1'*OY1);
  Vtest = OX1'*OY1/(NOX1*NOY1);
  if abs(Vtest) > 1.e-02
     Text = strvcat('Attention les axes ne sont pas vraiment perpendiculaires !',...
        ['Angle entre les deux axes : ', num2str(acos(Vtest)*180/pi),'degré(s)'],...
        'Voulez-vous recommencer ?');
     butzoom = questdlg(Text,'Question','Oui','Non','Non');
     if strcmp(butzoom,'Oui')
        ISTOP=1;
     elseif strcmp(butzoom,'Non')
        ISTOP=0;
     end
  else
     ISTOP=0;
  end
end
%
prompt  = {'Echelle suivant (x):','Echelle suivant (y):'};
title   = 'Entrée des échelles';
lines= 1;
def     = {'1','1'};
answer  = inputdlg(prompt,title,lines,def);
XEch = str2num(answer{1,1});
YEch = str2num(answer{2,1});
% Changement de repère
L = NOX1;
co = OX1(1)/L;% = cos teta
si = OX1(2)/L;% = sin teta
PX1 = [co,si;si,-co];
% fin changement de repère
ECHX = XEch/NOX1;
ECHY = YEch/NOY1;
%
Text = strvcat('Sélectionnez maintenant autant de points que vous le désirez',...
  'pour construire votre courbe!',...
  'Tapez sur "Return" pour arrêter');
h = msgbox(Text,'Message');
uiwait(h);
[X,Y]=ginput;
%
% Construction et tracé de la courbe
%
[nl,nc]= size(X);
PTS=(PX1*([X,Y]-repmat(O,nl,1))')'.*repmat([ECHX,ECHY],nl,1);
fen = figure('Name','Courbe numérisée','Numbertitle','off','Visible','Off');
plot(PTS(:,1),PTS(:,2),'r-o');
set(fen,'Visible','On');
%
% Ecriture dans un fichier texte'
%
[newfile,newpath] = uiputfile('courbe.txt','Save file name');
NomfichierSauv = [newpath,newfile];
fid = fopen(NomfichierSauv,'wt');
fprintf(fid,'%22.15e   %22.15e \n',(PTS'));
fclose(fid);
%
return
%
function ouvrir
[fname,pname] = uigetfile('*.jpeg','Choisir le fichier *.jpeg de la courbe à numériser ...');
NomFichier = [pname,fname];
info = imfinfo(NomFichier);
[Img,map] = imread(NomFichier,info.Format);
courbimg = imagesc(Img);
colormap(map);
axis equal;
f  = uimenu('Label','Action');
f1 = uimenu(f,'Label','Zoom');
uimenu(f1,'Label','On','Callback','zoom on');
uimenu(f1,'Label','Off','Callback','zoom off');
return
%