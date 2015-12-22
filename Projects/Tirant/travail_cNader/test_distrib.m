clear all

nval = 1000;

% - Pour la résistance

data = importdata('./n_fiss.txt');
Xf = data(:,1);
Yf = data(:,2);

v = random('Unif',0,1,nval,1);
id = find(Yf~=0);

id1 = id(1);
id0 = id1-1;
id = [id0;id];

xp =Xf(id);
yp = Yf(id) / max(Yf(id));
for i=1:length(v)
    R(i) = interp1(yp,xp,v(i));
end

figure
interv = 100;
dPlay = 1;
if dPlay~=0
    %mx = floor(max(eta))+1;
    %mn = floor(min(eta));
    mx = max(R);
    mn = min(R);
    dm = (mx - mn)/interv;
    rgt = [mn:dm:mx];
    nEch_prec = 0;
    for i=1:(interv+1)
        nEch=sum(R<rgt(i));
        dist(i)=nEch-nEch_prec;
        nEch_prec=nEch;
    end
    dist = dist/sum(dist);
    bar(rgt,dist,'b');
    
    distcum(1) = dist(1);
    for i=2:length(dist)
        distcum(i) = distcum(i-1)+dist(i);
    end
    %bar(rgt,distcum,'b')
end
RTmin = min(R);
RTmax = max(R);
alp = (R-RTmin)/(RTmax-RTmin);


hold on
plot(xp,yp,'r-o')

% - Pour le module

nelt = 50;
nf = 30;
nfmin = 1;

Eb = 35000;
Es = 210000;
Ls = 1.75;

A = 0.8*.072;
As = 5* pi*(.012)^2/4;

Emba = (Eb*(A-As)+Es*As)/A;
Ems  = Es*As/A;
Emf = nf*Emba*Ems/(nelt*Emba-(nelt-nf)*Ems);
Emfm = nfmin*Emba*Ems/(nelt*Emba-(nelt-nfmin)*Ems);

%lambda = 1/(Emf-Emfm);

for i=1:length(v)
    %E(i) = Emfm-1/lambda*log(1-v(i));
    E(i) = 1/(1/Emfm-alp(i)*(1/Emfm-1/Emf));
end

figure
interv = 100;
dPlay = 1;
if dPlay~=0
    %mx = floor(max(eta))+1;
    %mn = floor(min(eta));
    mx = max(E);
    mn = min(E);
    dm = (mx - mn)/interv;
    rgt = [mn:dm:mx];
    nEch_prec = 0;
    for i=1:(interv+1)
        nEch=sum(E<rgt(i));
        dist(i)=nEch-nEch_prec;
        nEch_prec=nEch;
    end
    dist = dist/sum(dist);
    bar(rgt,dist,'b');
    
    distcum(1) = dist(1);
    for i=2:length(dist)
        distcum(i) = distcum(i-1)+dist(i);
    end
    %bar(rgt,distcum,'b');
end
