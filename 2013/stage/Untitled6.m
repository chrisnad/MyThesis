clear all
%total fit
figure(1)
clf
%part fit
figure(2)
clf
%parameter histograms
figure(3)
clf


x = [40; 80; 160];%; 320];
y = [2.031252; 1.808726; 1.780513];%;1.76818];


dy = y*0.02;

%formula converted to
%The inline version
func = inline('0.5*log(403.11*p(1)./(170.38*p(2)))+1./p(3)*log(1-(-985.036*p(3)*p(2))./(4*x*170.38))','p','x');  
%initial parameter guess
p0 = [42.6 25.7 0.11];

%To detect the sensitivity of the fit to starting parameter guess,
%the fit is run a number of times.
%each fit is plotted and each parameter plotted as a histogram
Nrepeat=100;
%each parameter is varied by a normal distribution with
%mean equal to the starting guess and std.dev. equal to
%sd*mean
sd = 0.3;
%histogram zoom factor (how many std dev to show)
zfactor = 2;
%parameter outlier cuttoff: lowest and highest N estimates are removed
outcut=10;
%========================================================
%END settable inputs
%========================================================

%list of all parameter outputs to use in histogram
pList=zeros(Nrepeat,3);

for rep =1:Nrepeat
    rep
    
    %form the new randomized start vector
    p = [p0(1)*(1+sd*randn), p0(2)*(1+sd*randn), p0(3)*(1+sd*randn)];
    %do the fit
    [p,r,j] = nlinfit(x,y,func,p);
    %copy fit to list
    pList(rep,:) = p';
    
    %get parameter errors
    c95 = nlparci(p,r,j);
    %conductance errors
    [yp, ci] = nlpredci(func,x,p,r,j);
    
    %plot the fit
    figure(1)
    errorbar(x,func(p,x),ci,ci,'b-');
    hold on
    errorbar(x,y,dy,dy,'ro')
    
    %plot the separated fits
    figure(2)
    subplot(2,1,1)
    hold on
    %errorbar(x, y-func(p,x)+ p(5)*exp((x-45)/p(6)),dy,dy,'rx')
    %plot(x, (y-func(p,x)+ p(5)*exp((x-45)/p(6))),'ro')
    %errorbar(x, p(5)*exp((x-45)/p(6)), 2*ci, 2*ci,'bx-')
    title('Exponential fit')
    
    subplot(2,1,2)
    hold on
    %plot(x, (y-func(p,x)+ p(3)./(1+exp((x-p(1))/p(2)))),'ro')
    errorbar(x, y-func(p,x)+ p(3)./(1+exp((x-p(1))/p(2))),dy,dy,'rx')
    errorbar(x, p(3)./(1+exp((x-p(1))/p(2))), 2*ci, 2*ci,'bx-')
    title('Boltzmann fit')
    
    %drawnow
end

figure(3)
%plot and print parameter table
fprintf('\r\rFit parameters and 95percent confidence range\r')
for i=1:6
    subplot(6,1,i)
    lowerLimit = mean(pList(:,i))-zfactor*std(pList(:,i));
    upperLimit = mean(pList(:,i))+zfactor*std(pList(:,i));
    hist(pList(:,i),linspace(lowerLimit,upperLimit,30))
    %
    fprintf('%7.3f\t +/- %7.3f \r',...
        mean(pList(:,i)),...
        max(2*std(pList(:,i)),mean(pList(:,i))-c95(i,1)));
end

fprintf('\r\rFit parameters omitting outliers\r')
for i=1:6
    %get rid of outliers
    pup = sort(pList(:,i));
    pup = pup(outcut:end-outcut);
    %print again
    fprintf('%7.3f\t +/- %7.3f \r',...
        mean(pup),...
        max(2*std(pup),mean(pup)-c95(i,1)));
    pbest(i)=mean(pup);
end

%print conductance table
%based on best parameters
v = [-30:5:45];
clear yp ci
[yp,ci] = nlpredci(func,x,pbest,r,j);
fprintf('\rVolt \t Total g\t Boltz\t Exp \r')
for i=1:length(v)
    fprintf('%7.3f\t%7.3f\t%7.3f\t%7.3f\r',...
        v(i),...
        yp(i),...
        pbest(3)./(1+exp((v(i)-pbest(1))/pbest(2))),...
        pbest(5)*exp((v(i)-45)/pbest(6)));
end