xa=[1:60]';
ya=zeros(60,1);
xb=[63:130]';
yb=zeros(67,1);

for i=1:60
    ya(i)=V(60,i,5);
end

for i=1:68
    j=i+62;
    yb(i)=V(60,j,5);
end

pa = polyfit(xa,ya,1);

pb = polyfit(xb,yb,1);


for i=1:60
    ya(i)=pa(1)*xa(i)+pa(2);
end

for i=1:68
    yb(i)=pb(1)*xb(i)+pb(2);
end

hold on

plot(V(60,:,5))
plot(xa,ya,'r')
plot(xb,yb,'m')

hold off

