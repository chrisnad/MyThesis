D=[40 80 160];
sigma=[7.623625 6.102667 5.9329];

%plot(D,sigma,'*');

y = zeros(1,180);
t = 21:200;
b=[2.05 ; 0.165 ; 108.413];

for i=1:180
y(i) = (((1+(1./(b(2)+i/2./b(3)))))*b(1));
end

plot(D,sigma,'*');
hold all
plot(t,y);
grid on
hold off