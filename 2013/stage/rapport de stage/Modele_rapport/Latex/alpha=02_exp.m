D=[40 80 160];
sigma=[3.489 2.83667 2.53];

%plot(D,sigma,'*');

y = zeros(1,180);
t = 21:200;

for i=1:180
y(i) = 6*1.7045*10^5/sqrt(3.40579*(0.0203+0.1797*(i+20)/1000))/1000000;
end

plot(D,sigma,'*');
hold all
plot(t,y);
grid on
hold off