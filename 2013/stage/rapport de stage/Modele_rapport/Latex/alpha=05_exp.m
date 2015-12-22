D=[40 80 160];
sigma=[1.9473 1.2539 0.96689];

%plot(D,sigma,'*');

y = zeros(1,180);
t = 21:200;

for i=1:180
y(i) = 6*1.1556*10^5/sqrt(19.8577*(0.0011+0.1582*(i+20)/1000))/1000000;
end

plot(D,sigma,'*');
hold all
plot(t,y);
grid on
hold off