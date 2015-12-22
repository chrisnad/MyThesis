D=[40 80 160];
sigma=[5.17 3.696 2.886];

%plot(D,sigma,'*');

y = zeros(1,180);
t = 21:200;

for i=1:180
y(i) = 5.019/sqrt(22.085*(i+20)/1000+0.2535151);
end

plot(D,sigma,'*');
hold all
plot(t,y);
grid on
hold off