%modelfun = @(b,x)();

%modelfun = @(b,x)(0.5*log(40311*b(1)./(170.38*b(2)))+1./b(3)*log(1-(-985.036*b(3)*b(2)*b(2))./(4*(b(4)+x)*170.38*b(2))))

modelfun = @(b,x)(0.5*log(40311*b(1)./(170.38*b(2)))+1./b(3)*log(1-(-985.036*b(3)*b(2))./(4*x*170.38)))

%modelfun = @(b,x)(0.5*log(40311*b(1)./(170.38*b(2)))+1./0.11*log(1-(-985.036*0.11*b(2))./(4*x*170.38)))

x = [40; 80; 160]
y = [2.031252; 1.808726; 1.780513]

opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';


beta0 = [42.6;25.7;0.11];%;12.9];
beta = nlinfit(x,y,modelfun,beta0,opts)

 %y=modelfun(b,x)