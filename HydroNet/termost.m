function [] = termost()
Tmin = 360;
Tmax = 350;
y_end = 99;
x0 = (Tmax+Tmin)/2;
k = - log(y_end) / (Tmin - x0);
temp=300:400;
opening = 1./(1 + exp(-k * (temp - x0)));

figure
hold on
% plot(temp,opening,'r')
% 
% opening = 1-1./(1 + exp(-k * (temp - x0)));

plot(temp,opening,'b')

Tmin = 360;
Tmax = 350;
y_end = 99;
x0 = (Tmax+Tmin)/2;
k = - log(y_end) / (Tmin - x0);
opening = 0.001+1./(1 + exp(-k * (temp - x0)));

plot(temp,opening,'k')

% legend('Original','1-original','Teps invertidas')
legend('Original','1-original','Teps invertidas')
end