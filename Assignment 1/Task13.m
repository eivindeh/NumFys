t = sym('t','positive');
syms x
L=20*10^-6;
gamma_i=6*pi*12*10^-9*1*10^-3;
D = 0.026*1.6*10^-19/gamma_i/(L*L);
n = 1/sqrt(4*pi*D*t)*exp(-x^2/(4*D*t));

P = int(n,x,0.2,inf)/(t)*3/4;
t=linspace(0,2,100);
plot(t,eval(P));
%ezplot(P);
hold on
fileID2 = fopen('velData.txt');
rawData = textscan(fileID2,'%f');
vel = rawData{1,1};
tau = linspace(10,300,290/2);
plot(tau/141,vel*141,'o');