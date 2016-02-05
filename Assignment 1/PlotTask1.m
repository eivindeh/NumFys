k_bT    = 26;
delta_U = 80000;
alpha   = 0.2;
L       = 20;
r       = 12;
eta     = 1;
N       = 100000;
dt      = 1;
gamma_i = 6*pi*r*eta;
omega   = delta_U/(gamma_i*L^2)*1.6021*10^-4;
dt_hat  = omega*dt;
dt_hat  = 0.000000001;
D_hat   = k_bT/delta_U;


x   = linspace(-L+alpha*L,alpha*L,1000);
U_r = getU(x,alpha,L,delta_U);

r   = exp(-U_r./(k_bT))./(k_bT*(1-exp(-delta_U/(k_bT))));
% 
% fileID  = fopen('simOut.txt');
% rawData = textscan(fileID,'%f %f %f');
% time    = rawData{1,2};
% pos     = rawData{1,1}*L;
% U       = rawData{1,3}*delta_U;
% figure(1);

fileID2 = fopen('veldata2.txt');
rawData = textscan(fileID2,'%f');
vel = rawData{1,1};
tau = linspace(50,390,(400-50)/10);

plot(time,pos);
hold on;
[f,g]=hist(pos,50);

figure(2)
bar(g,f/trapz(g,f));
hold on;
plot(x,r/trapz(x,r));

figure(3)
tau=linspace(10,230,23);
plot(tau/141.47,vel*141,47);