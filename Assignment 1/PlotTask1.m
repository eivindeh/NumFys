k_bT    = 26;
delta_U = 260;
alpha   = 0.2;
L       = 20;
r       = 12;
eta     = 1;

x   = linspace(-L+alpha*L,alpha*L,1000);
U_r = getU(x,alpha,L,delta_U);

r   = exp(-U_r./(k_bT))./(k_bT*(1-exp(-delta_U/(k_bT))));

fileID  = fopen('simOut.txt');
rawData = textscan(fileID,'%f %f %f');
time    = rawData{1,2};
pos     = rawData{1,1}*L;
U       = rawData{1,3}*delta_U;
figure(1);

plot(time,pos);
hold on;
[f,g]=hist(pos,50);

figure(2)
bar(g,f/trapz(g,f));
hold on;
plot(x,r/trapz(x,r));

