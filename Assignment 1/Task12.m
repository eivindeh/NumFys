k_bT    = 26;
delta_U = 80000;
alpha   = 0.2;
L       = 20;
r       = 12;
eta     = 1;
N       = 100000;
dt      = 1;
gamma_i = 6*pi*r*eta;
omega   = delta_U/(gamma_i*L^2)*160.21;
dt_hat  = omega*dt;
dt_hat  = 0.000000001;
D_hat   = k_bT/delta_U;
M       = 200;
N       = 1000;

fileID  = fopen('simOut.txt');
rawData = textscan(fileID,'%f %f');
time    = rawData{1,2};
pos     = rawData{1,1};

for i = 1:N
    for j = 1:M
        time_small(j,i)=rawData{2}((i-1)*M+j);
        pos_small(j,i) =rawData{1}((i-1)*M+j);
        time_big(j,i)  =rawData{2}((i-1)*M+j+M*N);
        pos_big(j,i)   =rawData{1}((i-1)*M+j+M*N);
    end
end

%%
figure(2)
for i = 1:N
    [f,g]=hist(pos_small(:,i),50);
    %bar(g,f/trapz(g,f));
    %hold on;
    [a,b]=hist(pos_big(:,i),50);
    bar(b,a/trapz(b,a),'r');
    axis([-0.1 0.1 0 5])
    drawnow update;
    hold off
end
%%