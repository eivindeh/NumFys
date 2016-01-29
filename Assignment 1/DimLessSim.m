k_bT    = 26;
delta_U = 260;
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

%Simple simulation
M=1;
X=zeros(M,N);

for i = 1:N
    X(:,i+1)=(X(:,i)-getF(X(:,i)*L,alpha,L,delta_U)'/delta_U*dt_hat+sqrt(2*D_hat*dt_hat)*randn(M,1));
    %X(i);
    %F(i)=-getF(X(i)*L,alpha,L,delta_U);
end
X=X*L;
x   = linspace(-L+alpha*L,alpha*L,1000);
U_r = getU(x,alpha,L,delta_U);

r   = exp(-U_r./(k_bT))./(k_bT*(1-exp(-delta_U/(k_bT))));
[f,g]=hist(X,100);
figure(2)
bar(g,f/trapz(g,f));
hold on;
%plot(x,r/trapz(x,r));