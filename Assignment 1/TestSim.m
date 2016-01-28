k_bT    = 26;
delta_U = 260;
alpha   = 0.2;
L       = 20;
r       = 12;
eta     = 1;
dt      = 0.001;
N       = 100000;

gamma_i = 6*pi*r*eta;
D_hat   = k_bT/delta_U;

%Simple simulation
X=zeros(1,N);

for i = 1:N
    X(i+1)=X(i)-getF(X(i),alpha,L,delta_U)*dt/gamma_i+sqrt(2*k_bT/(gamma_i)*dt)*randn();
    X(i);
    F(i)=getF(X(i),alpha,L,delta_U);
end
x   = linspace(-L+alpha*L,alpha*L,1000);
U_r = getU(x,alpha,L,delta_U);

r   = exp(-U_r./(k_bT))./(k_bT*(1-exp(-delta_U/(k_bT))));
[f,g]=hist(X,100);
figure(2)
bar(g,f/trapz(g,f));
hold on;
plot(x,r/trapz(x,r));