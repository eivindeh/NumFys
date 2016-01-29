x=1;
x_hat=x/L;
dt=0.0001;
omega   = delta_U/(gamma_i*L^2)*1.6021*10^-4;
dt_hat=omega*dt;

-getF(x,alpha,L,delta_U)*dt/gamma_i*1.6021*10^-4;
-getF(x_hat*L,alpha,L,delta_U)*L/delta_U*dt_hat*L;

-getF(x,alpha,L,delta_U)*dt/gamma_i*1.6021*10^-4+sqrt(2*k_bT/(gamma_i)*dt*1.6021*10^-4)*0.2
(-getF(x_hat*L,alpha,L,delta_U)/delta_U*dt_hat+sqrt(2*D_hat*dt_hat)*0.2)*L