function [ U ] = getU( x, alpha, L, delta_U )
U=zeros(1,length(x));
for i = 1:length(x);
    
xmod = abs(x(i)-floor((x(i))/L)*L);
if xmod >= 0 && xmod < alpha*L
    U(i)=xmod/(alpha*L)*delta_U;
    
elseif xmod >= alpha*L && xmod < L
    U(i) = (L-xmod)/(L*(1-alpha))*delta_U;
end

end