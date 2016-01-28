function [ U ] = getF( x, alpha, L, delta_U )
U=zeros(1,length(x));
for i = 1:length(x);
    
xmod = abs(x(i)-floor((x(i))/L)*L);
if xmod >= 0 && xmod < alpha*L
    U(i)=1/(alpha*L)*delta_U;
    
else 
    U(i) = -(1)/(L*(1-alpha))*delta_U;
end

end
