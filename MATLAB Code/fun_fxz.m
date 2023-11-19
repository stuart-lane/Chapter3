function fxz = fun_fxz(J,x,z,alpha,n,delta)

fxz = 1;

for j = 2:J
    fxz = fxz + (sqrt(0.2)/((j-1)^(alpha/2)))*(sqrt(2)*cos((j-1)*pi*x))*(sqrt(2)*cos((j-1)*pi*z))/(n^(delta/2));
end


