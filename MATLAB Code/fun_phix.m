function phix = fun_phix(J,x,beta)

phix = 0.5;

for j = 2:J
    phix = phix + sqrt(2)*((j-1)*pi)^(-beta/2)*cos((j-1)*pi*x);
end