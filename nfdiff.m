function df = nfdiff(x,t,n)

N = length(x);
L = N*(t(2)-t(1));
kappa = (2*pi/L)*(-N/2:N/2-1);
kappa = fftshift(kappa);
kappa = 1i*diag(kappa);
F = dftmtx(N);
I = F^-1;


for iter = 1:n
   df = real(I*kappa*F*x);
   x = df;
end


end