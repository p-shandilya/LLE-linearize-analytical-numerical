function output = jverify(u)

global N
maxiter = 10;

fin0 = u;

r = [1e-4*rand(N,1);1e-4*rand(N,1);1e-4*rand];
fin1 = fin0+r;
feval1 = calF(fin1) - calF(fin0) - calJ(fin0)*(1*r);
iarr = linspace(1,maxiter,maxiter)';
rat = zeros(length(iarr),1);
sprintf("Beginning verification iterations...");
for iter = 1:maxiter
    disp(iter)
    fin = fin0 + iter*r;
    feval = calF(fin) - calF(fin0) - calJ(fin0)*(iter*r);
    rat(iter) = norm(feval)/norm(feval1);

end

plot(iarr,sqrt(rat))
p = polyfit(iarr,sqrt(rat),1);
output = p;
end