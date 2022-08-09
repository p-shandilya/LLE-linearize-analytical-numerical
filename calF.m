function fout = calF(u)
global N
global theta
param


u0r = u(1:N);
u0m = u(N+1:2*N);
ts = u(2*N + 1);


D1 = nfdiff(eye(N),theta,1);
D2 = nfdiff(eye(N),theta,2);
D3 = nfdiff(eye(N),theta,3);

F1 = (1/2)*beta2*D2*u0m - (1/6)*beta3*D3*u0r - (u0r.^2 + u0m.^2).*u0m - u0r + alpha1*u0m + F + ts*D1*u0r;
F2 = -(1/2)*beta2*D2*u0r - (1/6)*beta3*D3*u0m + (u0r.^2 + u0m.^2).*u0r - u0m - alpha1*u0r + ts*D1*u0m;

fout = [F1;F2];
end