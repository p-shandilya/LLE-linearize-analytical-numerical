function jout = calJ(u)
global N
global theta
param
u0r = u(1:N);
u0m = u(N+1:2*N);
ts = u(2*N + 1);
D1 = nfdiff(eye(N),theta,1);
D2 = nfdiff(eye(N),theta,2);
D3 = nfdiff(eye(N),theta,3);

J11 = -(1/6)*beta3*D3 - diag(2*u0r.*u0m) - eye(N) + ts*D1;
J12 = (1/2)*beta2*D2 - diag(u0r.^2) - diag(3*u0m.^2) + alpha1*eye(N);
J21 = -(1/2)*beta2*D2 + diag(u0m.^2) + diag(3*u0r.^2) - alpha1*eye(N);
J22 = -(1/6)*beta3*D3 + diag(2*u0r.*u0m) - eye(N) + ts*D1;
J13 = D1*u0r;
J23 = D1*u0m;

jout = [J11,J12,J13;J21,J22,J23];


end