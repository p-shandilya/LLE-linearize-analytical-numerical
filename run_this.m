close all
clearvars
clc

global N
global theta
param
N = 2^10; 
naxis = (-N/2:N/2-1).';
dth = 2*pi/N;
theta = naxis*dth; 
dw = 1; 
w = fftshift(dw*naxis); %"Mode number" axis
%u, phi and t_s needs to be calculated

%u0 = 1 + 1.1*exp(-(theta/0.02).^2);

%u0r = 5*sech(theta);
%u0m = zeros(N,1);
load uout4.mat
u0 = uout;
clear uout;

u0r = real(u0);
u0m = imag(u0);
ts = 0;
% u0r = u0(1:N);
% u0m = u0(N+1:2*N);
% ts = u0(end);

u = [u0r;u0m;ts];
%u = u+1e-0*rand(size(u)); %Adding random noise to test if code is actually
%working

%JACOBIAN VERIFICATION - UNCOMMENT THIS SECTION TO VERIFY JACOBIAN
%disp("Starting jverify")
%jva = jverify(u);
%drawnow()
%disp("jverify complete")
%JACOBIAN VERIFICATION

ucp = u;
up = u;



options = optimset('Display','iter','Jacobian','on','TolFun',1e-8,'TolX',...
                1e-8,'Algorithm','levenberg-marquardt','ScaleProblem','Jacobian',...
                'MaxIter',100,'LargeScale','on','UseParallel',true);

% optimset('Display','off','Jacobian','on','TolFun',1e-16,'TolX',...
%                 1e-14,'Algorithm','levenberg-marquardt','ScaleProblem','Jacobian',...
%                 'MaxIter',300,'LargeScale','on');

inp_fun = @(u) calFJ(u);
tic
[uout,fval] = fsolve(inp_fun,u,options);
toc
res = sum(sum(fval.*fval));
disp(res)

urout = uout(1:N);
uiout = uout(N+1:N*2);
u_ans = urout + 1i*uiout;
u_abs = abs(u_ans).^2;

figure(100)
plot(theta,u_abs)
drawnow()
save("uout_prev.mat",'uout')
disp(norm(calF(uout)))
[~] = linearize(uout);
uout_nl = uout;
uout = u_ans;
figure;
uf = fft(uout);
%uf = lin2dbm(uout,1);
yyaxis left
u_int_fft = abs(fftshift(uf)).^2;
u_norm = u_int_fft / max(u_int_fft);
u_db = 10*log10(u_norm);
figure
u_sm = stem(fftshift(w),u_db,'Marker','None','Color','b','LineWidth',2);
u_sm.BaseValue = -400;
%xlim([-40,40]);
%ylim([-70,0]);
xlabel("Relative mode number")
ylabel("Power (dB)")

uout = uout_nl;
save("uout_nl.mat","uout");