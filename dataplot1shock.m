clear all
close all
clc
clear
% One gas plots 
Array=csvread('premetive.csv');
extra=csvread('extra.csv');
x = Array(:, 1);
rho = Array(:, 2);
u = Array(:, 3);
p = Array(:, 4);
Gamma = Array(:, 6);
T = p./(287.15*rho);
[m,n] = size(Array);
m
dx = extra(1,1);
t = extra(1,2);

% p2sim = 
for i = 2:m-1
	d2u(i) = (u(i+1)-(2*u(i))+u(i-1))/(dx*dx);
end

[Max,I] = max(d2u);
Fact = 0.01;
rho21Ratiosim = rho(I-floor(Fact*m))/rho(I+floor(Fact*m)) ; 
u21Ratiosim = u(I-floor(Fact*m))/u(I+floor(Fact*m)) ;
p21Ratiosim = p(I-floor(Fact*m))/p(I+floor(Fact*m)) ;	
T21Ratiosim = (p(I-floor(Fact*m))*rho(I+floor(Fact*m)))/(rho(I-floor(Fact*m))*p(I+floor(Fact*m))) ;

% Theory
rho4 = extra(1,3) ; 
u4 = extra(1,4) ;
p4 = extra(1,5) ;

rho1 = extra(1,6) ; 
u1 = extra(1,7) ;
p1 = extra(1,8) ;

T41Ratiotheo = p4*rho1/(p1*rho4);

fun = @(p21Ratiotheo) (p4/p1) - p21Ratiotheo*(1-(0.4*(sqrt(1/T41Ratiotheo))*(p21Ratiotheo-1))/sqrt(2.8*(2.8+(2.4)*(p21Ratiotheo-1))))^(-7.0); % function
p21Ratiotheo0 = [0 100]; % initial interval

% options = optimset('Display','iter'); % show iterations
% [p21Ratiotheo fval exitflag output] = fzero(fun,p21Ratiotheo0,options);

[p21Ratiotheo fval exitflag output] = fzero(fun,p21Ratiotheo0);

T21Ratiotheo = p21Ratiotheo*((6.0+p21Ratiotheo)/(1.0+6.0*p21Ratiotheo));
rho21Ratiotheo = p21Ratiotheo/T21Ratiotheo;

% simulation Shock Spped
ussim = ((I-floor(m/2))*dx)/t;
% theoretical Shock speed
usTheo = sqrt( (1.4*p1/rho1)*( 1 + ((1.4+1)/(2*1.4))*(p21Ratiotheo -1) ) )	;

p21RatioErr = (p21Ratiotheo- p21Ratiosim)*100/p21Ratiotheo
T21RatioErr = (T21Ratiotheo- T21Ratiosim)*100/T21Ratiotheo;
rho21RatioErr = (rho21Ratiotheo- rho21Ratiosim)*100/rho21Ratiotheo;
usErr = (usTheo- ussim)*100/usTheo;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable gamma
Array_gamma=csvread('premetive_gamma.csv');
extra_gamma=csvread('extra_gamma.csv');


x_gamma = Array_gamma(:, 1);
rho_gamma = Array_gamma(:, 2);
u_gamma = Array_gamma(:, 3);
p_gamma = Array_gamma(:, 4);
Gamma_gamma = Array_gamma(:, 6);
T_gamma = p_gamma./(287.15*rho_gamma);
[m_gamma,n_gamma] = size(Array_gamma);

dx_gamma = extra_gamma(1,1);
t_gamma = extra_gamma(1,2);

% p2sim = 
for i = 2:m_gamma-1
	d2u_gamma(i) = (u_gamma(i+1)-(2*u_gamma(i))+u_gamma(i-1))/(dx_gamma*dx_gamma);
end

[Max_gamma,I_gamma] = max(d2u_gamma);
Fact = 0.01;
rho21Ratiosim_gamma = rho_gamma(I_gamma-floor(Fact*m_gamma))/rho_gamma(I_gamma+floor(Fact*m_gamma)) ; 
u21Ratiosim_gamma = u_gamma(I_gamma-floor(Fact*m_gamma))/u_gamma(I_gamma+floor(Fact*m_gamma)) ;
p21Ratiosim_gamma = p_gamma(I_gamma-floor(Fact*m_gamma))/p_gamma(I_gamma+floor(Fact*m_gamma)) ;	
T21Ratiosim_gamma = (p_gamma(I_gamma-floor(Fact*m_gamma))*rho_gamma(I_gamma+floor(Fact*m_gamma)))/(rho_gamma(I_gamma-floor(Fact*m_gamma))*p_gamma(I_gamma+floor(Fact*m_gamma))) ;

% Theory
rho4_gamma = extra_gamma(1,3) ; 
u4_gamma = extra_gamma(1,4) ;
p4_gamma = extra_gamma(1,5) ;

rho1_gamma = extra_gamma(1,6) ; 
u1_gamma = extra_gamma(1,7) ;
p1_gamma = extra_gamma(1,8) ;

T41Ratiotheo_gamma = p4_gamma*rho1_gamma/(p1_gamma*rho4_gamma);

fun = @(p21Ratiotheo_gamma) (p4_gamma/p1_gamma) - p21Ratiotheo_gamma*(1-(0.4*(sqrt(1/T41Ratiotheo_gamma))*(p21Ratiotheo_gamma-1))/sqrt(2.8*(2.8+(2.4)*(p21Ratiotheo_gamma-1))))^(-7.0); % function
p21Ratiotheo0_gamma = [0 100]; % initial interval

% options = optimset('Display','iter'); % show iterations
% [p21Ratiotheo fval exitflag output] = fzero(fun,p21Ratiotheo0,options);

[p21Ratiotheo_gamma fval exitflag output] = fzero(fun,p21Ratiotheo0_gamma);

T21Ratiotheo_gamma = p21Ratiotheo_gamma*((6.0+p21Ratiotheo_gamma)/(1.0+6.0*p21Ratiotheo_gamma));
rho21Ratiotheo_gamma = p21Ratiotheo_gamma/T21Ratiotheo_gamma;

% simulation Shock Spped
ussim_gamma = ((I_gamma-floor(m_gamma/2))*dx_gamma)/t_gamma;
% theoretical Shock speed
usTheo_gamma = sqrt( (1.4*p1_gamma/rho1_gamma)*( 1 + ((1.4+1)/(2*1.4))*(p21Ratiotheo_gamma -1) ) )	;

p21RatioErr_gamma = (p21Ratiotheo_gamma- p21Ratiosim)*100/p21Ratiotheo_gamma
T21RatioErr_gamma = (T21Ratiotheo_gamma- T21Ratiosim_gamma)*100/T21Ratiotheo_gamma;
rho21RatioErr_gamma = (rho21Ratiotheo_gamma- rho21Ratiosim_gamma)*100/rho21Ratiotheo_gamma;
usErr_gamma = (usTheo_gamma- ussim_gamma)*100/usTheo_gamma;















%Only Temperature and gamma variation

% f=figure('units','normalized','outerposition',[0 0 1 1])
% subplot(2,1,1)
% plot(x, T_gamma, 'b--o','linewidth', 2)
% grid on 
% axis tight
% set(gca,'FontSize',20)
% xlabel('X')
% ylabel('Temperature')


% subplot(2,1,2)
% plot(x, Gamma_gamma, 'b--o','linewidth', 2)
% grid on 
% axis tight
% set(gca,'FontSize',20)
% xlabel('X')
% ylabel('Gamma')
% saveas(f,'GammaTemperature.png')

% With gamma change and without gamma change comparision
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,1,1)
plot(x, T, 'r--o',x_gamma, T_gamma, 'b--x','linewidth', 2)
legend('\gamma = 1.4', '\gamma = f(T)')
grid on 
axis tight
set(gca,'FontSize',20)
xlabel('X')
ylabel('Temperature')
