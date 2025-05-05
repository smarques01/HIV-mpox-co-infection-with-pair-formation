close all
clearvars

set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

% parameters
mu = 1/18250;
d_m = 1/17520;
d_h = 1/4745;
d_c = 1/4015;
tau_m = 0.785;
tau_h = 0.0021;
phi = 1; % increased risk for HIV-infectious ppl
del_m = 1/27;
del_c = 1/32;
psi = 1/3;


P = 412000;
Lambda = mu*P;

tf = 100*365;
t = 0:1:tf; % days

y0 = zeros(1,6);
y0(1, 1) = 378981; %S
y0(1, 2) = 60; %Im
y0(1, 3) = 32959; %Ih

odeopts = odeset('Nonnegative', 1:6, 'RelTol',1e-10,'AbsTol',1e-12);
[t,y] = ode23t(@(t,y) coinf_base_modelODE(t,y,Lambda, mu, d_m, d_h, d_c, tau_m,... 
tau_h, phi, del_m, del_c, psi), [0 tf], y0, odeopts); 

S = y(:,1);
I_m = y(:,2);
I_h = y(:,3);
I_mh = y(:,4);
I_Rh = y(:,5);
R = y(:,6);

I = I_m + I_h + I_mh + I_Rh;
N = S + I_m + I_h + I_mh + I_Rh + R;

figure(1)
plot(t/365, I_mh, 'LineWidth', 2)
xlabel('Time (years)')
ylabel('Co-infecteds')
fontsize(14, 'points')
xlim([0 50])

figure(2)
plot(t/365, (I_mh+I_m), 'LineWidth', 2)
xlabel('Time (years)')
ylabel('Mpox Infecteds (all)')
fontsize(14, 'points')
xlim([0 50])


figure(3)
plot(t/365, (I_mh+I_h+I_Rh), 'LineWidth', 2)
xlabel('Time (years)')
ylabel('HIV Infecteds (all)')
xlim([0 50])
fontsize(14, 'points')


