clearvars
close all
clc

set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

T = readtable("canadampox_data");
T(1:5,1:5)


figure(1)
scatter(T, "month_lab", 'cases', 'Marker', '.')

%%
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

sig = 1/3; %[1/14, 1/10, 1/7, 1/5, 1/3]; 
rho = 1/3; % 1/7, 1/5, 1/3

P = 412000;
Lambda = mu*P;

tf = 365;
t = 0:1:tf; % days

for i=1:length(sig)
    sigma = sig(i);

    y0 = zeros(1,27);
    y0(1,1) = 227388; % S
    y0(1,2) = 36; %I_m
    y0(1,3) = 19776; %I_h
    y0(1,7) = 73807; % P_SS
    y0(1,8) = 24; % P_SIm
    y0(1,9) = 3955; % P_SIh
    y0(1,18) = 4614; % P_IhIh

    
    odeopts = odeset('NonNegative', (1:27),'RelTol',1e-8,'AbsTol',1e-9);
    sol = ode45(@(t,y) coinf_pair_modelODE(t,y,mu, d_m, d_h, d_c, tau_m,...
        tau_h, phi, del_m, del_c, psi, rho, sigma  ), [0 tf], y0, odeopts);

    [y,dy] = deval(sol,t);

    S = y(1,:);
    I_m = y(2,:);
    I_h = y(3,:);
    I_mh = y(4,:);
    I_Rh = y(5,:);
    R = y(6,:);

    P_SS = y(7,:);
    P_SIm = y(8,:);
    P_SIh = y(9,:);
    P_SImh = y(10,:);
    P_SIRh = y(11,:);
    P_SR = y(12,:);

    P_ImIm = y(13,:);
    P_ImIh = y(14,:);
    P_ImImh = y(15,:);
    P_ImIRh = y(16,:);
    P_ImR = y(17,:);

    P_IhIh = y(18,:);
    P_IhImh = y(19,:);
    P_IhIRh = y(20,:);
    P_IhR = y(21,:);

    P_ImhImh = y(22,:);
    P_ImhIRh = y(23,:);
    P_ImhR = y(24,:);

    P_IRhIRh = y(25,:);
    P_IRhR = y(26,:);

    P_RR = y(27,:);

    %total mpox only infected individuals
    T_m = I_m + P_SIm + 2*P_ImIm + P_ImIh + P_ImImh + P_ImIRh + P_ImR;

    %total co-infected individuals
    T_c = I_mh + P_SImh + P_ImImh + P_IhImh + 2*P_ImhImh + P_ImhIRh + P_ImhR;


    % mpox (all)
    figure(2)
    plot(t, (T_c+T_m), 'LineWidth', 2)
    hold on
end


legendStr = "$\sigma =  $" + string(sig);

figure(2)
%plotting pair model
xlabel('Time (days)')
ylabel('Mpox Infecteds (all)')
legend(legendStr)

hold on

%plotting data
scatter(30*(0:height(T)-1), table2array(T(:, 'cases')), 'filled') %convert months to days (approx)
xlim([0 200])

hold on

%plotting base model
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

plot(t, (I_mh+I_m), 'LineWidth', 2)
ylim([0 5000])

% vax in july 2022
xline(90)

fontsize(14, 'points')

