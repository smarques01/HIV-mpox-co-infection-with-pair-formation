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
tau_m = 0.745;
tau_h = 0.0021;
phi = 1; % increased risk for HIV-infectious ppl
del_m = 1/27;
del_c = 1/32;
psi = 1/3;

sig = [1/14, 1/10, 1/7, 1/5, 1/3]; 
rho = 1/7; % 1/7, 1/5, 1/3
P = 412000;
Lambda = mu*P;

tf = 100*365;
t = 0:1:tf; % days

%%
for i=1:length(sig)
    sigma = sig(i);

    % Reproduction number
    % R0 of submodels
    Km = psi*tau_m*(rho + del_m + mu + d_m) + (del_m + d_m + mu)*(rho + del_m + d_m + sigma  + 2*mu);
    R0m_den = (2*del_m + 2*d_m + sigma  + 2*mu)*(del_m + d_m + sigma  + 2*mu)*Km;
    R0m_num = 2*psi*tau_m*rho*((2*mu + sigma  + d_m)*(mu + sigma  + d_m) + del_m*(2*mu + 2*sigma  + d_m));
    R0m(i) = R0m_num / R0m_den;

    R0h_den = (2*d_h + sigma  + 2*mu)*(psi*tau_h*(rho + mu + d_h) + (mu + d_h)*(rho + d_h + sigma  + 2*mu));
    R0h(i) = 2*psi*tau_h*rho*(sigma  + mu + d_h) / R0h_den;
     
    % % analytical R0c (lambda* in overleaf)
    tm = psi*tau_m;
    th = psi*tau_h;

    am = rho + mu + del_m + d_m;
    ac = rho + mu + del_c + d_c;
    ah = rho + mu + d_h;

    b = sigma  + mu;
    bm = sigma  + mu + d_m;
    bc = sigma  + mu + d_c;
    bh = sigma  + mu + d_h;

    fm = del_m + d_m + sigma  + 2*mu;
    fc = del_c + d_c + sigma  + 2*mu;
    fh = d_h + sigma  + 2*mu;
    % 
    % A = (th + fm + del_c + d_c)*(phi*tm + fc + d_h)*(th + fc)*(ac*(tm + th + fc) - rho*b);
    % 
    % B = 2*(th + fm + d_h)*(bc*(fc + d_h) + bh*del_c)*(phi*(th + fm + del_c + d_c) + phi*tm + fc + d_h) + ...
    %     bh*del_c*(phi*tm + fc + d_h)*(fc + del_c + d_c);
    % 
    % p1 = rho*(bh*th*(th + fm + del_c + d_c)*(th + fc) + ...
    %     tm*(phi*tm + fc + d_h)*(b*del_m + bm*(th + fc)))/A;
    % 
    % p0_n = (th + fc)*B + bh*del_m*(th + fm + d_h)*(phi*tm + fc + d_h)*(fc + del_c + d_c);
    % p0_d = A*(th + fm + d_h)*(fc + del_c + d_c)*(fc + d_h);
    % p0 = rho*tm*th*p0_n/p0_d;
    % 
    % R0c(i) = (p1 + sqrt(p1^2 + 4*p0))/2;
    % 
    % R0(i) = max([R0m(i), R0h(i), R0c(i)]);
    % 
    % for base model
    base_R0m = tau_m*psi/(del_m + mu + d_m);
    base_R0h = tau_h*psi/(d_h + mu);
    base_R0 = max(base_R0m, base_R0h);


    % for naive population
    S_dfe = Lambda*(sigma + 2*mu)/(mu*(2*mu + rho + sigma));
    PSS_dfe = Lambda*rho/(2*mu*(2*mu + rho + sigma));

    y0 = zeros(1,27);

    y0(1,1) = 227388; % S
    y0(1,2) = 36; %I_m
    y0(1,3) = 19776; %I_h
    y0(1,7) = 73807; % P_SS
    y0(1,8) = 24; % P_SIm
    y0(1,9) = 3955; % P_SIh
    y0(1,18) = 4614; % P_IhIh


    % for MFE scenario
    % % if R0h(i) > 1
    % %     K = -mu^2*psi^2*rho*tau_h^2 + 2*mu^2*psi^2*d_h*tau_h^2 - 2*mu*psi^2*rho*sigma*tau_h^2 - ...
    % %         mu*psi^2*rho*d_h*tau_h^2 + 4*mu*psi^2*sigma*d_h*tau_h^2 + 2*mu*psi^2*d_h^2*tau_h^2 -...
    % %         psi^2*rho*sigma^2*tau_h^2 + psi^2*sigma^2*d_h*tau_h^2 + 2*psi^2*sigma*d_h^2*tau_h^2 +...
    % %         4*mu^3*psi*d_h*tau_h+ 2*mu^2*psi*rho*d_h*tau_h+ 10*mu^2*psi*sigma*d_h*tau_h+ ...
    % %         4*mu^2*psi*d_h^2*tau_h+ 4*mu*psi*rho*sigma*d_h*tau_h+ 2*mu*psi*rho*d_h^2*tau_h+ ...
    % %         6*mu*psi*sigma^2*d_h*tau_h+ 10*mu*psi*sigma*d_h^2*tau_h+ psi*rho*sigma^2*d_h*tau_h+ ...
    % %         2*psi*rho*sigma*d_h^2*tau_h+ psi*sigma^3*d_h*tau_h+ 3*psi*sigma^2*d_h^2*tau_h+ ...
    % %         2*psi*sigma*d_h^3*tau_h- 4*mu^3*d_h^2 - mu^2*rho*d_h^2 - 2*mu^2*sigma*d_h^2 - ...
    % %         6*mu^2*d_h^3 - mu*rho*d_h^3 - 2*mu*sigma*d_h^3 - 2*mu*d_h^4;
    % % 
    % %     S_num = 412000*mu*(sigma + 2*mu)*(psi*tau_h*(mu + sigma + d_h) - mu*d_h-d_h^2)*(psi*tau_h+ ...
    % %         2*mu + sigma + d_h);
    % %     Ih_num = 412000*mu*(sigma+2*mu)*(2*(mu^2)*psi*tau_h+mu*psi*sigma*tau_h+4*mu*psi*d_h*tau_h...
    % %         - psi*rho*sigma*tau_h+ psi*sigma*d_h*tau_h+ 2*psi*(d_h^2)*tau_h+ 4*mu^3 + ...
    % %         2*(mu^2)*rho + 4*(mu^2)*sigma + 10*(mu^2)*d_h+ mu*rho*sigma + 4*mu*rho*d_h+ ...
    % %         mu*sigma^2 + 7*mu*sigma*d_h+ 8*mu*d_h^2 + rho*sigma*d_h+ 2*rho*d_h^2 + (sigma^2)*d_h+ ...
    % %         3*sigma*d_h^2 + 2*d_h^3)*(mu*psi*tau_h+ psi*sigma*tau_h+ psi*d_h*tau_h- mu*d_h- d_h^2);
    % %     Ih_den = K*(mu + d_h)*(2*d_h+ 2*mu + rho + sigma);
    % %     PSS_num = (mu + d_h)*(psi*tau_h+ 2*mu + sigma + d_h)^2*(2*d_h+ 2*mu + rho + sigma);
    % %     PSIh_num = (sigma + 2*mu)*(psi*tau_h*(2*mu^2 + mu*sigma + 4*mu*d_h- rho*sigma + sigma*d_h+ 2*d_h^2) +...
    % %         (mu + d_h)*(2*d_h+ sigma + 2*mu)*(d_h+ 2*mu + rho + sigma));
    % %     PIhIh_num = (sigma + 2*mu)*(-psi*tau_h*(rho + mu + d_h) + (mu + d_h)*(d_h+ 2*mu + rho + ...
    % %         sigma))*(psi*tau_h*(2*mu^2 + mu*sigma + 4*mu*d_h- rho*sigma + sigma*d_h+ 2*d_h^2) +...
    % %         (mu + d_h)*(2*d_h+ sigma + 2*mu)*(d_h+ 2*mu + rho + sigma));
    % % 
    % %     S_eq = -S_num/K;
    % %     Ih_eq = Ih_num/Ih_den;
    % %     PSS_eq =  -412000*mu*PSS_num/(2*K);
    % %     PSIh_eq = 412000*mu*PSIh_num/K;
    % %     PIhIh_eq = -412000*mu*PIhIh_num/(2*K*(mu + d_h)*(2*d_h+ 2*mu + rho + sigma));

        % MFE perturbation
        % y0= zeros(1,27);
        % y0(1,1) = S_eq*0.999;
        % y0(1,2) = S_eq*0.001; % I_m
        % y0(1,3) = Ih_eq;
        % y0(1,7) = PSS_eq;
        % y0(1,9) = PSIh_eq;
        % y0(1,18) = PIhIh_eq;
      
  
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

        %total HIV only infected individuals (never had mpox)
        T_h = I_h + P_SIh + P_ImIh + 2*P_IhIh + P_IhImh + P_IhIRh + P_IhR;

        %total co-infected individuals
        T_c = I_mh + P_SImh + P_ImImh + P_IhImh + 2*P_ImhImh + P_ImhIRh + P_ImhR;

        %total HIV only infected individuals (mpox recovered)
        T_rh = I_Rh + P_SIRh + P_ImIRh + P_IhIRh + P_ImhIRh + 2*P_IRhIRh + P_IRhR;

        %total infected individuals
        T_I = T_m + T_h + T_rh + T_c;

        % total susceptibles
        T_s = S + 2*P_SS + P_SIm + P_SIh + P_SImh + P_SIRh + P_SR;

        % total recovereds (not including I_Rh)
        T_r = R + P_SR + P_ImR + P_IhR + P_ImhR + P_IRhR + 2*P_RR;

        % population size
        N_tot = T_I + T_s + T_r;

        % singles
        singles = S + I_m + I_h + I_mh + I_Rh + R;

        % pairs
        pairs = N_tot - singles;

        % plot coinfecteds
        figure(1)
        plot(t/365, T_c, 'LineWidth', 2)
        hold on

        % mpox (all)
        figure(2)
        plot(t/365, (T_m + T_c), 'LineWidth', 2)
        hold on

        % HIV (all)
        figure(3)
        plot(t/365, (T_h+T_rh+T_c), 'LineWidth', 2)
        hold on

      
        % %R0mh
        % eq_den = rho*(psi*tau_h*(sigma + mu + d_h) - d_h*(d_h + mu));
        % Seq = ((psi*tau_h + 2*mu + sigma + d_h)*(mu + d_h)*(rho + sigma + 2*d_h + 2*mu))/eq_den;
        % Iheq = (-psi*tau_h*(sigma*(mu+d_h-rho)+2*(mu+d_h)^2) - (mu+d_h)*(2*d_h+sigma+2*mu)*(rho+sigma+d_h+2*mu))/eq_den;
        % 
        % F = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        %     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        %     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        %     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        %     0, 0, tm, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        %     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        %     0, 0, 0, tm, 0, phi*tm, 0, 0, 0, 0, 0, 0, 0;
        %     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        %     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        %     0, 0, 0, th, 0, th, 0, 0, 0, 0, 0, 0, 0;
        %     0, 0, 0, 0, 0, 0, th, 0, 0, phi*tm, 0, 0, 0;
        %     0, 0, 0, 0, 0, 0, 0, th, 0, 0, 0, 0, th;
        %     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
        % 
        % V = [am, 0, -b, 0, -2*bm, -bh, -bc, -bh, -b, 0, 0, 0, 0;
        %     0, ac, 0, -b, 0, 0, -bm, 0, 0, -bh, -2*bc, -bh, -b;
        %     -rho*Seq, 0, tm + fm, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        %     0, -rho*Seq, 0, tm + th + fc, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        %     0, 0, 0, 0, fm + del_m + d_m, 0, 0, 0, 0, 0, 0, 0, 0;
        %     -rho*Iheq, 0, 0, 0, 0, th + phi*tm + fm + d_h, 0, 0, 0, 0, 0, 0, 0;
        %     0, 0, 0, 0, 0, 0, th + fm + del_c + d_c, 0, 0, 0, 0, 0, 0;
        %     0, 0, 0, 0, 0, 0, -del_c, th + fm + d_h, 0, 0, 0, 0, 0;
        %     0, 0, 0, 0, -2*del_m, 0, 0, 0, fm, 0, 0, 0, 0;
        %     0, -rho*Iheq, 0, 0, 0, 0, 0, 0, 0, phi*tm + fc + d_h, 0, 0, 0;
        %     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, fc + del_c + d_c, 0, 0;
        %     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2*del_c, fc + d_h, 0;
        %     0, 0, 0, 0, 0, 0, -del_m, 0, 0, 0, 0, 0, th + fc];
        % 
        % ngm = F*inv(V);
        % R0mh(i) = max(real(eig(ngm)));

    end
%end
%%
legendStr = "$\sigma =  $" + string(sig);

figure(1)
xlabel('Time (years)')
ylabel('Co-infecteds')
legend(legendStr)
fontsize(14, 'points')
%xlim([0 50])

figure(2)
xlabel('Time (years)')
ylabel('Mpox Infecteds (all)')
legend(legendStr)
fontsize(14, 'points')
%xlim([0 3])

figure(3)
xlabel('Time (years)')
ylabel('HIV Infecteds (all)')
legend(legendStr)
fontsize(14, 'points')
% xlim([0 50])



