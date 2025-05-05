clc  

alpha = 0.05;
runs = 10750;

coinf_Parameter_settings_LHS;

% LHS matrix
dm_LHS = LHS_Call(1/(50*365), dm, 0.0001, 0, runs,'unif'); 
dh_LHS = LHS_Call(1/(15*365), dh, 1/(5*365), 0, runs,'unif'); 
dc_LHS = LHS_Call(1/(13*365), dc, 0.0005, 0, runs,'unif'); 
taum_LHS = LHS_Call(0.1, taum, 0.8, 0, runs,'unif');
tauh_LHS = LHS_Call(0.0006, tauh, 0.0338, 0, runs,'unif'); 
phi_LHS = LHS_Call(1, phi, 1.2, 0, runs, 'unif');
delm_LHS = LHS_Call(1/28, delm, 1/14, 0, runs, 'unif'); 
delc_LHS = LHS_Call(1/50, delc, 1/21, 0, runs, 'unif'); 
psi_LHS = LHS_Call(0.25, psi, 1, 0, runs, 'unif');
rho_LHS = LHS_Call(1/(3*365), r, 1/3, 0, runs, 'unif');
sig_LHS = LHS_Call(1/(3*365), sig, 1/3, 0, runs, 'unif');


LHSmatrix=[dm_LHS dh_LHS dc_LHS taum_LHS tauh_LHS ...
   phi_LHS delm_LHS delc_LHS psi_LHS rho_LHS sig_LHS];
 
%% R0 and R0mh PRCC (from DFE and MFE, resp.)
for i=1:size(LHSmatrix,1)
    mu = 1/(365*50);
    dm = LHSmatrix(i,1);
    %dh = LHSmatrix(i,2);
    dh = 1/4745; 
    dc = LHSmatrix(i,3);
    taum = LHSmatrix(i,4);
    tauh = LHSmatrix(i,5);
    phi = LHSmatrix(i,6);
    delm = LHSmatrix(i,7);
    delc = LHSmatrix(i,8);
    psi= LHSmatrix(i,9);
    rho = LHSmatrix(i,10);
    %sig = LHSmatrix(i,11);
    sig = 0.005;
    

    % Km = psi*taum*(rho + delm + mu + dm) + (delm + dm + mu)*(rho + delm + dm + sig + 2*mu);
    % R0m_den = (2*delm + 2*dm + sig + 2*mu)*(delm + dm + sig + 2*mu)*Km;
    % R0m_num = 2*psi*taum*rho*((2*mu + sig + dm)*(mu + sig + dm) + delm*(2*mu + 2*sig + dm));
    % 
    % R0m(i) = R0m_num / R0m_den;


    R0h_den = (2*dh + sig + 2*mu)*(psi*tauh*(rho + mu + dh) + (mu + dh)*(rho + dh + sig + 2*mu));
    R0h(i) = 2*psi*tauh*rho*(sig + mu + dh) / R0h_den;
    
    % analytical R0c (lambda* in overleaf)
    tm = psi*taum;
    th = psi*tauh;

    am = rho + mu + delm + dm;
    ac = rho + mu + delc + dc;
    ah = rho + mu + dh;

    b = sig + mu;
    bm = sig + mu + dm;
    bc = sig + mu + dc;
    bh = sig + mu + dh;

    fm = delm + dm + sig + 2*mu;
    fc = delc + dc + sig + 2*mu;
    fh = dh + sig + 2*mu;

    % A = (th + fm + delc + dc)*(phi*tm + fc + dh)*(th + fc)*(ac*(tm + th + fc) - rho*b);
    % 
    % B = 2*(th + fm + dh)*(bc*(fc + dh) + bh*delc)*(phi*(th + fm + delc + dc) + phi*tm + fc + dh) + ...
    %     bh*delc*(phi*tm + fc + dh)*(fc + delc + dc);
    % 
    % p1 = rho*(bh*th*(th + fm + delc + dc)*(th + fc) + ...
    %     tm*(phi*tm + fc + dh)*(b*delm + bm*(th + fc)))/A;
    % 
    % p0_n = (th + fc)*B + bh*delm*(th + fm + dh)*(phi*tm + fc + dh)*(fc + delc + dc);
    % p0_d = A*(th + fm + dh)*(fc + delc + dc)*(fc + dh);
    % p0 = rho*tm*th*p0_n/p0_d;
    % 
    % R0c = (p1 + sqrt(p1^2 + 4*p0))/2;
    % 
    % R0(i) = max([R0m(i), R0h(i), R0c]);
   
    %filtering for R0m
    % v1 = 1 + (delm*(delm+dm+mu)/fm^2);
    % v2 = ((b*delm)/fm) + bm;
    % v3 = 1/(fm + delm + dm);
    % v4 = (delm + dm + mu)/Km;
    % 
    % if v1 < v2*(v3 + v4)
    %     R0m_dec(i) = R0m(i);
    %     R0m_inc(i) = -1; %indicator 
    % else
    %     R0m_dec(i) = -1;
    %     R0m_inc(i) = R0m(i);
    % end


    %R0mh
    eq_den = rho*(psi*tauh*(sig + mu + dh) - dh*(dh + mu));
    Seq = ((psi*tauh + 2*mu + sig + dh)*(mu + dh)*(rho + sig + 2*dh + 2*mu))/eq_den;
    Iheq = (-psi*tauh*(sig*(mu+dh-rho)+2*(mu+dh)^2) - (mu+dh)*(2*dh+sig+2*mu)*(rho+sig+dh+2*mu))/eq_den;

    F = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, tm, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, tm, 0, phi*tm, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, th, 0, th, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, th, 0, 0, phi*tm, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, th, 0, 0, 0, 0, th;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];

    V = [am, 0, -b, 0, -2*bm, -bh, -bc, -bh, -b, 0, 0, 0, 0;
        0, ac, 0, -b, 0, 0, -bm, 0, 0, -bh, -2*bc, -bh, -b;
        -rho*Seq, 0, tm + fm, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, -rho*Seq, 0, tm + th + fc, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, fm + delm + dm, 0, 0, 0, 0, 0, 0, 0, 0;
        -rho*Iheq, 0, 0, 0, 0, th + phi*tm + fm + dh, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, th + fm + delc + dc, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, -delc, th + fm + dh, 0, 0, 0, 0, 0;
        0, 0, 0, 0, -2*delm, 0, 0, 0, fm, 0, 0, 0, 0;
        0, -rho*Iheq, 0, 0, 0, 0, 0, 0, 0, phi*tm + fc + dh, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, fc + delc + dc, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2*delc, fc + dh, 0;
        0, 0, 0, 0, 0, 0, -delm, 0, 0, 0, 0, 0, th + fc];

    ngm = F*inv(V);
    R0mh(i) = max(real(eig(ngm)));
end

%% only used for R0mh
idx = find(R0h > 1);
LHS_filt = LHSmatrix(idx, :);
LHS_filt_sig = LHS_filt(:, [1 3:10]);
% R0h_filt = R0h(idx);
R0mh_filt = R0mh(idx);

%%
% sigma increasing wrt R0m
idx_R0m_dec = find(R0m_dec >= 0);
idx_R0m_inc = find(R0m_inc >= 0);

R0m_dec_filt = R0m_dec(idx_R0m_dec);
R0m_inc_filt = R0m_inc(idx_R0m_inc);

R0m_dec_LHS = LHSmatrix(idx_R0m_dec, :);
R0m_inc_LHS = LHSmatrix(idx_R0m_inc, :);

sig_R0m_dec = sig_LHS(idx_R0m_dec);
sig_R0m_inc = sig_LHS(idx_R0m_inc);


[prcc{1},sign{1},sign_label{1}] = PRCC(R0m_dec_LHS,R0m_dec_filt,1,PRCC_var,alpha);
[prcc{2},sign{2},sign_label{2}] = PRCC(R0m_inc_LHS,R0m_inc_filt,1,PRCC_var,alpha);

%%
sig_star = -(dh_LHS+mu)+sqrt((dh_LHS+mu).^2+psi_LHS.*tauh_LHS.*(rho_LHS+mu+dh_LHS)+(mu+dh_LHS).*(rho_LHS-dh_LHS));

%sigma monotonically increasing wrt R0h
idx_mon = find(sig_LHS < sig_star);
pos_mon_LHS = LHSmatrix(idx_mon, :);
sig_mon = sig_LHS(idx_mon);
[prcc{3},sign{3},sign_label{3}] = PRCC(pos_mon_LHS,R0h(idx_mon),1,PRCC_var,alpha);

%sigma monotonically descreasing
idx_sig_dec = find(sig_LHS > sig_star);
dec_mon_LHS = LHSmatrix(idx_sig_dec, :);
sig_mon_dec = sig_LHS(idx_sig_dec);
[prcc{4},sign{4},sign_label{4}] = PRCC(dec_mon_LHS,R0h(idx_sig_dec),1,PRCC_var,alpha);

%%
PRCC_var = {'$d_m$', '$d_c$','$\tau_m$'	,'$\tau_h$'	,'$\phi$'	,'$\delta_m$',	'$\delta_c$',	'$\psi$',	'$\rho$'};
[prcc{1}, sign{1}, sign_label{1}] = PRCC(LHS_filt_sig, R0mh_filt, 1, PRCC_var, alpha);
