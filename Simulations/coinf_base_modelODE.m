function dy = coinf_base_modelODE(~, y, Lambda, mu, d_m, d_h, d_c, tau_m, ...
    tau_h, phi, del_m, del_c, psi)

dy = zeros(6,1);

S = y(1);
I_m = y(2);
I_h = y(3);
I_mh = y(4);
I_Rh = y(5);
R = y(6);

N = S + I_m + I_h + I_mh + I_Rh + R;

lam_m = psi*tau_m*(I_m + I_mh)/N;
lam_h = psi*tau_h*(I_h + I_mh + I_Rh)/N;

dy(1) = Lambda - (lam_m + lam_h + mu)*S;
dy(2) = lam_m*S - (lam_h + del_m + mu + d_m)*I_m;
dy(3) = lam_h*S - (phi*lam_m + mu + d_h)*I_h;
dy(4) = lam_h*I_m + phi*lam_m*I_h - (del_c + mu + d_c)*I_mh;
dy(5) = del_c*I_mh + lam_h*R - (mu + d_h)*I_Rh;
dy(6) = del_m*I_m - (lam_h + mu)*R;

end