function dy = coinf_pair_modelODE(~, y, mu, d_m, d_h, d_c, tau_m,... 
tau_h, phi, del_m, del_c, psi, rho, sigma)
%Lambda removed
dy = zeros(27,1);

S = y(1);
I_m = y(2);
I_h = y(3);
I_mh = y(4);
I_Rh = y(5);
R = y(6);

P_SS = y(7);
P_SIm = y(8);
P_SIh = y(9);
P_SImh = y(10);
P_SIRh = y(11);
P_SR = y(12);

P_ImIm = y(13);
P_ImIh = y(14);
P_ImImh = y(15);
P_ImIRh = y(16);
P_ImR = y(17);

P_IhIh = y(18);
P_IhImh = y(19);
P_IhIRh = y(20);
P_IhR = y(21);

P_ImhImh = y(22);
P_ImhIRh = y(23);
P_ImhR = y(24);

P_IRhIRh = y(25);
P_IRhR = y(26);

P_RR = y(27);

N = S + I_m + I_h + I_mh + I_Rh + R;
% N_t = N + P_SS + P_SIm + P_SIh + P_SImh + P_SIRh + P_SR + P_ImIm + P_ImIh + ...
%     P_ImImh + P_ImIRh + P_ImR + P_IhIh + P_IhImh + P_IhIRh + P_IhR + ...
%     P_ImhImh + P_ImhIRh + P_ImhR + P_IRhIRh + P_IRhR + P_RR;

%S
dy(1)= mu*412000 - (mu + rho)*S + (sigma + mu)*(2*P_SS + P_SIm + P_SIh + P_SImh + P_SIRh +  P_SR) + d_h*(P_SIRh + P_SIh) + d_c*P_SImh + d_m*P_SIm;
%I_m
dy(2) = -(rho + del_m + mu + d_m)*I_m + (sigma + mu)*(2*P_ImIm + P_SIm + P_ImIh + P_ImImh + P_ImIRh + P_ImR) + d_h*(P_ImIh + P_ImIRh) + d_c*P_ImImh + 2*d_m*P_ImIm;
%I_h
dy(3) = -(rho + mu + d_h)*I_h + (sigma + mu)*(2*P_IhIh + P_SIh + P_ImIh + P_IhImh + P_IhIRh + P_IhR) + d_h*(2*P_IhIh + P_IhIRh) + d_c*P_IhImh + d_m*P_ImIh;
%I_mh
dy(4) = -(rho + del_c + mu + d_c)*I_mh + (sigma + mu)*(2*P_ImhImh + P_SImh + P_ImImh + P_IhImh + P_ImhIRh + P_ImhR) + d_h*(P_IhImh + P_ImhIRh) + 2*d_c*P_ImhImh + d_m*P_ImImh;
%I_Rh
dy(5) = del_c*I_mh - (rho + mu + d_h)*I_Rh + (sigma + mu)*(2*P_IRhIRh + P_SIRh + P_ImIRh + P_IhIRh + P_ImhIRh + P_IRhR) + d_h*(P_IhIRh + 2*P_IRhIRh) + d_c*P_ImhIRh + d_m*P_ImIRh;
%R
dy(6) = del_m*I_m -(rho + mu)*R + (sigma + mu)*(2*P_RR + P_SR + P_ImR + P_IhR + P_ImhR + P_IRhR) + d_h*(P_IhR + P_IRhR) + d_c*P_ImhR + d_m*P_ImR;

%P_SS
dy(7) = rho*S^2/(2*N) - (sigma + 2*mu)*P_SS;
%P_SIm
dy(8) = rho*S*I_m/N - (psi*tau_m + del_m +d_m + sigma + 2*mu)*P_SIm;
%P_SIh
dy(9) = rho*S*I_h/N - (psi*tau_h + d_h + sigma + 2*mu)*P_SIh;
%P_SImh
dy(10) = rho*S*I_mh/N - (psi*(tau_h + tau_m) + del_c + d_c + sigma + 2*mu)*P_SImh;
%P_SIRh
dy(11) = rho*S*I_Rh/N + del_c*P_SImh - (psi*tau_h + d_h + sigma + 2*mu)*P_SIRh;
%P_SR
dy(12) = rho*S*R/N + del_m*P_SIm - (sigma + 2*mu)*P_SR;

%P_ImIm
dy(13) = rho*I_m^2/(2*N) + psi*tau_m*P_SIm - (2*del_m + 2*d_m + sigma + 2*mu)*P_ImIm;
%P_ImIh
dy(14) = rho*I_m*I_h/N - (psi*(tau_h + phi*tau_m) + del_m + d_m + d_h + sigma + 2*mu)*P_ImIh;
%P_ImImh
dy(15) = rho*I_m*I_mh/N + psi*tau_m*(P_SImh + phi*P_ImIh) - (psi*tau_h + del_m + del_c + d_m + d_c + sigma + 2*mu)*P_ImImh;
%P_ImIRh
dy(16) = rho*I_m*I_Rh/N + del_c*P_ImImh - (psi*tau_h + del_m + d_m + d_h + sigma + 2*mu)*P_ImIRh;
%P_ImR
dy(17) = rho*I_m*R/N + 2*del_m*P_ImIm - (del_m + d_m + sigma + 2*mu)*P_ImR;

%P_IhIh
dy(18) = rho*I_h^2/(2*N) + psi*tau_h*P_SIh - (2*d_h + sigma + 2*mu)*P_IhIh;
%P_IhImh
dy(19) = rho*I_h*I_mh/N + psi*tau_h*(P_SImh + P_ImIh) - (psi*phi*tau_m + del_c + d_h + d_c + sigma + 2*mu)*P_IhImh;
%P_IhIRh
dy(20) = rho*I_h*I_Rh/N + psi*tau_h*(P_SIRh + P_IhR) + del_c*P_IhImh - (2*d_h + sigma + 2*mu)*P_IhIRh;
%P_IhR
dy(21) = rho*I_h*R/N + del_m*P_ImIh - (psi*tau_h + d_h + sigma + 2*mu)*P_IhR;

%P_ImhImh
dy(22) = rho*I_mh^2/(2*N) + psi*tau_h*P_ImImh + psi*phi*tau_m*P_IhImh - (2*del_c + 2*d_c + sigma + 2*mu)*P_ImhImh;
%P_ImhIRh
dy(23) = rho*I_mh*I_Rh/N + psi*tau_h*(P_ImhR + P_ImIRh) + 2*del_c*P_ImhImh - (del_c + d_c + d_h + sigma + 2*mu)*P_ImhIRh;
%P_ImhR
dy(24) = rho*I_mh*R/N + del_m*P_ImImh - (psi*tau_h + del_c + d_c + sigma + 2*mu)*P_ImhR;

%P_IRhIRh
dy(25) = rho*I_Rh^2/(2*N) + psi*tau_h*P_IRhR + del_c*P_ImhIRh - (2*d_h + sigma + 2*mu)*P_IRhIRh;
%P_IRhR
dy(26) = rho*I_Rh*R/N + del_c*P_ImhR + del_m*P_ImIRh - (psi*tau_h + d_h + sigma + 2*mu)*P_IRhR;

%P_RR
dy(27) = rho*R^2/(2*N) + del_m*P_ImR - (sigma + 2*mu)*P_RR;



end