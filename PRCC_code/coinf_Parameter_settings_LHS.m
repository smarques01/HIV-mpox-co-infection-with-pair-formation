% PARAMETER BASELINE VALUES
mu = 1/(50*365); 
dm = 1/(48*365);
dh = 1/(13*365);
dc = 1/(11*365);
taum = 0.745;
tauh = 0.0021;
phi = 1; % increased risk for HIV-infectious ppl
delm = 1/27;
delc = 1/32;
psi = 1/3;
r = 1/10; %rho
sig = 1/10;
%%
% Parameter Labels 
 PRCC_var = {'d_m', 'd_h', 'd_c', '\tau_m', '\tau_h',...
     '\phi', '\delta_m', '\delta_c', '\psi','\rho', '\sigma'};

%% TIME SPAN OF THE SIMULATION
t_end = 365*40; % length of the simulations %50 years for MFE, 100 for other
tspan = 0:1:t_end;   % time points where the output is calculated

% INITIAL CONDITION FOR unpaired, disease free
% y0 = zeros(1,27);
% 
% unpaired, disease free
% y0(1,1) = 0.98; % S
% y0(1,2) = 0.01; % I_m
% y0(1,3) = 0.01; % I_h



% Variables Labels
y_var_label={'S';'I_m';'I_h';'I_{mh}';'I_{Rh}';'R'; 'P_{SS}'; 'P_{SI_m}';...
    'P_{SI_h}';'P_{SI_{mh}}';'P_{SI_{Rh}}';'P_{SR}';'P_{I_mI_m}';...
    'P_{I_mI_h}';'P_{I_mI_{mh}}';'_{I_mI_{Rh}}';'P_{I_mR}';'P_{I_hI_h}';...
    'P_{I_hI_{mh}}';'P_{I_hI_{Rh}}';'P_{I_hR}';'P_{I_{mh}I_{mh}}';...
    'P_{I_{mh}I_{Rh}}';'P_{I_{mh}R}';'P_{I_{Rh}I_{Rh}}';'P_{I_{Rh}R}';
    'P_{RR}'};
