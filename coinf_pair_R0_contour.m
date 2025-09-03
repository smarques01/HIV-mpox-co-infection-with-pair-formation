clearvars
close all
clc

set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

% baseline parameters
mu = 1/(50*365); 
dm = 1/(48*365);
dh = 1/(13*365);
dc = 1/(11*365);
taum = 0.745;
tauh = 0.0021;
phi = 1.2; %1; increased risk for HIV-infectious ppl
delm = 1/27;
delc = 1/32;
psi = 1/3;


R = linspace(1/(3*365), 1/3, 1000); % rho
S = linspace(1/(3*365), 1/3, 1000); % sigma

%% from DFE
% for base model
base_R0m = taum*psi/(delm + mu + dm);
base_R0h = tauh*psi/(dh + mu);
base_R0 = max(base_R0m, base_R0h);

for i = 1 : length(R)
    rho = R(i);
    for j = 1 : length(S)
        sig = S(j);

        % R0 of submodels
        Km = psi*taum*(rho + delm + mu + dm) + (delm + dm + mu)*(rho + delm + dm + sig + 2*mu);
        R0m_den = Km*(2*delm + 2*dm + sig + 2*mu)*(delm + dm + sig + 2*mu);

        R0m_num = 2*psi*taum*rho*((sig + mu)*delm + (sig + mu + dm)*(delm + dm + sig + 2*mu));

        R0m(i,j) = R0m_num / R0m_den;

        R0h_den = (2*dh + sig + 2*mu)*(psi*tauh*(rho + mu + dh) + (mu + dh)*(rho + dh + sig + 2*mu));

        R0h(i,j) = (2*psi*tauh*rho*(sig + mu + dh)) / R0h_den;

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

        A = (th + fm + delc + dc)*(phi*tm + fc + dh)*(th + fc)*(ac*(tm + th + fc) - rho*b);

        B = 2*(th + fm + dh)*(bc*(fc + dh) + bh*delc)*(phi*(th + fm + delc + dc) + phi*tm + fc + dh) + ...
            bh*delc*(phi*tm + fc + dh)*(fc + delc + dc);

        p1 = rho*(bh*th*(th + fm + delc + dc)*(th + fc) + ...
            tm*(phi*tm + fc + dh)*(b*delm + bm*(th + fc)))/A;

        p0_n = (th + fc)*B + bh*delm*(th + fm + dh)*(phi*tm + fc + dh)*(fc + delc + dc);
        p0_d = A*(th + fm + dh)*(fc + delc + dc)*(fc + dh);
        p0 = rho*tm*th*p0_n/p0_d;

        R0c(i,j) = (p1 + sqrt(p1^2 + 4*p0))/2;

        R0(i,j) = max([R0m(i,j), R0h(i,j), R0c(i,j)]);
    end
end

%%
figure(1)
contourf(S, R, R0, 'ShowText', 'on', 'linewidth',1.25)
%title('$\mathcal{R}_0$', 'Interpreter', 'latex')
ylabel('$\rho$', 'Interpreter', 'latex')
xlabel('$\sigma$', 'Interpreter', 'latex')
fontsize(12,"points")
hold on
contour(S,R,R0,[1 1],'linewidth',1.25,'color','red')

range = 0:0.1:1.6; 
contourcmap('parula', range);

figure(2)
contourf(S, R, R0m, 'ShowText', 'on', 'linewidth',1.25)
%title('$\mathcal{R}_0^m$', 'Interpreter', 'latex')
ylabel('$\rho$', 'Interpreter', 'latex')
xlabel('$\sigma$', 'Interpreter', 'latex')
fontsize(12,"points")
hold on
contour(S,R,R0m,[1 1],'linewidth',1.25,'color','red')

range = 0:0.1:1.6; 
contourcmap('parula', range);

figure(3)
contourf(S, R, R0h, 'ShowText', 'on', 'linewidth',1.25)
%title('$\mathcal{R}_0^h$', 'Interpreter', 'latex')
ylabel('$\rho$', 'Interpreter', 'latex')
xlabel('$\sigma$', 'Interpreter', 'latex')
fontsize(12,"points")
hold on
contour(S,R,R0h,[1 1],'linewidth',1.25,'color','red')

range = 0:0.1:1.6; 
contourcmap('parula', range);


figure(4)
contourf(S, R, R0c, 'ShowText', 'on', 'linewidth',1.25)
%title('$\lambda^*$', 'Interpreter', 'latex')
ylabel('$\rho$', 'Interpreter', 'latex')
xlabel('$\sigma$', 'Interpreter', 'latex')
fontsize(12,"points")
hold on
contour(S,R,R0c,[1 1],'linewidth',1.25,'color','red')

range = 0:0.1:1.6; 
contourcmap('parula', range);

%% From MFE (HIV present)
% for base model
q0 = taum*(mu+dh)*(psi*tauh - mu - dh)/(tauh*(psi*tauh - dh + dm + delm)*(delc + mu + dc));
q1 = phi*taum*(psi*tauh - mu - dh)/(tauh*(delc + mu + dc)) + ...
    taum*(mu + dh)/(tauh*(psi*tauh - dh + delm + dm));

base_R0mh = (q1 + sqrt(q1^2 - 4*q0))/2 ; 

% %double checking with NGM
% S_b = (mu+dh)/(psi*tauh);
% Ih_b = (psi*tauh - mu - dh)/(psi*tauh);
% 
% F_b = [psi*taum*S_b, psi*taum*S_b;
%     psi*(tauh + phi*taum)*Ih_b, phi*taum*psi*Ih_b];
% V_b = [psi*tauh*Ih_b + delm + mu + dm, 0;
%     0, delc + mu + dc];
% b_ngm = F_b*inv(V_b);
% 
% b_rmh = max(real(eig(b_ngm)));

%%
% for pair model
for i = 1 : length(R)
    rho = R(i);
    for j = 1 : length(S)
        sig = S(j);

        % MFE
        eq_den = rho*(psi*tauh*(sig + mu + dh) - dh*(dh + mu));
        Seq(i,j) = ((psi*tauh + 2*mu + sig + dh)*(mu + dh)*(rho + sig + 2*dh + 2*mu))/eq_den;
        Iheq(i,j) = (-psi*tauh*(sig*(mu+dh-rho)+2*(mu+dh)^2) - (mu+dh)*(2*dh+sig+2*mu)*(rho+sig+dh+2*mu))/eq_den;

        %notation
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

        % NGM 
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
            -rho*Seq(i,j), 0, tm + fm, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            0, -rho*Seq(i,j), 0, tm + th + fc, 0, 0, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, fm + delm + dm, 0, 0, 0, 0, 0, 0, 0, 0;
            -rho*Iheq(i,j), 0, 0, 0, 0, th + phi*tm + fm + dh, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0, th + fm + delc + dc, 0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0, -delc, th + fm + dh, 0, 0, 0, 0, 0;
            0, 0, 0, 0, -2*delm, 0, 0, 0, fm, 0, 0, 0, 0;
            0, -rho*Iheq(i,j), 0, 0, 0, 0, 0, 0, 0, phi*tm + fc + dh, 0, 0, 0;
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, fc + delc + dc, 0, 0;
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2*delc, fc + dh, 0;
            0, 0, 0, 0, 0, 0, -delm, 0, 0, 0, 0, 0, th + fc];
        
        ngm = F*inv(V);

        R0mh(i,j) = max(real(eig(ngm)));
        
        R0h_den = (2*dh + sig + 2*mu)*(psi*tauh*(rho + mu + dh) + (mu + dh)*(rho + dh + sig + 2*mu));

        R0h(i,j) = (2*psi*tauh*rho*(sig + mu + dh)) / R0h_den;

    end
end

mat = R0h >= 1;
R0mh_new = mat .* R0mh;
R0mh_new(R0mh_new==0) = NaN;
%%
figure(5)
contourf(S, R, R0mh_new, 'ShowText', 'on', 'LineWidth',1.25)
hold on
contour(S,R,R0mh_new,[1 1],'linewidth',1,'color','red')
%title('$\mathcal{R}_i^{mh}$', 'Interpreter', 'latex')
ylabel('$\rho$', 'Interpreter', 'latex')
xlabel('$\sigma$', 'Interpreter', 'latex')
fontsize(12,"points")

range = 0:0.2:1.6; 
contourcmap('parula', range);

%%
figure(6)
new_mat = R0m ./ R0mh_new;
contourf(S, R, new_mat, 'ShowText', 'on', 'LineWidth',1.25)
hold on
contour(S,R,new_mat,[1 1],'linewidth',1,'color','red')
%title('$\mathcal{R}_0^m / \mathcal{R}_i^{mh}$', 'Interpreter', 'latex')
ylabel('$\rho$', 'Interpreter', 'latex')
xlabel('$\sigma$', 'Interpreter', 'latex')
fontsize(12,"points")

annotation("textbox",'interpreter','latex','string', '$ 0.9 < \mathcal{R}_0^m / \mathcal{R}_i^{mh} < 1 $', 'LineStyle','none')
xlim([1/(3*365) 0.05])

%% Looking at all of the regions
figure(7)
% [c1] = contour(S,R,new_mat,[1 1]);
% plot(c1(1,2:end), c1(2, 2:end), 'linewidth', 1, 'color', 'red')
% hold on
% patch([c1(1,2:end), fliplr(c1(1,2:end))], [c1(2, 2:end), 0.34*ones(1, length(c1(2, 2:end)))], 'r', 'EdgeColor', 'none')
% alpha(0.3)

[c3] = contour(S, R, R0h, [1 1]);
plot(c3(1,2:end), c3(2, 2:end), 'linewidth', 1.25, 'color', 'black')
hold on

[c4] = contour(S, R, R0mh_new, [1 1]);
plot(c4(1,2:end), c4(2, 2:end), 'linewidth', 1.25, 'color', 'cyan')
hold on
patch([c4(1,2:end), fliplr(c4(1,2:end))], [c4(2, 2:end), 0.34*ones(1, length(c4(2, 2:end)))], 'c', 'EdgeColor', 'b')
%alpha(0.3)
hold on

% patch([c3(1, 432:end), fliplr(c3(1, 432:end))], [c3(2, 432:end), 0.34*ones(1, length(c3(2, 432:end)))], 'c', 'EdgeColor', 'none')
% alpha(0.3)
% hold on
rectangle('Position',[0.209738,0.129544,0.14,0.28],'FaceColor','cyan','EdgeColor','none');
hold on

[c2] = contour(S, R, R0m, [1 1]);
plot(c2(1,2:end), c2(2, 2:end), 'linewidth', 1.25, 'color', [0.4660 0.6740 0.1880])
hold on
patch([c2(1,2:end), fliplr(c2(1,2:end))], [c2(2, 2:end), 0.34*ones(1, length(c2(2, 2:end)))], [0.4660 0.6740 0.1880], 'EdgeColor', [0.144 0.238 0.144])
%alpha(0.3)
hold on 

area(c3(1,2:end), c3(2, 2:end), 'FaceColor',[.7 .7 .7])
hold on 

ylabel('$\rho$', 'Interpreter', 'latex')
xlabel('$\sigma$', 'Interpreter', 'latex')
axis([1/(3*365) 1/3 1/(3*365) 1/3])

annotation("textbox",'interpreter','latex','string', '$\mathcal{R}_0^m < 1$,\\ $\mathcal{R}_i^{mh} > 1$', 'LineStyle','none')
annotation("textbox",'interpreter','latex','string', '$1 < \mathcal{R}_0^m < \mathcal{R}_i^{mh}$', 'LineStyle','none')
% annotation("textbox",'interpreter','latex','string', '$\mathcal{R}_0^m > \mathcal{R}_i^{mh} > 1$', 'LineStyle','none')
% annotation("textbox",'interpreter','latex','string', '$\mathcal{R}_0^m < 1$,\\ $\mathcal{R}_i^{mh} > 1$', 'LineStyle','none')
annotation("textbox",'interpreter','latex','string', '$ \mathcal{R}_0^m < \mathcal{R}_i^{mh} < 1$', 'LineStyle','none')
annotation("textbox",'interpreter','latex','string', '$\mathcal{R}_0^h < 1$', 'LineStyle','none')
fontsize(12,"points")

