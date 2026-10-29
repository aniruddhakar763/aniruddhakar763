%% Blood-Fe3O4 nanofluid in a squeezing microchannel (electro-osmotic, variable viscosity)
%  Solves eqs. (21)-(34) of the model with bvp4c:
%    (21)  momentum            -> f''''
%    (22)-(23) electro-osmotic potential psi, psi'
%    (24)  pressure            -> dp/deta  (integrated with p(0) = 0)
%    (25)  temperature         -> theta''
%    (26)  concentration       -> phi''
%    (27)-(34) boundary conditions
%  for Bi = [0.001 0.01 0.1 1] at omega = 0.02 and plots
%  f, f', theta, phi, p in figures 1-5.
%  Plain script with anonymous functions only (no function definitions), so it
%  runs in any MATLAB version, under any file name, and section by section.
clear; clc; close all;

%% ---------------- Parameters ----------------
P.St    = -1;     % squeeze number (positive or negative)
P.lam   = 5;      % mixed convection (Gr/Re)
P.N     = 0.1;    % buoyancy ratio
P.delta = 0.05;   % aspect ratio
P.Pr    = 21.0;   % Prandtl number of blood
P.Nt    = 0.2;    % thermophoresis
P.Nb    = 0.2;    % Brownian motion
P.gam   = 0.5;    % Joule heating
P.Sc    = 10.0;   % Schmidt number
Bi      = [0.001, 0.01, 0.1, 1];   % Biot number
P.K     = 10;     % electro-osmotic parameter m (= K)
P.Ue    = 1.0;    % electric field parameter
P.R     = 1.5;    % zeta-potential ratio zeta_u/zeta_l
P.beta  = 0.1;    % variable-viscosity parameter
omega   = 0.02;   % nanoparticle volume fraction

%% ---------------- Thermophysical properties (eqs. 1-20) ----------------
% Base fluid: blood         Nanoparticle: Fe3O4
rho_f = 1050;    rho_p = 5180;     % kg/m^3
cp_f  = 3900;    cp_p  = 670;      % J/(kg K)
k_f   = 0.492;   k_p   = 6.0;      % W/(m K)
sig_f = 0.8;     sig_p = 2500;     % S/m
bt_f  = 3.8e-4;  bt_p  = 1.3e-5;   % 1/K (thermal expansion)

n = 2;  w = omega;                 % shape factor, volume fraction
rho_nf   = (1-w)*rho_f + w*rho_p;
rhocp_nf = (1-w)*rho_f*cp_f + w*rho_p*cp_p;
bt_nf    = (1-w)*bt_f + w*bt_p;
k_nf   = k_f  *(k_p + (n-1)*k_f - (n-1)*w*(k_p-k_f)) ...
              /(k_p + (n-1)*k_f + w*(k_p-k_f));
sig_nf = sig_f*(sig_p + (n-1)*sig_f - (n-1)*w*(sig_p-sig_f)) ...
              /(sig_p + (n-1)*sig_f + w*(sig_p-sig_f));

P.muf = (1-w)^(-2.5);              % Brinkman factor 1/(1-omega)^2.5
P.L1  = rho_nf/rho_f;              % density ratio
P.L2  = bt_nf/bt_f;                % thermal-expansion ratio
P.L3  = rhocp_nf/(rho_nf*cp_f);    % specific-heat ratio
P.L4  = k_nf/k_f;                  % conductivity ratio
P.L5  = sig_nf/sig_f;              % electrical-conductivity ratio
fprintf('omega = %4.2f :  L1=%.4f  L2=%.4f  L3=%.4f  L4=%.4f  L5=%.4f\n', ...
        omega, P.L1, P.L2, P.L3, P.L4, P.L5);

%% ---------------- ODE system and boundary conditions ----------------
% y = [f, f', f'', f''', theta, theta', phi, phi', p]
% (25) temperature -> theta''
thppF = @(x,y) ( P.L1*P.L3*P.Pr*P.St*(x*y(6) - y(1)*y(6)) ...
               - P.Pr*(P.Nt*y(6)^2 + P.Nb*y(8)*y(6)) - P.L5*P.gam ) / P.L4;
% (26) concentration -> phi''
phppF = @(x,y,thpp) P.St*P.Sc*(x*y(8) - y(1)*y(8)) - (P.Nt/P.Nb)*thpp;
% (23) electro-osmotic potential gradient psi'
psipF = @(x) P.K*( P.R*cosh(P.K*x) - cosh(P.K*(1-x)) )/sinh(P.K);
% exp(-beta*theta)/(1-omega)^2.5
viscF = @(y) exp(-P.beta*y(5))*P.muf;
% (21) momentum -> f''''
f4F   = @(x,y,thpp) ( P.L1*P.St*(3*y(3) + x*y(4) + y(2)*y(3) - y(1)*y(4)) ...
                    - P.L1*P.L2*P.lam*P.delta*(y(6) + P.N*y(8)) ...
                    - P.delta*P.Ue*P.K^2*psipF(x) )/viscF(y) ...
                    + 2*P.beta*y(6)*y(4) + P.beta*thpp*y(3) - P.beta^2*y(6)^2*y(3);
% (24) pressure -> dp/deta
dpF   = @(x,y) -P.L1*P.St*(-y(1) - x*y(2) + y(1)*y(2)) + viscF(y)*(2*P.beta*y(6)*y(2) - y(3));

odeR   = @(x,y,thpp) [y(2); y(3); y(4); f4F(x,y,thpp); y(6); thpp; y(8); phppF(x,y,thpp); dpF(x,y)];
odefun = @(x,y) odeR(x, y, thppF(x,y));

% (27)-(34) + reference pressure p(0) = 0
slipF = @(th) P.delta*P.Ue*exp(P.beta*th)/P.muf;   % delta*Ue*e^{beta theta}(1-omega)^2.5
bcfun = @(ya,yb,BiV) [ ya(2) - 1 - slipF(ya(5))*ya(3);    % (27) f'(0)
                       ya(1);                             % (28) f(0) = 0
                       P.L4*ya(6) + BiV*(1 - ya(5));      % (29) L4 theta'(0) = -Bi(1-theta(0))
                       ya(7);                             % (30) phi(0) = 0
                       yb(2) - slipF(yb(5))*yb(3);        % (31) f'(1)
                       yb(1) - 1;                         % (32) f(1) = 1
                       yb(5);                             % (33) theta(1) = 0
                       yb(7) - 1;                         % (34) phi(1) = 1
                       ya(9) ];                           % p(0) = 0

%% ---------------- Solve for each Bi ----------------
eta   = linspace(0,1,801);
cols  = lines(numel(Bi));
opts  = bvpset('RelTol',1e-6,'AbsTol',1e-8,'NMax',20000);
sol   = cell(size(Bi));
ok    = false(size(Bi));
yinit = @(x) [3*x^2-2*x^3; 6*x-6*x^2; 6-12*x; -12; 0.5*(1-x); -0.5; x; 1; 0];
s     = bvpinit(linspace(0,1,101), yinit);

% Pass 1: increasing Bi, previous Bi's solution is the initial guess for the next one
for i = 1:numel(Bi)
    lastwarn('');
    try
        s1   = bvp4c(odefun, @(ya,yb) bcfun(ya,yb,Bi(i)), s, opts);
        conv = isempty(lastwarn);          % bvp4c warns when it cannot meet RelTol
    catch
        conv = false;
    end
    if conv
        sol{i} = s1;  ok(i) = true;  s = s1;
    else
        fprintf('Bi = %g not reached from the previous solution (retried below)\n', Bi(i));
        s = bvpinit(linspace(0,1,101), yinit);   % restart the next Bi from the simple guess
    end
end

% Pass 2: the low-theta branch has a turning point at Bi ~ 0.099, so a Bi past it
% cannot be reached from below. Continue down from the nearest converged larger Bi.
for i = fliplr(find(~ok))
    j = find(ok & Bi > Bi(i), 1);
    if isempty(j), continue; end
    fprintf('Bi = %g -> continuation from Bi = %g\n', Bi(i), Bi(j));
    s = sol{j};
    for Bstep = logspace(log10(Bi(j)), log10(Bi(i)), 21)
        lastwarn('');
        try
            s    = bvp4c(odefun, @(ya,yb) bcfun(ya,yb,Bstep), s, opts);
            conv = isempty(lastwarn);
        catch
            conv = false;
        end
        if ~conv, break; end
    end
    if conv, sol{i} = s;  ok(i) = true; end
end

for i = 1:numel(Bi)
    if ~ok(i), warning('Bi = %g: not converged', Bi(i)); continue; end
    s = sol{i};
    fprintf('Bi = %5.3f :  f''''(0) = %10.6f   theta(0) = %8.5f   phi''(0) = %8.5f   p(1) = %10.5f   (mesh %d pts)\n', ...
            Bi(i), s.y(3,1), s.y(5,1), s.y(8,1), s.y(9,end), numel(s.x));
end

%% ---------------- Plots (one figure per quantity) ----------------
names = {'$f(\eta)$','$f''(\eta)$','$\theta(\eta)$','$\phi(\eta)$','$p(\eta)$'};
rows  = [1 2 5 7 9];
for k = 1:numel(rows)
    figure(k); hold on; box on; grid on;
    h = [];  labs = {};
    for i = find(ok)
        Y = deval(sol{i}, eta);
        h(end+1) = plot(eta, Y(rows(k),:), 'LineWidth', 1.8, 'Color', cols(i,:)); %#ok<SAGROW>
        labs{end+1} = sprintf('$Bi = %g$', Bi(i));                                %#ok<SAGROW>
    end
    xlabel('$\eta$','Interpreter','latex','FontSize',14);
    ylabel(names{k},'Interpreter','latex','FontSize',14);
    if ~isempty(h)
        legend(h, labs, 'Interpreter','latex','Location','best');
    end
    set(gca,'FontSize',12);
end
