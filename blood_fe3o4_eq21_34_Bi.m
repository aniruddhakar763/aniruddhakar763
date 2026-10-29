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
%  Written as a function file (helpers below are subfunctions) so it runs in any
%  MATLAB version. Run the whole file (F5 / Run, or type its name in the Command
%  Window), not section by section.
function blood_fe3o4_eq21_34_Bi
clc; close all;

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
TP.rho_f = 1050;    TP.rho_p = 5180;     % kg/m^3
TP.cp_f  = 3900;    TP.cp_p  = 670;      % J/(kg K)
TP.k_f   = 0.492;   TP.k_p   = 6.0;      % W/(m K)
TP.sig_f = 0.8;     TP.sig_p = 2500;     % S/m
TP.bt_f  = 3.8e-4;  TP.bt_p  = 1.3e-5;   % 1/K (thermal expansion)

P = setOmega(P, omega, TP);
fprintf('omega = %4.2f :  L1=%.4f  L2=%.4f  L3=%.4f  L4=%.4f  L5=%.4f\n', ...
        P.omega, P.L1, P.L2, P.L3, P.L4, P.L5);

%% ---------------- Solve for each Bi ----------------
eta   = linspace(0,1,801);
cols  = lines(numel(Bi));
opts  = bvpset('RelTol',1e-6,'AbsTol',1e-8,'NMax',20000);
sol   = cell(size(Bi));
ok    = false(size(Bi));
% y = [f, f', f'', f''', theta, theta', phi, phi', p]
yinit = @(x) [3*x.^2-2*x.^3; 6*x-6*x.^2; 6-12*x; -12; 0.5*(1-x); -0.5; x; 1; 0];
s     = bvpinit(linspace(0,1,101), yinit);

% Pass 1: increasing Bi, previous Bi's solution is the initial guess for the next one
for i = 1:numel(Bi)
    [s1, conv] = solveBi(Bi(i), s, P, opts);
    if conv
        sol{i} = s1;  ok(i) = true;  s = s1;
    else
        s = bvpinit(linspace(0,1,101), yinit);   % restart the next Bi from the simple guess
    end
end

% Pass 2: the low-theta branch has a turning point at Bi ~ 0.099, so a Bi past it
% cannot be reached from below. Continue down from the nearest converged larger Bi.
for i = fliplr(find(~ok))
    j = find(ok & Bi > Bi(i), 1);
    if isempty(j), continue; end
    fprintf('Bi = %g not reached from below -> continuation from Bi = %g\n', Bi(i), Bi(j));
    s = sol{j};
    for Bstep = logspace(log10(Bi(j)), log10(Bi(i)), 21)
        [s, conv] = solveBi(Bstep, s, P, opts);
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
    for i = find(ok)
        Y = deval(sol{i}, eta);
        plot(eta, Y(rows(k),:), 'LineWidth', 1.8, 'Color', cols(i,:), ...
             'DisplayName', sprintf('$Bi = %g$', Bi(i)));
    end
    xlabel('$\eta$','Interpreter','latex','FontSize',14);
    ylabel(names{k},'Interpreter','latex','FontSize',14);
    legend('Interpreter','latex','Location','best','FontSize',12);
    set(gca,'FontSize',12);
end

end

%% ---------------- One bvp4c solve at a given Bi ----------------
function [s, conv] = solveBi(BiVal, s0, P, opts)
    P.Bi = BiVal;
    ws = warning('off','all');            % failed attempts are reported by the caller
    try
        s = bvp4c(@(x,y) odefun(x,y,P), @(ya,yb) bcfun(ya,yb,P), s0, opts);
        conv = s.stats.maxres <= opts.RelTol;
    catch
        s = s0;  conv = false;
    end
    warning(ws);
end

%% ---------------- omega-dependent quantities -> L1..L5, Brinkman factor ----------------
function P = setOmega(P, w, TP)
    n = 2;                                        % shape factor
    rho_nf   = (1-w)*TP.rho_f + w*TP.rho_p;
    rhocp_nf = (1-w)*TP.rho_f*TP.cp_f + w*TP.rho_p*TP.cp_p;
    bt_nf    = (1-w)*TP.bt_f + w*TP.bt_p;
    k_nf   = TP.k_f  *(TP.k_p + (n-1)*TP.k_f - (n-1)*w*(TP.k_p-TP.k_f)) ...
                     /(TP.k_p + (n-1)*TP.k_f + w*(TP.k_p-TP.k_f));
    sig_nf = TP.sig_f*(TP.sig_p + (n-1)*TP.sig_f - (n-1)*w*(TP.sig_p-TP.sig_f)) ...
                     /(TP.sig_p + (n-1)*TP.sig_f + w*(TP.sig_p-TP.sig_f));

    P.omega = w;
    P.muf   = (1-w)^(-2.5);                      % Brinkman factor 1/(1-omega)^2.5
    P.L1 = rho_nf/TP.rho_f;                      % density ratio
    P.L2 = bt_nf/TP.bt_f;                        % thermal-expansion ratio
    P.L3 = rhocp_nf/(rho_nf*TP.cp_f);            % specific-heat ratio
    P.L4 = k_nf/TP.k_f;                          % conductivity ratio
    P.L5 = sig_nf/TP.sig_f;                      % electrical-conductivity ratio
end

%% ---------------- ODE system ----------------
% y = [f, f', f'', f''', theta, theta', phi, phi', p]
function dy = odefun(x,y,P)
    f=y(1); fp=y(2); fpp=y(3); fppp=y(4); th=y(5); thp=y(6); php=y(8);
    St=P.St; b=P.beta; K=P.K;

    % (25) temperature -> theta''
    thpp = ( P.L1*P.L3*P.Pr*St*(x*thp - f*thp) ...
           - P.Pr*(P.Nt*thp^2 + P.Nb*php*thp) - P.L5*P.gam ) / P.L4;

    % (26) concentration -> phi''
    phpp = St*P.Sc*(x*php - f*php) - (P.Nt/P.Nb)*thpp;

    % (23) electro-osmotic potential gradient psi'
    psip = K*( P.R*cosh(K*x) - cosh(K*(1-x)) )/sinh(K);

    % (21) momentum -> f''''
    LHS  = P.L1*St*( 3*fpp + x*fppp + fp*fpp - f*fppp );
    buoy = P.L1*P.L2*P.lam*P.delta*( thp + P.N*php );
    EO   = P.delta*P.Ue*K^2*psip;
    visc = exp(-b*th)*P.muf;              % = exp(-beta*theta)/(1-omega)^2.5
    fpppp = (LHS - buoy - EO)/visc ...
          + 2*b*thp*fppp + b*thpp*fpp - b^2*thp^2*fpp;

    % (24) pressure -> dp/deta
    dp = -P.L1*St*( -f - x*fp + f*fp ) + visc*( 2*b*thp*fp - fpp );

    dy = [fp; fpp; fppp; fpppp; thp; thpp; php; phpp; dp];
end

%% ---------------- Boundary conditions (27)-(34) ----------------
function res = bcfun(ya,yb,P)
    slip0 = P.delta*P.Ue*exp(P.beta*ya(5))/P.muf;   % delta*Ue*e^{beta theta}(1-omega)^2.5
    slip1 = P.delta*P.Ue*exp(P.beta*yb(5))/P.muf;
    res = [ ya(2) - 1 - slip0*ya(3);                % (27) f'(0) = 1 + ... f''(0)
            ya(1);                                  % (28) f(0) = 0
            P.L4*ya(6) + P.Bi*(1 - ya(5));          % (29) L4 theta'(0) = -Bi(1-theta(0))
            ya(7);                                  % (30) phi(0) = 0
            yb(2) - slip1*yb(3);                    % (31) f'(1) = ... f''(1)
            yb(1) - 1;                              % (32) f(1) = 1
            yb(5);                                  % (33) theta(1) = 0
            yb(7) - 1;                              % (34) phi(1) = 1
            ya(9) ];                                % reference pressure p(0) = 0
end
