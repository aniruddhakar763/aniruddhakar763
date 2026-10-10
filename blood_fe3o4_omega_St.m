%% Blood-Fe3O4 nanofluid flow in a microchannel: effect of omega and St
%  Solves the non-dimensional momentum, temperature and concentration
%  equations with bvp4c for omega = [0 0.01 0.02 0.04] and St = [-1 1].
%  Figures 1-4 : f, f', theta, phi
%  Figure  5   : dp/deta, built from f, f', f'', theta, theta'
%  Prints L1..L5 for every omega and the wall values (eta = 0 and eta = 1)
%  of f, f', f'', dp/deta, theta, theta', phi, phi' for every case.
clear; clc; close all;

%% ---------------- Parameters ----------------
P.beta  = 0.1;   % variable-viscosity parameter
P.lam   = 5.0;   % mixed convection (Gr/Re)
P.N     = 0.1;   % buoyancy ratio
P.delta = 0.05;  % aspect ratio (multiplies buoyancy, EO force and slip)
P.Pr    = 21.0;  % Prandtl number of blood
P.Nt    = 0.2;   % thermophoresis
P.Nb    = 0.3;   % Brownian motion
P.gam   = 0.5;   % Joule heating
P.Sc    = 10.0;  % Schmidt number
P.Bi    = 0.1;   % Biot number
P.K     = 10;    % electro-osmotic parameter
P.Ue    = 1.0;   % electric field parameter
P.R     = 1.5;   % zeta-potential ratio zeta_u/zeta_l

omegaVec = [0 0.01 0.02 0.04];   % nanoparticle volume fraction
StVec    = [-1 1];               % squeezing parameter

%% ---------------- Thermophysical properties ----------------
% Base fluid: blood         Nanoparticle: Fe3O4
TP.rho_bf = 1050;    TP.rho_n = 5180;     % kg/m^3
TP.cp_bf  = 3900;    TP.cp_n  = 670;      % J/(kg K)
TP.k_bf   = 0.492;   TP.k_n   = 6.0;      % W/(m K)
TP.sig_bf = 0.8;     TP.sig_n = 2500;     % S/m
TP.bt_bf  = 3.8e-4;  TP.bt_n  = 1.3e-5;   % 1/K
TP.n      = 2;                            % shape factor (Hamilton-Crosser)

%% ---------------- L1..L5 ----------------
fprintf('Thermophysical ratios\n');
fprintf('%8s %10s %10s %10s %10s %10s\n','omega','L1','L2','L3','L4','L5');
for i = 1:numel(omegaVec)
    Q = setOmega(P, omegaVec(i), TP);
    fprintf('%8.2f %10.6f %10.6f %10.6f %10.6f %10.6f\n', ...
            Q.omega, Q.L1, Q.L2, Q.L3, Q.L4, Q.L5);
end
fprintf('\n');

%% ---------------- Solve every (St, omega) case ----------------
eta   = linspace(0,1,801);
opts  = bvpset('RelTol',1e-6,'AbsTol',1e-8,'NMax',50000);
yinit = @(x) [3*x.^2-2*x.^3; 6*x-6*x.^2; 6-12*x; -12; 0.5*(1-x); -0.5; x; 1];
nO = numel(omegaVec);  nS = numel(StVec);
sol = cell(nS,nO);  ok = false(nS,nO);  Q = cell(1,nO);

for j = 1:nS
    s = bvpinit(linspace(0,1,201), yinit);
    for i = 1:nO
        Q{i} = setOmega(P, omegaVec(i), TP);
        Q{i}.St = StVec(j);
        try
            % previous omega's solution is the initial guess for the next one
            s = bvp4c(@(x,y) odefun(x,y,Q{i}), @(ya,yb) bcfun(ya,yb,Q{i}), s, opts);
        catch ME
            warning('St = %g, omega = %g: %s', StVec(j), omegaVec(i), ME.message);
            s = bvpinit(linspace(0,1,201), yinit);
            continue
        end
        if s.stats.maxres > opts.RelTol
            warning('St = %g, omega = %g: residual %.2e -- not converged', ...
                    StVec(j), omegaVec(i), s.stats.maxres);
        end
        sol{j,i} = s;  ok(j,i) = true;
    end
end

%% ---------------- Wall values ----------------
hdr = {'f','f''','f''''','dp/deta','theta','theta''','phi','phi'''};
for w = 1:2
    fprintf('Values at eta = %d\n', w-1);
    fprintf('%5s %6s', 'St', 'omega');
    fprintf(' %12s', hdr{:});  fprintf('\n');
    for j = 1:nS
        for i = 1:nO
            if ~ok(j,i), continue, end
            y  = deval(sol{j,i}, w-1);
            dp = dpdeta(w-1, y, Q{i}, StVec(j));
            fprintf('%5g %6.2f', StVec(j), omegaVec(i));
            fprintf(' %12.6f', [y(1) y(2) y(3) dp y(5) y(6) y(7) y(8)]);
            fprintf('\n');
        end
    end
    fprintf('\n');
end

%% ---------------- Plots ----------------
cols   = [0.00 0.45 0.74;  0.85 0.33 0.10;  0.47 0.67 0.19;  0.49 0.18 0.56;
          0.93 0.69 0.13;  0.30 0.75 0.93;  0.64 0.08 0.18;  0.00 0.00 0.00];
styles = {'-','--',':','-.','-','--',':','-.'};
names  = {'$f(\eta)$','$f''(\eta)$','$\theta(\eta)$','$\phi(\eta)$', ...
          '$\partial p/\partial\eta$'};
rows   = [1 2 5 7 0];          % 0 -> pressure gradient

for k = 1:5
    figure(k); hold on; box on; grid on;
    c = 0;
    for j = 1:nS
        for i = 1:nO
            c = c + 1;
            if ~ok(j,i), continue, end
            Y = deval(sol{j,i}, eta);
            if rows(k) > 0
                v = Y(rows(k),:);
            else
                v = dpdeta(eta, Y, Q{i}, StVec(j));
            end
            plot(eta, v, 'LineStyle', styles{c}, 'Color', cols(c,:), ...
                 'LineWidth', 1.8, 'DisplayName', ...
                 sprintf('$St = %g,\\ \\omega = %g$', StVec(j), omegaVec(i)));
        end
    end
    xlabel('$\eta$','Interpreter','latex','FontSize',14);
    ylabel(names{k},'Interpreter','latex','FontSize',14);
    legend('Interpreter','latex','Location','best','FontSize',11);
    set(gca,'FontSize',12);
end

%% ---------------- omega-dependent quantities -> L1..L5, Brinkman factor ----------------
function P = setOmega(P, w, TP)
    n        = TP.n;
    rho_nbf  = (1-w)*TP.rho_bf + w*TP.rho_n;
    rhocp    = (1-w)*TP.rho_bf*TP.cp_bf + w*TP.rho_n*TP.cp_n;
    bt_nbf   = (1-w)*TP.bt_bf + w*TP.bt_n;
    k_nbf    = TP.k_bf*(TP.k_n+(n-1)*TP.k_bf-(n-1)*w*(TP.k_n-TP.k_bf)) ...
                      /(TP.k_n+(n-1)*TP.k_bf+w*(TP.k_n-TP.k_bf));
    sig_nbf  = TP.sig_bf*(TP.sig_n+(n-1)*TP.sig_bf-(n-1)*w*(TP.sig_n-TP.sig_bf)) ...
                        /(TP.sig_n+(n-1)*TP.sig_bf+w*(TP.sig_n-TP.sig_bf));

    P.omega = w;
    P.B  = (1-w)^2.5;                     % Brinkman: mu_nbf = mu_bf/(1-omega)^2.5
    P.L1 = rho_nbf/TP.rho_bf;             % density ratio
    P.L2 = bt_nbf/TP.bt_bf;               % thermal-expansion ratio
    P.L3 = rhocp/(rho_nbf*TP.cp_bf);      % specific-heat ratio
    P.L4 = k_nbf/TP.k_bf;                 % thermal-conductivity ratio
    P.L5 = sig_nbf/TP.sig_bf;             % electrical-conductivity ratio
end

%% ---------------- ODE system ----------------
% y = [f, f', f'', f''', theta, theta', phi, phi']
function dy = odefun(x,y,P)
    f=y(1); fp=y(2); fpp=y(3); fppp=y(4); th=y(5); thp=y(6); php=y(8);
    St=P.St; b=P.beta; K=P.K;

    % temperature -> theta''
    thpp = ( P.L1*P.L3*P.Pr*St*(x*thp - f*thp) ...
           - P.Pr*(P.Nt*thp^2 + P.Nb*php*thp) - P.L5*P.gam ) / P.L4;

    % concentration -> phi''
    phpp = St*P.Sc*(x*php - f*php) - (P.Nt/P.Nb)*thpp;

    % electro-osmotic potential derivative psi'(eta)
    psip = K*(P.R*cosh(K*x) - cosh(K*(1-x)))/sinh(K);

    % momentum -> f''''
    LHS  = P.L1*St*(3*fpp + x*fppp + fp*fpp - f*fppp);
    buoy = P.L1*P.L2*P.lam*P.delta*(thp + P.N*php);
    EO   = P.delta*P.Ue*K^2*psip;
    visc = exp(-b*th)/P.B;                % e^(-beta theta)/(1-omega)^2.5
    fpppp = (LHS - buoy - EO)/visc ...
          + 2*b*thp*fppp + b*thpp*fpp - b^2*thp^2*fpp;

    dy = [fp; fpp; fppp; fpppp; thp; thpp; php; phpp];
end

%% ---------------- Pressure gradient dp/deta ----------------
% works for a single point (y column) or a whole profile (Y matrix)
function dp = dpdeta(x,Y,P,St)
    f=Y(1,:); fp=Y(2,:); fpp=Y(3,:); th=Y(5,:); thp=Y(6,:);
    visc = exp(-P.beta*th)/P.B;
    dp = -( P.L1*St*(-f - x.*fp + f.*fp) - visc.*(2*P.beta*thp.*fp - fpp) );
end

%% ---------------- Boundary conditions ----------------
function res = bcfun(ya,yb,P)
    s0 = P.delta*P.Ue*exp(P.beta*ya(5))*P.B;
    s1 = P.delta*P.Ue*exp(P.beta*yb(5))*P.B;
    res = [ ya(2) - 1 - s0*ya(3);                   % f'(0) = 1 + dUe e^(b th) (1-w)^2.5 f''(0)
            ya(1);                                  % f(0) = 0
            P.L4*ya(6) + P.Bi*(1 - ya(5));          % L4 theta'(0) = -Bi(1-theta(0))
            ya(7);                                  % phi(0) = 0
            yb(2) - s1*yb(3);                       % f'(1) = dUe e^(b th) (1-w)^2.5 f''(1)
            yb(1) - 1;                              % f(1) = 1
            yb(5);                                  % theta(1) = 0
            yb(7) - 1 ];                            % phi(1) = 1
end
