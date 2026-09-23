%% Blood-Fe3O4 nanofluid: streamlines psi = X*f(eta) for St = -1 and St = 1
%  Solves the non-dimensional system (momentum, energy, concentration) with bvp4c
%  and plots the stream function psi(X,eta) = X*f(eta), X in [0,10], eta in [0,1],
%  in figure 1 (St = -1) and figure 2 (St = 1).
%
%  For St < 0 the energy equation is stiff (Pr = 21) and theta ~ 1.1 in the core,
%  so each case is solved with a homotopy in Pr, 7 -> 10 -> 14 -> 18 -> 21.
%  Solving directly at Pr = 21 does not converge for St = -1.
clear; clc; close all;

%% ---------------- Parameters ----------------
P.delta = 0.1;   % aspect ratio
P.alpha = 0.3;   % variable-viscosity parameter
P.lam   = 5.0;   % mixed convection (Gr/Re)
P.Res   = 0.5;   % stretching Reynolds number
P.N     = 0.5;   % buoyancy ratio
P.R     = 1.0;   % zeta-potential ratio
P.Re    = 1.0;   % electro-osmotic Reynolds number
P.K     = 10;    % electro-osmotic parameter
P.Pr    = 21.0;  % Prandtl number of blood
P.Nt    = 0.3;   % thermophoresis
P.Nb    = 0.3;   % Brownian motion
P.gam   = 0.5;   % Joule heating
P.Sc    = 10.0;  % Schmidt number
P.Bi    = 1.0;   % Biot number

omega = 0.02;             % nanoparticle volume fraction
StVec = [-1 1];           % squeeze number: figure 1 -> St = -1, figure 2 -> St = 1

%% ---------------- Thermophysical properties ----------------
% Base fluid: blood         Nanoparticle: Fe3O4 (magnetite)
TP.rho_f = 1063;    TP.rho_p = 5180;     % kg/m^3
TP.cp_f  = 3594;    TP.cp_p  = 670;      % J/(kg K)
TP.k_f   = 0.492;   TP.k_p   = 9.7;      % W/(m K)
TP.sig_f = 0.8;     TP.sig_p = 0.74e6;   % S/m
TP.bt_f  = 0.18e-5; TP.bt_p  = 1.3e-5;   % 1/K (thermal expansion)

%% ---------------- Solve for each St ----------------
P = setOmega(P, omega, TP);
fprintf('omega = %4.2f :  L1=%.4f  L2=%.4f  L3=%.4f  L4=%.4f  L5=%.4f\n', ...
        P.omega, P.L1, P.L2, P.L3, P.L4, P.L5);
PrPath   = [7 10 14 18 P.Pr];    % homotopy in Pr (needed for St < 0)
PrTarget = P.Pr;
opts     = bvpset('RelTol',1e-6,'AbsTol',1e-8,'NMax',20000);
sol      = cell(size(StVec));
ok       = false(size(StVec));

for i = 1:numel(StVec)
    P.St  = StVec(i);
    yinit = @(x) [P.St*x; P.St; 0; 0; 0.5*(1-x); -0.5; x; 1];
    s     = bvpinit(linspace(0,1,101), yinit);
    success = true;
    for Pr = PrPath
        P.Pr = Pr;
        try
            s = bvp4c(@(x,y) odefun(x,y,P), @(ya,yb) bcfun(ya,yb,P), s, opts);
        catch ME
            success = false;
            warning('St = %g: %s at Pr = %g', P.St, ME.message, Pr);
            break
        end
        if s.stats.maxres > opts.RelTol
            success = false;
            warning('St = %g: residual %.2e at Pr = %g -- not converged', P.St, s.stats.maxres, Pr);
            break
        end
    end
    P.Pr = PrTarget;
    if success
        sol{i} = s;  ok(i) = true;
        fprintf('St = %4.1f :  f''''(0) = %10.6f   theta(0) = %8.5f   phi''(0) = %8.5f   (mesh %d pts)\n', ...
                P.St, s.y(3,1), s.y(5,1), s.y(8,1), numel(s.x));
    end
end

%% ---------------- Streamlines psi = X*f(eta) ----------------
Xv     = linspace(0,10,201);
etav   = linspace(0,1,401);
[XX,EE] = meshgrid(Xv, etav);
nLev   = 12;                     % contour levels on each side of zero

for i = find(ok)
    f   = deval(sol{i}, etav);
    f   = f(1,:).';              % column: f(eta)
    PSI = XX .* repmat(f, 1, numel(Xv));
    tol = 1e-3*max(abs(PSI(:)));

    figure(i); hold on; box on;
    h = gobjects(0);
    if max(PSI(:)) > tol                                        % psi > 0 : solid
        lev = linspace(0, max(PSI(:)), nLev+1);  lev = lev(2:end);
        [~,hp] = contour(XX, EE, PSI, lev, 'LineStyle','-', 'LineWidth',1.3);
        hp.DisplayName = '$\psi > 0$';  h(end+1) = hp;
    end
    if min(PSI(:)) < -tol                                       % psi < 0 : dashed
        lev = linspace(min(PSI(:)), 0, nLev+1);  lev = lev(1:end-1);
        [~,hn] = contour(XX, EE, PSI, lev, 'LineStyle','--', 'LineWidth',1.3);
        hn.DisplayName = '$\psi < 0$';  h(end+1) = hn;
    end
    fin = f(2:end-1);                                           % interior sign change -> dividing streamline
    if any(fin > tol/10) && any(fin < -tol/10)
        [~,h0] = contour(XX, EE, PSI, [0 0], 'LineColor','k', 'LineWidth',2);
        h0.DisplayName = '$\psi = 0$';  h(end+1) = h0;
    end
    colormap(gca, parula);
    cb = colorbar;  cb.Label.String = '$\psi$';  cb.Label.Interpreter = 'latex';  cb.Label.FontSize = 14;
    xlabel('$X$','Interpreter','latex','FontSize',14);
    ylabel('$\eta$','Interpreter','latex','FontSize',14);
    title(sprintf('Streamlines $\\psi = X f(\\eta)$, $St = %g$', StVec(i)), ...
          'Interpreter','latex','FontSize',14);
    legend(h, 'Interpreter','latex','Location','northwest','FontSize',12);
    set(gca,'FontSize',12);  axis([0 10 0 1]);
end

%% ---------------- omega-dependent quantities -> L1..L5, Brinkman factor ----------------
function P = setOmega(P, w, TP)
    rho_nf   = (1-w)*TP.rho_f + w*TP.rho_p;
    rhocp_nf = (1-w)*TP.rho_f*TP.cp_f + w*TP.rho_p*TP.cp_p;
    rhobt_nf = (1-w)*TP.rho_f*TP.bt_f + w*TP.rho_p*TP.bt_p;
    k_nf   = TP.k_f  *(TP.k_p+2*TP.k_f-2*w*(TP.k_f-TP.k_p)) ...
                     /(TP.k_p+2*TP.k_f+w*(TP.k_f-TP.k_p));                         % Maxwell
    sig_nf = TP.sig_f*(TP.sig_p+2*TP.sig_f-2*w*(TP.sig_f-TP.sig_p)) ...
                     /(TP.sig_p+2*TP.sig_f+w*(TP.sig_f-TP.sig_p));                 % Maxwell

    P.omega = w;
    P.muf   = (1-w)^(-2.5);                      % Brinkman factor 1/(1-omega)^2.5
    P.L1 = rho_nf/TP.rho_f;                      % density ratio
    P.L2 = (rhobt_nf/rho_nf)/TP.bt_f;            % thermal-expansion ratio
    P.L3 = (rhocp_nf/rho_nf)/TP.cp_f;            % specific-heat ratio
    P.L4 = k_nf/TP.k_f;                          % conductivity ratio
    P.L5 = sig_nf/TP.sig_f;                      % electrical-conductivity ratio
end

%% ---------------- ODE system ----------------
% y = [f, f', f'', f''', theta, theta', phi, phi']
function dy = odefun(x,y,P)
    f=y(1); fp=y(2); fpp=y(3); fppp=y(4); th=y(5); thp=y(6); php=y(8);
    St=P.St; a=P.alpha; K=P.K;

    % energy -> theta''
    thpp = ( P.L1*P.L3*P.Pr*(St/2*x*thp - thp*f) ...
           - P.Pr*(P.Nt*thp^2 + P.Nb*php*thp) - P.L5*P.gam ) / P.L4;

    % concentration -> phi''
    phpp = St*P.Sc*x/2*php - P.Sc*f*php - (P.Nt/P.Nb)*thpp;

    % electro-osmotic body force
    EO = P.delta*P.Re*K^3*( -cosh(K*(1-x)) + P.R*cosh(K*x) )/sinh(K);

    % momentum -> f''''
    LHS  = P.L1*( St/2*(3*fpp + x*fppp) + fp*fpp - f*fppp );
    buoy = P.L1*P.lam*P.Res*( P.L2*thp + P.N*php );
    visc = exp(-a*th)*P.muf;              % = exp(-a*theta)/(1-omega)^2.5
    fpppp = (LHS - buoy - EO)/visc ...
          + 2*a*thp*fppp + a*thpp*fpp - a^2*thp^2*fpp;

    dy = [fp; fpp; fppp; fpppp; thp; thpp; php; phpp];
end

%% ---------------- Boundary conditions ----------------
function res = bcfun(ya,yb,P)
    slip0 = P.Res + P.delta*P.Re*exp(P.alpha*ya(5))/P.muf;
    slip1 = P.Res + P.delta*P.Re*exp(P.alpha*yb(5))/P.muf;
    res = [ ya(2) - slip0;                          % f'(0)
            ya(1);                                  % f(0) = 0
            P.L4*ya(6) + P.Bi*(1 - ya(5));          % L4 theta'(0) = -Bi(1-theta(0))
            ya(7);                                  % phi(0) = 0
            yb(2) - slip1;                          % f'(1)
            yb(1) - P.St;                           % f(1) = St
            yb(5);                                  % theta(1) = 0
            yb(7) - 1 ];                            % phi(1) = 1
end
