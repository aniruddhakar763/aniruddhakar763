function blood_fe3o4_St
%% Blood-Fe3O4 nanofluid (modified model): profiles vs squeeze number St
%  Solves the non-dimensional system of modified_model.tex (momentum, energy,
%  concentration) with bvp4c for St = [-1 -0.5 0 0.5 1] and plots f, f',
%  theta, phi in figures 1-4.
%
%  Each St is solved with a homotopy in Pr, 1 -> 5 -> 10 -> 15 -> 21, starting
%  from a simple guess; the energy equation is stiff at Pr = 21.
%
%  Written as a function file (main function + local functions) so it also
%  runs on MATLAB releases older than R2016b, which do not allow functions
%  inside scripts. Run it by typing  blood_fe3o4_St  at the prompt.
clc; close all;

%% ---------------- Parameters ----------------
P.beta  = 0.1;   % variable-viscosity parameter
P.lam   = 5.0;   % mixed convection (Gr/Re)
P.N     = 0.1;   % buoyancy ratio
P.delta = 0.05;  % aspect ratio
P.Pr    = 21.0;  % Prandtl number of blood
P.Nt    = 0.2;   % thermophoresis
P.Nb    = 0.2;   % Brownian motion
P.gam   = 0.5;   % Joule heating
P.Sc    = 10.0;  % Schmidt number
P.Bi    = 0.1;   % Biot number
P.K     = 10;    % electro-osmotic parameter m (= K)
P.Ue    = 1.0;   % electric field parameter
P.R     = 1.0;   % zeta-potential ratio zeta_u/zeta_l

omega = 0.02;                   % nanoparticle volume fraction (fixed)
StVec = [-1 -0.5 0 0.5 1];      % squeeze number

%% ---------------- Thermophysical properties (modified_model.tex) ----------------
% Base fluid: blood         Nanoparticle: Fe3O4 (magnetite)
TP.rho_f = 1050;    TP.rho_p = 5180;     % kg/m^3
TP.cp_f  = 3900;    TP.cp_p  = 670;      % J/(kg K)
TP.k_f   = 0.492;   TP.k_p   = 6.0;      % W/(m K)
TP.sig_f = 0.8;     TP.sig_p = 2500;     % S/m
TP.bt_f  = 3.8e-4;  TP.bt_p  = 1.3e-5;   % 1/K (thermal expansion)

%% ---------------- Solve for each St ----------------
P = setOmega(P, omega, TP);
fprintf('omega = %4.2f :  L1=%.4f  L2=%.4f  L3=%.4f  L4=%.4f  L5=%.4f\n', ...
        P.omega, P.L1, P.L2, P.L3, P.L4, P.L5);
PrPath   = [1 5 10 15 P.Pr];      % homotopy in Pr
PrTarget = P.Pr;
opts     = bvpset('RelTol',1e-6,'AbsTol',1e-8,'NMax',50000);
sol      = cell(size(StVec));
ok       = false(size(StVec));

for i = 1:numel(StVec)
    P.St  = StVec(i);
    yinit = @(x) [x-x^2+P.St*(3*x^2-2*x^3); 1-2*x+P.St*(6*x-6*x^2); -2+P.St*(6-12*x); -12*P.St; ...
                  0.5*(1-x); -0.5; x; 1];
    s     = bvpinit(linspace(0,1,201), yinit);
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

%% ---------------- Plots ----------------
eta   = linspace(0,1,401);
cols  = [0.165 0.471 0.839; 0.922 0.408 0.204; 0.106 0.686 0.478; ...
         0.929 0.631 0.000; 0.290 0.227 0.655];
lsty  = {'-','--','-.',':','-'};
mrk   = {'o','s','^','d','v'};
names = {'$f(\eta)$','$f''(\eta)$','$\theta(\eta)$','$\phi(\eta)$'};
files = {'St_f','St_fp','St_theta','St_phi'};
rows  = [1 2 5 7];
mIdx  = @(i) 1+mod(8*(i-1),40):40:numel(eta);   % marker positions, staggered per curve
for k = 1:4
    figure(k); set(gcf,'Color','w','Units','inches','Position',[1 1 6.4 4.8], ...
                      'PaperUnits','inches','PaperPosition',[0 0 6.4 4.8]);
    hold on; box on; grid on;
    hLeg = [];  labels = {};
    for i = find(ok)
        Y  = deval(sol{i}, eta);
        yk = Y(rows(k),:);
        idx = mIdx(i);
        plot(eta, yk, lsty{i}, 'LineWidth', 2.0, 'Color', cols(i,:));
        plot(eta(idx), yk(idx), mrk{i}, 'LineStyle', 'none', 'Color', cols(i,:), ...
             'MarkerSize', 7, 'MarkerFaceColor', 'w', 'LineWidth', 1.4);
        % invisible line + marker used only for the legend entry
        hLeg(end+1) = plot(NaN, NaN, lsty{i}, 'LineWidth', 2.0, 'Color', cols(i,:), ...
                           'Marker', mrk{i}, 'MarkerSize', 7, 'MarkerFaceColor', 'w'); %#ok<AGROW>
        labels{end+1} = sprintf('$S_t = %g$', StVec(i));                               %#ok<AGROW>
    end
    xlabel('$\eta$','Interpreter','latex','FontSize',16);
    ylabel(names{k},'Interpreter','latex','FontSize',16);
    lg = legend(hLeg, labels, 'Location','northoutside','Orientation','horizontal');
    set(lg,'Interpreter','latex','FontSize',11,'Box','off');
    set(gca,'FontSize',13,'TickDir','in','XMinorTick','on','YMinorTick','on','LineWidth',1.0);
    xlim([0 1]);
    print(gcf, files{k}, '-dpng', '-r400');
end
end   % main function blood_fe3o4_St

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
    P.B  = (1-w)^2.5;                            % (1-omega)^2.5
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
    St=P.St; b=P.beta; K=P.K;

    % energy -> theta''
    thpp = ( P.L1*P.L3*P.Pr*St*(x*thp - thp*f) ...
           - P.Pr*(P.Nt*thp^2 + P.Nb*php*thp) - P.L5*P.gam ) / P.L4;

    % concentration -> phi''
    phpp = St*P.Sc*(x*php - f*php) - (P.Nt/P.Nb)*thpp;

    % electro-osmotic body force delta*Ue*K^2*psi'(eta)
    psip = K*( P.R*cosh(K*x) - cosh(K*(1-x)) )/sinh(K);
    EO   = P.delta*P.Ue*K^2*psip;

    % momentum -> f''''
    LHS  = P.L1*St*( (3*fpp + x*fppp) + fp*fpp - f*fppp );
    buoy = P.L1*P.lam*St*( P.L2*thp + P.N*php );
    visc = exp(-b*th)/P.B;                % = exp(-beta*theta)/(1-omega)^2.5
    fpppp = (LHS - buoy - EO)/visc ...
          + 2*b*thp*fppp + b*thpp*fpp - b^2*thp^2*fpp;

    dy = [fp; fpp; fppp; fpppp; thp; thpp; php; phpp];
end

%% ---------------- Boundary conditions ----------------
function res = bcfun(ya,yb,P)
    dU = P.delta*P.Ue;
    res = [ ya(2) - 1 - dU*exp(P.beta*ya(5))*P.B*ya(3);   % f'(0) = 1 + slip
            ya(1);                                        % f(0) = 0
            P.L4*ya(6) + P.Bi*(1 - ya(5));                % L4 theta'(0) = -Bi(1-theta(0))
            ya(7);                                        % phi(0) = 0
            yb(2) - dU*exp(P.beta*yb(5))*P.B*yb(3);       % f'(1) = slip
            yb(1) - P.St;                                 % f(1) = St
            yb(5);                                        % theta(1) = 0
            yb(7) - 1 ];                                  % phi(1) = 1
end
