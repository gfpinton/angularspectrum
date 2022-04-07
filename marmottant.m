function [time,rad,pscat,psurf] = marmottant(distance,bubbleType,R0,freq,pm,taxis)
%%%  Computes  the radial response and scattered pressure from an encapsulated microbubble using 
%%%  the Marmottant model.
%%%  Code modified from Boulder Summer School on Ultrasound Contrast Agents, May 2019.
%%%
%%%  Input variables: 
%%%  distance   = distance from bubble to calculate scattered pressure in [m]
%%%  bubbleType = 'def','mic','div' --- type of bubble to simulate for shell & core parameters
%%%  R0         = initial bubble radius in [um]
%%%  freq       = excitation frequency in [MHz]
%%%  pm         = incident pressure waveform in [Pa]
%%%  taxis      = time axis for pressure waveform (pm) in [s]
%%%
%%%  Output variables:
%%%  time       = time scale for radius data in [s]
%%%  rad        = radius data in [m]
%%%  pscat      = scattered pressure at specified distance in [Pa]
%%%  psurf      = scattered pressure at bubble surface (10 um) in [Pa]


%% Initial parameters

% ------- Medium Parameters ------------------------ %
par.p0 = 99562        ; % [Pa]        ambient pressure in Chapel Hill (486' altitude)
par.Sw = 0.073        ; % [N/m]       surface tension of water
par.c = 1540          ; % [m/s]       speed of sound in water at 20deg C
par.mu = 0.001        ; % [Pas]       medium dynamic viscosity at 20 deg C
par.rho = 1000        ; % [kg/m^3]    water density

% ------- Ultrasound transmission parameters ------- %
par.f = freq*1e6      ; % [Hz]        driving frequency

% ------- Bubble parameters ------------------------ %
par.R0 = R0*1e-6      ; % [m]         bubble radius at equilibrium

switch bubbleType
    
    case 'div'  % DSPC & PEG40S shell, C4F10 core
        
        % core parameters
        par.Dth = 2e-6                             ; % [m^2/s]  thermal diffusivity of gas core (C4F10 = 2e-6)
        par.peclet = (par.R0^2)*2*pi*par.f/par.Dth ; %          Peclet number: Pe >> 1 --> adiabatic, else isothermal
        if par.peclet > 2
            par.gamma = 1.07                       ; %          adiabatic polytropic exponent (air = 1.4, C4F10 = 1.07, C3F8 = 1.06, SF6 = 1.095)
        else
            par.gamma = 1                          ; %          isothermal polytropic exponent
        end
        
        % shell parameters
        if par.R0 > 2
            par.chi = 0.77         ; % [N/m]       shell elasticity
            par.kappas = 5.8*10^-9 ; % [kg/s]      shell viscosity
        elseif par.R0 < 2
            par.chi = 0.88         ; % [N/m]       shell elasticity
            par.kappas = 1.6*10^-9 ; % [kg/s]      shell viscosity
        end
        par.sigma_R0 = 0.04        ; % [N/m]       initial surface tension of shell

    case 'def'  % Definity: DPPC, DPPA, & mPEG5000 shell, C3F8 core
        par.Dth = 2e-6                             ; % [m^2/s]  thermal diffusivity of gas core (C4F10 = 2e-6)
        par.peclet = (par.R0^2)*2*pi*par.f/par.Dth ; %          Peclet number: Pe >> 1 --> adiabatic, else isothermal
        if par.peclet > 2
            par.gamma = 1.06                       ; %          adiabatic polytropic exponent (air = 1.4, C4F10 = 1.07, C3F8 = 1.06, SF6 = 1.095)
        else
            par.gamma = 1                          ; %          isothermal polytropic exponent
        end
        
        % shell parameters
        par.chi = 0.81        ; % [N/m]       shell elasticity
        par.kappas = 6*10^-9  ; % [kg/s]      shell viscosity
        par.sigma_R0 = 0.04   ; % [N/m]       initial surface tension of shell
        
    case 'mic'  % Micromarker: phospholipid shell, C4F10/N2 core
        par.Dth = 2e-6                             ; % [m^2/s]  thermal diffusivity of gas core (C4F10 = 2e-6)
        par.peclet = (par.R0^2)*2*pi*par.f/par.Dth ; %          Peclet number: Pe >> 1 --> adiabatic, else isothermal
        if par.peclet > 2
            par.gamma = 1.07                       ; %          adiabatic polytropic exponent (air = 1.4, C4F10 = 1.07, C3F8 = 1.06, SF6 = 1.095)
        else
            par.gamma = 1                          ; %          isothermal polytropic exponent
        end
        
        % shell parameters
        par.chi = 4.1         ; % [N/m]       shell elasticity
        par.kappas = 5*10^-9  ; % [kg/s]      shell viscosity
        par.sigma_R0 = 0.04   ; % [N/m]       initial surface tension of shell
        
end

%% Initialization

% Driving pulse
pdrivt = taxis; 
pdriv = pm; 

% Computation time
T_start = taxis(1); % start time
T_end = taxis(end); % end time

% Initial conditions for R and Rdot
Initial_con(1,1) = par.R0; % Bubble size
Initial_con(2,1) = 0;      % Bubble wall velocity

%% Solve ODE with Marmottant

% Compute buckling and rupture radii
par.Rbuck = par.R0./sqrt((par.sigma_R0/par.chi)+1); % buckling radius
par.Rrupt = par.Rbuck*sqrt(1+par.Sw/par.chi);  % upper limit radius after rupture
par.Sbreakup = par.chi*(((1.3*par.R0/par.Rbuck)^2)-1); % breakup tension for breakup radius = 1.3R0
par.Rbreakup = par.Rbuck*sqrt(1+par.Sbreakup/par.chi);

% Numerically solve the Rayleigh-Plesset with Marmottant using ode45-solver
% Output: T_m=time, R_m=radius
options  = odeset('RelTol', 1e-6, 'AbsTol', [1e-9 1e-9]', 'Refine', 1);
[T_m, R_m] = ode45(@(t,y) RP(t, y, pdrivt, pdriv, par), [T_start T_end], Initial_con, options);

%% Calculate scattered pressure

% time step for derivatives
time = T_m; 
for i = 2:length(T_m), tstep(i-1) = T_m(i) - T_m(i-1); end

% radius from Marmottant solution
rad = R_m(:,1);                   % [m]

% wall velocity from Marmottant solution
vel = R_m(:,2);                   % [m/s]

% approximate acceleration (Rdotdot) from velocity
acc = diff(vel)./tstep';             % [m/s^2]
acc = padarray(acc,1,0,'pre');  % match size of rad and vel

% compute scattered pressure at surface (10 um from bubble)
psurf = (par.rho/10e-6)*(((rad.^2).*acc)+(rad.*vel.^2));   % [Pa]

% compute scattered pressure at distance
pscat = (par.rho/distance)*(((rad.^2).*acc)+(2.*rad.*vel.^2));   % [Pa]

%% Plot results

if 0
    % plot excitation and scattered pressures
    figure; hold on;
    h = plot(time*(1e6),pdriv/(1e3),'r-','LineWidth',1);
    h2 = plot(time*(1e6),pscat/(1e3),'b-','LineWidth',1);
    xlabel('Time (\mus)'); ylabel('Pressure (kPa)'); box on;
    xlim([time(1) time(end)]*(1e6))
    legend([h, h2],'Excitation','Scattered')
    set(gca,'FontSize',14);
    print('-dpng','-r300',[saveName(1:end-4),'_pressures.png']);
end
  
if 0
    % plot excitation and scattered pressures
    figure; hold on;
    yyaxis left; 
    h = plot(time,pdriv/(1e3),'k-','LineWidth',1);
    ax = gca; ax.YColor = 'k';
    ylabel('Pressure (kPa)');
    set(ax,'FontSize',14);
    ylim([-600 600]);
    yyaxis right;
    h2 = plot(time,pscat/(1e3),'r-','LineWidth',1);
    ax = gca; ax.YColor = 'r';
    ylim([-50 50]);
    xlabel('Time (s)'); ylabel('Pressure (kPa)'); box on;
    xlim([time(1) time(end)])
    legend([h, h2],'Excitation','Scattered')
    set(ax,'FontSize',14);
    print('-dpng','-r300',[saveName(1:end-4),'_pressures.png']);
end

if 0
    % to visualize R, Rdot, and Rdotdot
    figure; plot(time,rad,'k-','LineWidth',1); xlabel('Time (s)'); ylabel('Radius (m)');  
    set(gca, 'FontSize', 14); axis tight
%     print('-dpng','-r300',[saveName(1:end-4),'_rad.png']);
    figure; plot(time,vel,'k-','LineWidth',1); xlabel('Time (s)'); ylabel('Wall Velocity (m/s)'); 
    set(gca, 'FontSize', 14); axis tight
%     print('-dpng','-r300',[saveName(1:end-4),'_vel.png']);
    figure; plot(time,acc,'k-','LineWidth',1); xlabel('Time (s)'); ylabel('Acceleration (m/s^{2})');
    set(gca, 'FontSize', 14); axis tight
%     print('-dpng','-r300',[saveName(1:end-4),'_acc.png']);
end    

if 0
    % plot radius with buckling and rupture radii
    figure;
    plot(time,rad,'k-','LineWidth',1); xlabel('Time (s)'); ylabel('Radius (m)');
    axis tight
    hold on; plot(time,par.Rbuck*ones(1,length(time)),'r-','LineWidth',0.5)
    hold on; plot(time,par.Rrupt*ones(1,length(time)),'b-','LineWidth',0.5)
    legend('R','R_{buckling}','R_{ruptured}')
    set(gca, 'FontSize', 14);
%     print('-dpng','-r300',[saveName(1:end-4),'_rad_buckling.png']);
end    
   

