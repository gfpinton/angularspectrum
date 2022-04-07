%% Marmottant Rayleigh-Plesset model

function dy = RP_v2(t,y,pdrivt,pdriv,par)

pdriv = interp1(pdrivt,pdriv,t); % Interpolate the data set (ft,f) at time t

% create variable that remembers previous radius value, e.g. y(1)
persistent yr
if isempty(yr)
    yr(1) = y(1); 
else
    yr = [yr y(1)];
end

% check for shell rupture
yrind = find(yr > par.Rbreakup,1);
if isempty(yrind)
    rupture = 0;
else   % existing index means bubble has passed breakup radius
    rupture = 1; 
end

% define dynamic surface tension at shell
if y(1) < par.Rbuck % Buckled regime
    sigma_R = 0;
elseif rupture == 0 && y(1) >= par.Rbuck && y(1) <= par.Rbreakup % Elastic regime before rupture
    sigma_R = par.chi*((y(1)^2/par.Rbuck^2)-1);
elseif rupture == 1 && y(1) >= par.Rbuck && y(1) <= par.Rrupt % Elastic regime after rupture
    sigma_R = par.chi*((y(1)^2/par.Rbuck^2)-1);
elseif rupture == 1 && y(1) > par.Rrupt % Ruptured regime
    sigma_R = par.Sw;
else
    sigma_R = 0;
end

dy = zeros(2,1); % y(1) = R and y(2) = Rdot
dy(1) = y(2);
dy(2) = ((((par.p0+(2*par.sigma_R0/par.R0))*((par.R0/y(1))^(3*par.gamma))*(1-(3*par.gamma/par.c)*y(2))-(2*sigma_R/y(1))-(4*par.mu*(y(2)/y(1)))-(4*par.kappas*(y(2)/(y(1)^2)))-(par.p0+pdriv))/par.rho)-(3/2)*y(2)^2)/y(1);
end
