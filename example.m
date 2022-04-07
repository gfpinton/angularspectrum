% example script for marmottant simulations

%% incident pressure

% define waveform parameters
freq = 2;    % incident frequency in [MHz]
ncyc = 3;    % number of cycles
pnp = 10;   % desired peak negative pressure in [kPa]
pol = 1;    % = 1 or -1, polarity of pulse

% generate waveform
[p0,p0t] = generate_waveform(freq,ncyc,pnp,pol);
dT=p0t(2)-p0t(1); 
f=(0:length(p0t)-1)/(length(p0t)-1)/dT;
plot(f,dbzero(abs(fft(p0))))

%% marmottant simulation

% define input parameters 
distance = .010;     % distance to compute scattered pressure in [m]
bubbleType = 'def';  % type of bubble to determine shell parameters
R0 = 1.75;           % initial bubble radius in [um]
freq = 2;            % incident frequency in [MHz]

tic
[time,rad,pscat,psurf] = marmottant(distance,bubbleType,R0,freq,p0,p0t);
toc

%% plot results

% radius vs time
figure(1);
plot(time*(1e6),rad*1e6/R0,'k-','LineWidth',1);
xlabel('Time (\mus)'); ylabel('R/R_{0}'); box on;
xlim([time(1) time(end)]*(1e6)); title('Radius');
set(gca,'FontSize',12);

% scattered pressure vs time
figure(2);
plot(time*(1e6),pscat,'k-','LineWidth',1);
xlabel('Time (\mus)'); ylabel('Pressure (Pa)'); box on;
xlim([time(1) time(end)]*(1e6)); title('Scattered Pressure');
set(gca,'FontSize',12);

% surface pressure vs time
figure(3);
plot(time*(1e6),psurf,'k-','LineWidth',1);
xlabel('Time (\mus)'); ylabel('Pressure (Pa)'); box on;
xlim([time(1) time(end)]*(1e6)); title('Surface Pressure');
set(gca,'FontSize',12);


%% interpolate for even time steps
pscat2=interp1(time,pscat,p0t);
dT=p0t(2)-p0t(1); 
f=(0:length(p0t)-1)/(length(p0t)-1)/dT;
plot(f,dbzero(abs(fft(p0)))), hold on
plot(f,dbzero(abs(fft(pscat2)))), hold off
grid on, xlim([0 5*freq*1e6])


