function [p0,p0t] = generate_waveform(freq,ncyc,pnp,pol)
%%% function to generate analytical waveform for input to Marmottant or AS
%%% outputs: 
%%% p0 = pressure waveform
%%% p0t = time for pressure waveform
%%% inputs: 
%%% freq = frequency in [MHz]
%%% ncyc = number of cycles (three is ~= to experimental "single-cycle")
%%% pnp  = peak negative pressure in [kPa]
%%% pol  = polarity of waveform (1 = pos first, -1 = neg first)

% convert to base units
freq = freq*1e6;  % convert to [Hz]
pnp = pnp*1e3;    % convert to [Pa]

% define time
dt = 1/(20*50e6);          % time step = 10x Nyquist for 50 MHz max freq in spectrum
Tedge_start = (10)/freq;   % time before pulse
Tedge_end = (20)/freq;     % time after pulse
Tref = (ncyc)/freq+Tedge_end;         % reference time
n = (Tref+Tedge_start)/dt;            % number of samples
p0t = linspace(-Tedge_start,Tref,n);  % time values

% find indices for beginning and end of waveform
i_start = find(p0t == 0 | p0t>0, 1, 'first');    % index for start of pulse (time = 0)
i_end = find(p0t == Tref-Tedge_end | p0t > Tref-Tedge_end, 1, 'first'); % index for end of pulse

% define waveform
p0 = zeros(size(p0t)); % initialize pressure vector
if pol == -1
    p0(i_start:i_end) = pnp.*(sin(2*pi*freq*p0t(i_start:i_end)+pi)); % pressure values
elseif pol == 1
    p0(i_start:i_end) = pnp.*(sin(2*pi*freq*p0t(i_start:i_end))); % pressure values
end

% apply Gaussian window to waveform
p0(i_start:i_end) = p0(i_start:i_end).*gausswin(i_end-i_start+1)';

% % plot waveform
% figure; plot(p0t,p0,'k'); xlabel('Time (s)'); ylabel('Pressure (Pa)'); 

