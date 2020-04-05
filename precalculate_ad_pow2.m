function [afilt3d,f,attenuation,dispersion] = precalculate_ad_pow2(alpha0,nX,nY,nT,dZ,dT,c0,f0,pow)
% Explanation of Outputs
% afilt3d : the frequency multiplicand e^(-(a+ja*)dx)
% alpha, a : the attenuation in Np/m
% alphaStar, a*: \frac{\omega}{c(\omega)}-\frac{\omega}{c(\omega_0)} because this is a retarded time simulation
% dispersion: c(\omega), the speed of frequency \omega in m/s
% f : frequency vector

disp(['Gianmarco Pinton, written on 2017-05-25'])
disp(['Precalculating attenuation/dispersion filter with power law...'])
f=(0:nT-1)/nT/dT;

alphaUnitConv=alpha0/(1e6^pow)*1e2/(20*log10(exp(1))); % dB/MHz^{pow}/cm * MHz^{pow}/Hz^{pow} * cm/m * Np/dB = Np/(Hz^{pow} m)

alpha=alphaUnitConv*f.^pow;

if mod(pow,2)==1
    alphaStar0=(-2*alphaUnitConv/((2*pi)^pow)/pi).*(log(2*pi*f)-log(2*pi*f0));  
else
    alphaStar0=(alphaUnitConv/(2*pi)^pow)*tan(pi*pow/2).*((2*pi.*f).^(pow-1)-(2*pi.*f0)^(pow-1));
end

alphaStar=2*pi.*alphaStar0.*f;
dispersion=((1/c0)+(alphaStar/2/pi./f)).^-1;
alphaStar(1)=0;

afilt = (exp((-alpha-1j*alphaStar).*dZ));% water
%afilt(find(afilt<0))=0;
afilt3d=ones(nX,nY).*permute(afilt,[1,3,2]);

% for i=1:nX
%     for j=1:nY
%         afilt3d(i,j,:)=afilt;
%     end
% end



disp(['done.'])
end

