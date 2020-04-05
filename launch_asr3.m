%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITTEN BY GIANMARCO PINTON 
% FIRST CREATED: 2017-05-01
% LAST MODIFIED: 2020-04-05
% NONLINEAR PROPAGATION IN HOMOGENEOUS MEDIUM
% MODIFIED ANGULAR SPECTRUM 
% RUSANOV
% FREQUENCY DOMAIN ATTENUATION AND DISPERSION
% KRAMERS KRONIG ATTENUATION AND DISPERSION
% ABSORBING BOUNDARY LAYERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load hydrophone measurements of imaging transducer taken near transducer face
load source_pressure_4C1_30V_M1_9_100MHz_single
p0=max(max(max(abs(pressure))));

dTT=1/100e6;
dXX=(elPos(2)-elPos(1))*1e-3;

a=32e-3; % source radius
f0=2.35e6; % center frequency of pulse
omega0 = f0*2*pi;
c0=1500; % speed of sound
lambda=c0/f0; % wavelength
dX=lambda/5; dY=dX;% grid spacing

dZ=dX*8;
dT=dX/5/c0; %2*pi/omega0/20;
k=omega0/c0;
beta=3.5; % nonlinear coefficient
rho0=1000; % equilibrium density
N=beta/2/c0^3/rho0;
p0=1e6; % transmitter pressure
prop_dist=9e-2;

%% INTERPOLATE PRESSURE %%
tmp=interp1easy(pressure(:,round(end/2),round(end/2)),dTT/dT);
p2=zeros(length(tmp),size(pressure,2),size(pressure,3));
for j=1:size(pressure,2)
    for k=1:size(pressure,3)
        p2(:,j,k)=interp1easy(pressure(:,j,k),dTT/dT);
    end
end
p2=single(p2);

outimg = squeeze(p2(round(end/2),:,:));
outimg = interp2easy(outimg,dXX/dX,dXX/dY);
p3 = zeros(size(p2,1),size(outimg,1),size(outimg,2),'single');

for i=1:size(p2,1)
    i
  p3(i,:,:) = interp2easy(squeeze(p2(i,:,:)),dXX/dX,dXX/dY);
end
clear p2;
%p3=p3(100:end,:,:);
%p3=flipdim(p3,1);
p4=permute(p3,[2 3 1]);

    
% GRID SIZE
nX=size(p4,1); nY=size(p4,2); nT=size(p4,3)% grid size
if(mod(nX,2)-1) % keep the grid odd
    nX=nX+1;
    p4(nX-1:nX,:,:)=0;
end
if(mod(nY,2)-1)
    nY=nY+1;
    p4(:,nY-1:nY,:)=0;
end
if(mod(nT,2)-1)
    nT=nT+1;
    p4(:,:,nT-1:nT,:,:)=0;
end


xaxis=(0:nX-1)*dX;xaxis=xaxis-mean(xaxis);
yaxis=(0:nY-1)*dY;yaxis=yaxis-mean(yaxis);
taxis = (0:nT-1)*dT;

apa=p4;
figure(1), imagesc(t,yaxis,(squeeze(apa(:,round(nY/2),:)))), 
title('initial condition (Pa)'), xlabel('t (s)'), ylabel('y (m)'), colorbar

%% PRECALCULATE MODIFIED ANGULAR SPECTRUM 
[HH] = precalculate_mas(nX,nY,nT,dX,dY,dZ,dT,c0);
%% ABSORBING BOUNDARY LAYER
[abl] = precalculate_abl(nX,nY,nT);
%% PRECALCULATE ABSORPTION/DISPERSION FILTER
alpha0=2.17e-3; % dB/MHz^2/cm for water, note attenuation law proportional to f^2 here
[afilt3d] = precalculate_ad(alpha0,nX,nY,nT,dZ,dT);
alpha0=0.3; pow=1; % dB/MHz^pow/cm power law attenuation here
[afilt3d] = precalculate_ad_pow2(alpha0,nX,nY,nT,dZ,dT,c0,omega0/2/pi,pow); 

figure(1), plot(real(squeeze(afilt3d(round(end/2),round(end/2),:))))
figure(2), plot(imag(squeeze(afilt3d(round(end/2),round(end/2),:))))
figure(3), imagesc(real(abl(:,:,round(end/2))))


%% MARCH ALONG PROPAGATION AXIS %%
apaz=apa; zvec=dZ; dZa=dZ; cc=1;
while(sum(zvec)<prop_dist)
    %N*dZ/dT*max(max(max(apaz)))
    if(N*dZ/dT*max(max(max(apaz)))>0.1) %stability criterion
      disp('Stability criterion violated, retrying with smaller step size')
      dZ=0.075*dT/max(max(max(apaz)))/N
      cc=1;
      apaz=apa;
      zvec=dZ;
      
      [HH] = precalculate_mas(nX,nY,nT,dX,dY,dZ,dT,c0);
      [abl] = precalculate_abl(nX,nY,nT);
      [afilt3d] = precalculate_ad_pow2(alpha0,nX,nY,nT,dZ,dT,c0,omega0/2/pi,pow); 
    end
    zvec(cc)=dZ; sum(zvec);
    disp(['Propagation distance = ' num2str(sum(zvec)) ' m'])

    apaz=march_asr2(apaz,dZ,dT,N,HH,abl,afilt3d);
    figure(1), imagesc(t,yaxis,(squeeze(apaz(round(nX/2),:,:)))), xlabel('t (s)'), ylabel('m'), title([num2str(sum(zvec)) 'm']), colorbar
    %str = ['print -djpeg movie/yimg_'  sprintf('%0.4d',cc) '.jpg']; eval(str)
    figure(2), imagesc(t,xaxis,(squeeze(apaz(:,round(nY/2),:)))), xlabel('t (s)'), ylabel('m'), title([num2str(sum(zvec)) 'm']), colorbar
    %str = ['print -djpeg movie/ximg_'  sprintf('%0.4d',cc) '.jpg']; eval(str)
    figure(3), plot(taxis,squeeze(apaz(round(nX/2),round(nY/2),:))), xlabel('t (s)'), ylabel('Pa'), title( [num2str(sum(zvec)) 'm'])
    %str = ['print -djpeg movie/timg_'  sprintf('%0.4d',cc) '.jpg']; eval(str)
    drawnow
    cc=cc+1;
end
