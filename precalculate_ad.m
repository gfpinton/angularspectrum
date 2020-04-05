function [afilt3d] = precalculate_ad(alpha0,nX,nY,nT,dZ,dT);

disp(['Gianmarco Pinton, written on 2017-05-25'])
disp(['Precalculating attenuation/dispersion filter...'])
f=(0:nT-1)/nT/dT;
alpha=alpha0/1e12*1e2/(20*log10(exp(1)))*f.^2;
afilt = (exp(-alpha*dZ));% water 
afilt(find(afilt<0))=0;
afilt3d=zeros(nX,nY,nT);
for i=1:nX
    for j=1:nY
        afilt3d(i,j,:)=afilt;
    end
end


disp(['done.'])
