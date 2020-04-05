function [abl] = precalculate_abl(nX,nY,nT);

disp(['Precalculating absorbing boundary layer...'])
abl_tmp=ablvec(nX,round(nX/5))*ablvec(nY,round(nY/5))';
abl_vec=ablvec(nT,round(nT/5));
for k=1:nT
    abl(:,:,k)=abl_tmp*abl_vec(k);
end
%imagesc(squeeze(abl(:,round(end/2),:)))
disp(['done.'])
