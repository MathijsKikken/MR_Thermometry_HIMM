function m = m_update(imgsw,Zw,m,tes,prc,algp,mask,body_mask,signal_mask,fat_mask,R)
%% Multi-echo MR Thermometry in the upper leg at 7T using near-harmonic 2D reconstruction for initialization
%% HIMM: Harmonic Initialized Model-based Multi-echo
% m_update: function to fit PRF shift due to heating
%
% Creators: Megan Poorman, William Grissom (Vanderbilt University Institute of Imaging Science)
% Updated by: Mathijs Kikken (University Medical Center Utrecht)
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%
% Inputs:
%   imgsw:          Dynamic multi-echo images water component (Nx x Ny x Necho)
%   Zw:             Model image for each echo and each baseline (Nx x Ny x Necho x Nbaseline)
%   m:              initial guess
%   tes:            echo times
%   prc:            direction of precession
%   algp.lam:       algorithm params sparsity constraint
%   mask:           masked heat
%   signal_mask:    region without signal (so also no heat)
%   fat_mask:       region with high fat content
%   R:              roughness penalty
%
% Output: 
%   m:        dim x dim fit


% reshape images and model into vectors
innprod = Zw.*conj(imgsw);

% loop over echoes to calc derivatives and curvatures
t1 = abs(innprod);
t2 = angle(innprod);
t3 = sin(t2);
t4 = t1.*t3;
hesssum = 0;gradsum = 0;
for ii = 1:length(tes)
    hesssum = hesssum + t4(:,:,ii)./t2(:,:,ii)*tes(ii)^2;
    gradsum = gradsum - (1-2*prc)*t4(:,:,ii)*tes(ii);
end
if exist('R','var') % add roughness penalty
    if ~isempty(R)
        hesssum = hesssum + R.denom(R,m);
        gradsum = gradsum + R.cgrad(R,m);
    end
end
% if exist('lam','var') && exist('w2','var') % add sparsity penalty
%     hesssum = hesssum + lam*w2; 
%     gradsum = gradsum + lam*w2.*m;
% end
if exist('mask','var') % apply mask
    if ~isempty(mask)
        gradsum = gradsum.*mask;
    end
end
if exist('signal_mask','var') % apply secondary signal mask
    if ~isempty(signal_mask)
        gradsum = gradsum.*abs(1-signal_mask);
    end
end
if exist('signal_mask','var') % apply tertiary signal mask
    if ~isempty(body_mask)
        gradsum = gradsum.*body_mask;
    end
end

% update m
dm = -gradsum./hesssum;
dm(isnan(dm)) = 0;
m = m + dm;

% Enforce limits
m(m < -100) = 0;
m(m >  100) = 0;

% %Enforce negativity (temp increase -> freq decrease)
% m(m > 0) = 0;