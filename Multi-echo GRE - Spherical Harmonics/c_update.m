function c = c_update(imgs,Z,A,c,tes,prc)
%% Multi-echo MR Thermometry in the upper leg at 7T using near-harmonic 2D reconstruction for initialization
%% HIMM: Harmonic Initialized Model-based Multi-echo
% c_update: function to find the coefficients on the backgrond polynomial fit to
% background phase
%
% Creators: Megan Poorman, William Grissom (Vanderbilt University Institute of Imaging Science)
% Updated by: Mathijs Kikken (University Medical Center Utrecht)
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%
% NEWTON's METHOD GRADIENT DESCENT
%
% Inputs:
%   imgs: Dynamic multi-echo images (Nx x Ny x Necho)
%   Z:    Model image for each echo and each baseline (Nx x Ny x Necho x Nbaseline)
%   A:    polynomial basis functions
%   c:    previous coefficient estimate
%   tes:  echo times
%   prc:  direction of precession
%
% Output: 
%   c:    polynomial coefficents


% reshape images and model into vectors
innprod = Z.*conj(imgs);
innprod = permute(innprod,[3 1 2]); 
innprod = innprod(:,:).';

% loop over echoes to calc derivatives and curvatures
t1 = (abs(innprod));
t2 = (angle(innprod)); %wrapped w/ heat

t3 = sin(t2);
t4 = t1.*t3;

hesssum = 0;gradsum = 0;
for ii = 1:length(tes)    
    hesssum = hesssum + t4(:,ii)./t2(:,ii)*(tes(ii))^2;
    gradsum = gradsum - (1-2*prc)* t4(:,ii)*(tes(ii));%
end

hesssum(isnan(hesssum)) = 0;
% gradsum(isnan(gradsum)) = 0;
dc = -(A'*bsxfun(@times,hesssum,A))\(A'*gradsum);
dc(isnan(dc)) = 0;
c = c + dc;
