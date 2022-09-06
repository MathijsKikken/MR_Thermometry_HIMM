function coefficients = c_fit(mask, depvar, dims, A)
%% Multi-echo MR Thermometry in the upper leg at 7T using near-harmonic 2D reconstruction for initialization
%% HIMM: Harmonic Initialized Model-based Multi-echo
% c_fit: fits a general polynomial regression model in n dimensions, making
% use of QR decomposition.
%
% c_fit fits a polynomial regression model of one or more
% independent variables, of the general form:
%
%   z = f(x,y,...) + error
%
% Creator: Mathijs Kikken (University Medical Center Utrecht)
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%
% Inputs:
%   mask:             	num_dims-dimensional array in which points that are known have a value of 1 
%                       Used to trim N_full (located in A) down to N
%   depvar:             (N x 1 or 1 x N) vector - dependent variable
%                       length(depvar) must be N.
%   dims:               array of integers, defines the number of pixels in each 
%                       physical dimension (1 x N)
%   A:                  (N_full x G) Polynomial array
%                       N_full is the total number of physical data points in the image
%                       G is a number that directly follows from the order of the fit
%
% Output: 
%   coefficients:       Fitted regression coefficients


% Only 1 dependent variable allowed at a time
if size(depvar,2)~=1
    error 'Only 1 dependent variable allowed at a time.'
end

if sum(mask,'all')~=size(depvar,1)
    error 'indepvar (from mask) and depvar are of inconsistent sizes.'
end

% Only take pixels in A of which the off-resonance value is known and to be fitted by coefficients 
A_2d = reshape(A,[dims(1) dims(2) size(A,2)]);      % Reshape to physical dimensions
mask = repmat(mask,[1 1 size(A_2d,3)]);             % Repeat mask in dimensions of all orders
A_known = A_2d(mask == 1);                         	% Only keep all pixels that are within the masked region
A_known = reshape(A_known,[size(A_known,1)/size(A_2d,3),size(A_2d,3)]); % Reshape back to physical dimensions

% Estimate the model using QR: provide a covariance matrix when all done. 
% Use a pivoted QR for maximum stability
[Q,R,E] = qr(A_known,0);

% Use the QR decomposition together with the dependent variables to
% estimate the coefficients
coefficients(E) = R\(Q'*depvar);
coefficients = coefficients';
