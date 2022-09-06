function A = get_spherical_function(order,dims)
%% Multi-echo MR Thermometry in the upper leg at 7T using near-harmonic 2D reconstruction for initialization
%% HIMM: Harmonic Initialized Model-based Multi-echo
% get_polynomial_function: function that acquires an A matrix with in each 
% dimension a function of specific order and degrees. The spherical harmonic 
% consist of a combination of spherical functions and Bessel spherical functions.
%
% Creator: Mathijs Kikken (University Medical Center Utrecht)
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%
% Inputs:
%   order:    	order
%   dims:       array of integers, defines the number of pixels in each 
%               physical dimension (1 x N)
%
% Output: 
%   A:          polynomial array


% Create 2D cartesian and spherical grids for the spherical functions
th = linspace(0,pi,dims(2));    % inclination
phi = linspace(0,2*pi,dims(1)); % azimuth
[th,phi] = meshgrid(th,phi);

szn = 1;
x_cart = linspace(-szn,szn,dims(1));
y_cart = linspace(-szn,szn,dims(2));
[x_cart,y_cart] = meshgrid(x_cart,y_cart);
r = sqrt(x_cart.^2+y_cart.^2);

% Initialize the A matrix
A = [];
for n = 0:order
    for m = -n:n
        spherical_function_n_m = build_spherical_function_n_m(n,m,r,th,phi);
        A = [ A spherical_function_n_m(:) ];
    end
end
end


function f = build_spherical_function_n_m(n,m,r,th,phi)
% Acquire a combination of spherical harmonics and Bessel spherical
% functions. These functions of various order and degree can be combined to
% accurately approximate the desired off-resonance map
%
% Inputs:
%   n:    		order
%   m:          degree
%   th:         inclination
%   phi:        azimuth
%
% Output: 
%   f:          sphericalal harmonic function of n order and m degrees
    
    % Cartesian grid: get spherical bessel function of first kind with order m
    J = get_spherical_bessel_firstkind(n,r);
    
    % Spherical grid: get spherical harmomic of order m and degree n
    Y = get_spherical_harmonic(n,m,th,phi,'type','real');
    
    % Construct the spherical function consisting of the 2
    f = J.*Y;

end


function J = get_spherical_bessel_firstkind(m,x)
% Returns the spherical Bessel functions jnu(x)
% https://arxiv.org/pdf/2102.02634.pdf
% https://stackoverflow.com/questions/20515829/matlab-custom-spherical-bessel-function-nan-at-0
%   x is a vector or it may be a matrix if nu is a scalar
%   if nu is a row and x a column vector, the output js is a matrix

    if isscalar(m) && isscalar(x)
        J = 0;
    elseif isscalar(m)
        J = zeros(size(x));
        m = J+m;
    else
        J = zeros(size(m));
        x = J+x;
    end
    x0 = (abs(x) < realmin);
    x0mlt0 = (x0 & m < 0);
    x0nueq0 = (x0 & m == 0);
    J(x0mlt0) = Inf;          % Re(m) < 0, X == 0
    J(x0nueq0) = 1;            % Re(m) == 0, X == 0
    i = ~x0mlt0 & ~x0nueq0 & ~(x0 & m > 0) & (abs(x) < realmax);
    J(i) = sign(x(i)).*sqrt(pi./(2*x(i))).*besselj(m(i)+0.5,x(i));
    
end
    

function Y = get_spherical_harmonic(n,m,th,phi,varargin)
%HARMONICY  Spherical harmonic function.
%
%   Y = HARMONICY(N,M,TH,PHI) computes the surface spherical harmonic of
%   degree N and order M, evaluated at each element of inclination TH and
%   azimuth PHI. N and M must be scalar integers where M <= abs(N). TH and
%   PHI must be arrays of equal size.
%
%   Y = HARMONICY(N,M,TH,PHI,R) computes the solid spherical harmonic of
%   degree N and order M, evaluated at each element of inclination TH,
%   azimuth PHI and radius R. N and M must be scalar integers where M <=
%   abs(N). TH, PHI and R must be arrays of equal size.
%
%   Y = HARMONICY(__,Name,Value) specifies additional options using one or
%   more name-value pair arguments. Valid options are:
%
%     - 'type' specifies whether to compute the complex spherical harmonics
%       or their real part. Valid values are 'complex' (default) and
%       'real'. Real spherical harmonics are of cosine type for M > 0 and
%       of sine type for M < 0.
%
%     - 'norm' specifies whether the result of the computation is to be
%       normalized. The normalization coefficient is chosen so as to ensure
%       that the spherical harmonics are orthonormal. Valid values are true
%       (default), false.
%       
%     - 'phase' specifies whether to include the Condon-Shortley phase
%       term. This term is not strictly necessary but may be useful in
%       quantum mechanics applications. Valid values are true (default),
%       false.
%
%   See also LEGENDRE.
%   Copyright (C) 2018-2020 Javier Montalt Tordera.
% first argument in varargin to be part of name,value pairs
    s = 1;
    if nargin > 4
        % if more than four arguments were passed, check varargin{1}
        if ischar(varargin{1})
            % if char array, this is a name and r is empty
            r = [];
        else
            % if not, this should be r
            r = varargin{1};
            % and name,value pairs, if any, start at the second element
            s = 2;
        end
    else
        % if only four arguments were passed, r is empty
        r = [];
    end
    % if r is not empty, we're computing a solid harmonic
    issolid = ~isempty(r);
    % check inputs
    validateattributes(n,{'numeric'},{'scalar','nonempty','integer','nonnegative'},1);
    validateattributes(m,{'numeric'},{'scalar','nonempty','integer','>=',-n,'<=',n},2);
    validateattributes(th,{'numeric'},{},3);
    validateattributes(phi,{'numeric'},{'size',size(th)},4);
    if issolid
        validateattributes(r,{'numeric'},{'size',size(th)},5);
    end
    % parse optional inputs
    ip = inputParser;
    addParameter(ip,'type','complex',@(x)any(validatestring(x,{'complex','cplx','real'})));
    addParameter(ip,'norm',true,@(x)islogical(x));
    addParameter(ip,'phase',true,@(x)islogical(x));
    parse(ip,varargin{s:end});
    pp = ip.Results;
    % check if m is odd and negative
    isoddm = mod(m,2) == 1;
    isnegm = m < 0;
    % calculate for m >= 0 for now
    m = abs(m);
    % to linear arrays
    sz = size(th);
    th = th(:);
    phi = phi(:);
    r = r(:);
    % normalization factor
    if pp.norm
        P = legendre(n,cos(th),'norm');
        C = 1/sqrt(2*pi);
    else
        P = legendre(n,cos(th),'unnorm');
        C = 1;
    end
    % associated Legendre function
    P = P(abs(m)+1,:)';
    switch pp.type
        case {'real'}
            if isnegm
                E = sin(m*phi);
            else
                E = cos(m*phi);
            end
        case {'complex','cplx'}
            E = exp(1i*m*phi);
            if isnegm
                E = conj(E);
                if isoddm
                    E = -E;
                end
            end
    end

    % surface spherical harmonics
    Y = C * P .* E;
    % solid spherical harmonics
    if issolid
        Y = Y .* r.^n;
    end
    % include Condon-Shortley phase term
    % this term is already included in MATLAB's associated Legendre function
    % depending on whether it is normalized or not
    if isoddm && ((~pp.phase && ~pp.norm) || (pp.phase && pp.norm))
        Y = -Y;
    end
    % reshape to original size
    Y = reshape(Y,sz);
    
end