function [A,scalefact] = get_polynomial_function(order,dims,num_dims)
%% Multi-echo MR Thermometry in the upper leg at 7T using near-harmonic 2D reconstruction for initialization
%% HIMM: Harmonic Initialized Model-based Multi-echo
% get_polynomial_function: function that builds the exponent array recursively 
% and use those terms to construct the polynomial A array that contains various 
% polynomial orders in its dimensions
%
% Creator: Mathijs Kikken (University Medical Center Utrecht)
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%
% Inputs:
%   order:    	scalar integer, defines the total (maximum) order 
%   num_dims:   integer, defines the total number of dimensions in 
%               the image (excluding time dimensions)
%   dims:       array of integers, defines the number of pixels in each 
%               physical dimension (1 x N)
%
% Output: 
%   A:          polynomial array
%   scalefact:	scaling factor to scale the independent variables back from unit variance


% Build the exponent array (acquire modelterms)
terms = get_polynomial_terms(order,num_dims);

% Use modelterms to construct the polynomial array
[A,scalefact] = construct_polynomial_matrix(terms,dims,num_dims);

end

function [terms] = get_polynomial_terms(order,num_dims)
% Build the exponent array recursively

% Inputs:
%   order:    	scalar integer, defines the total (maximum) order 
%   num_dims:   integer, defines the total number of dimensions in 
%               the image (excluding time dimensions)
%
% Output: 
%   terms:    	exponent array for the model

% Build the exponent array recursively
if num_dims == 0                   	% terminal case
    terms = [];
elseif (order == 0)             % terminal case
    terms = zeros(1,num_dims);
elseif (num_dims==1)              	% terminal case
    terms = (order:-1:0)';
else                            % general recursive case
    terms = zeros(0,num_dims);
    for k = order:-1:0
        t = get_polynomial_terms(order-k,num_dims-1);
        nt = size(t,1);
        terms = [terms; [repmat(k,nt,1),t]];
    end
end
end

function [A,scalefact] = construct_polynomial_matrix(terms,dims,num_dims)
% Construct polynomial array A with varous polynomial orders in its dimenions
%
% Inputs:
%   terms:    	exponent array for the model
%   dims:       array of integers, defines the number of pixels in each 
%               physical dimension (1 x N)
%   num_dims:   integer, defines the total number of dimensions in 
%               the image (excluding time dimensions)
%
% Output: 
%   A:        	polynomial array
%   scalefact:	scaling factor to scale the independent variables back from unit variance

% Determine a grid with correct dimensions
[x_grid,y_grid] = meshgrid(linspace(0,1,dims(1)),linspace(0,1,dims(2)));
grid = [x_grid(:),y_grid(:)];

% Calculate scaling factor to scale the independent variables to unit variance
scale = sqrt(diag(cov(grid)));
if any(scale==0)
    warning 'Constant terms in the model must be processed differenrly'
    scale(scale==0) = 1;
end
grid_scaled = grid*diag(1./scale);     % Scaled variables

% Create polynomial matrix of given order (defined by modelterms)
A = ones(dims(1)*dims(2),size(terms,1));
scalefact = ones(1,size(terms,1));
for i = 1:size(terms,1)
    for j = 1:num_dims
        A(:,i) = A(:,i).*grid_scaled(:,j).^terms(i,j);
        scalefact(i) = scalefact(i)/(scale(j)^terms(i,j));
    end
end
end
