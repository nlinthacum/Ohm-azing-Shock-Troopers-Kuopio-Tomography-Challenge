function result = DBO(f,M,T ,fft_beta)
%==========================================================================
% This function computes the action on the input vector f of the operator
% in the Dbar equation, ( result = [I - PTR rho]f ).
%
% Authors:              Melody Alsaker, Jennifer Mueller. 
% Date Modified:        January 2013
%
% Inputs:   f =         Function on which the operater acts.
%                       Input as a vector of stacked real and complex
%                       parts, total size 2*M^2.
%                   
%           M =         Defines vector and matrix sizes
%
%           T =         tR ./(4*pi*k) .* exp(-1i*(k.*z + conj(k).*conj(z)))      
%              
%          fft_beta =   fft of the Green's function beta
%                   
% Outputs:  result =    [I - PT_r rho]f  is a vector of stacked real and
%                       imaginary parts of total size 2*M^2
%
% External files:       None                         
%==========================================================================

% Re-complexify f and convert to a matrix of size MxM (needed for fft)
f = f(1:M^2) + 1i * f(M^2+1:2*M^2);
matf = reshape(f,M,M);

% Construct conj(matf) .* T with zero padding to accommodate convolution
temp_fT = conj(matf) .* T;
temp_fT_Big  = zeros(2*M);
temp_fT_Big((M/2+1):(3*M/2),(M/2+1):(3*M/2)) = temp_fT;  % zero padding

% Compute action of operator on f. Note, the h^2 is already included in 
% beta, but it could go here instead
temp_fT_Big = ifft2(fft_beta.*fft2(temp_fT_Big)); 
tmp = temp_fT_Big((M/2+1):(3*M/2),(M/2+1):(3*M/2)); %Remove zero padding
f = f - tmp(:);

%Stack real and imaginary parts
result = [real(f); imag(f)];

