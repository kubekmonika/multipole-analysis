function [sca, scb] = coefForEveryWavelength(E, r, eneis, N, theta, phi, w)
% Calculate scattering coefficients.
%   
%   COEFFOREVERYWAVELENGTH(E, r, eneis, N, theta, phi, w)
% 
%   E - electric field given in spherical coordinates
%   r - radius of the sphere
%   eneis - wavelength
%   N - refractive index in the surrounding medium; 
%       real values for each wavelength (an array or a scalar);
%   theta - elevation, 0 <= theta <= pi
%   phi - azimuth, 0 <= phi <= 2pi
%   w - weights for each point

assert(isreal(N), 'N must be real - there is no theory provided for a different case')
% function which calculates the value of the wavevector
wavevector = @(n, enei) 2 * pi * n / ((enei) * 0.1^9);

if isscalar(N)
    N = ones(size(eneis)) * N;
end

n_eneis = length(eneis);
% result arrays
sca = zeros(n_eneis, 1);
scb = zeros(n_eneis, 1);
% value for every wavelength
for i = 1: n_eneis
    k = wavevector(N(i), eneis(i));
    e = E(:,:,i);
    % scattering coefficients
    sca(i) = scattCoefficientA(e, k, r, theta, phi, w);
    scb(i) = scattCoefficientB(e, k, r, theta, phi, w);
end
end
