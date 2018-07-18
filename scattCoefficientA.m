function [a, A] = scattCoefficientA(e, k, r, theta, phi, w)
% Scattering coefficient a(1).
%
%   SCATTCOEFFICIENTA(e, k, r, theta, phi, w)
% 
%   e - electric field; as a vector
%   k - wavenumber
%   r - radius
%   theta - elevation, 0 <= theta <= pi
%   phi - azimuth, 0 <= phi <= 2pi
%   w - weights

% l - degree, l >= |m|
l = 1;
a = 0;
A = ones(1, 3);
% removing NaNs
e(isnan(e)) = 0;
for m = -l : 1 : l
    assert(l >= abs(m), 'Blad: |m| > l')
    % harmonic
    [N, ~] = vsh(m, theta, phi, r, k);
    % replacing NaN with 0
    N(isnan(N)) = 0;
    % integrand in the nominator
    nominator = dot(N, e, 2) .* w;
    % integrand in the denominator
    denominator = dot(N, N, 2) .* w;
    % counting the integrals
    alm = sum(nominator(:)) / sum(denominator(:));
    % scattering coefficient
    A(m+2) = alm;
    a = a + k^2 * l * (l + 1) * (alm * conj(alm));
end
end