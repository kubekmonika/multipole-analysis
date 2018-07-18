function [b, B] = scattCoefficientB(e, k, r, theta, phi, w)
% Scattering coefficient b(1).
%
%   SCATTCOEFFICIENTB(e, k, r, theta, phi, w)
% 
%   e - electric field; as a vector
%   k - wavenumber
%   r - radius
%   theta - elevation, 0 <= theta <= pi
%   phi - azimuth, 0 <= phi <= 2pi
%   w - weights

% l - degree, l >= |m|
l = 1;
b = 0;
B = ones(1, 3);
% removing NaNs
e(isnan(e)) = 0;
for m = -l : 1 : l
    assert(l >= abs(m), 'Blad: |m| > l')
    % harmonic
    [~, M] = vsh(m, theta, phi, r, k);
    % replacing NaN with 0
    M(isnan(M)) = 0;
    % integrand in the nominator
    nominator = dot(M, e, 2) .* w;
    % integrand in the denominator
    denominator = dot(M, M, 3) .* w;
    % counting the integrals
    blm = sum(nominator(:)) / sum(denominator(:));
    % scattering coefficient
    B(m+2) = blm;
    b = b + k^2 * l * (l + 1) * (blm * conj(blm));
end
end