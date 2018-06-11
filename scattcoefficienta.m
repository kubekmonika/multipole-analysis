function [a, A] = scattcoefficienta(e, k, r, theta, phi)
% Oblicza wspolczynnik rozproszenia a(1).
%
%   SCATTCOEFFICIENTA(m, E, R, theta, phi)
% 
%   e - rozklad pola (w postaci wektora)
%   k - wartosc wektora falowego
%   r - promien sfery
%   phi - katy azymutalne
%   theta - katy zenitalne

% l - stopien, l >= |m|
l = 1;
a = 0;
A = ones(1, 3);
e(isnan(e)) = 0;
for m = -l : 1 : l
    assert(l >= abs(m), 'Blad: |m| > l')
    % liczymy harmonike
    N = vsh(m, theta, phi, r, k, 'N');
    % zastepujemy NaN wartoscia 0
    N(isnan(N)) = 0;
    % wartosc w liczniku
    licznik = dot(N, e, 2) .* sin(theta);
    % wartosc w mianowniku
    mianownik = dot(N, N, 2) .* sin(theta);
    % wartosc calej calki
    alm = sum(licznik(:)) / sum(mianownik(:));
    % liczymy wspolczynnik
    A(m+2) = alm;
    a = a + k^2 * l * (l + 1) * (alm * conj(alm));
end
end