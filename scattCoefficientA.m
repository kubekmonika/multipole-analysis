function [a, A] = scattCoefficientA(e, k, r, theta, phi, w)
% Oblicza wspolczynnik rozproszenia a(1).
%
%   SCATTCOEFFICIENTA(e, k, r, theta, phi, w)
% 
%   e - rozklad pola (w postaci wektora)
%   k - wartosc wektora falowego
%   r - promien sfery
%   phi - katy azymutalne
%   theta - katy zenitalne
%   w - wagi dla każdego punktu

% l - stopien, l >= |m|
l = 1;
a = 0;
A = ones(1, 3);
e(isnan(e)) = 0;
for m = -l : 1 : l
    assert(l >= abs(m), 'Blad: |m| > l')
    % liczymy harmonike
    [N, ~] = vsh(m, theta, phi, r, k);
    % zastepujemy NaN wartoscia 0
    N(isnan(N)) = 0;
    % wartosc w liczniku
    licznik = dot(N, e, 2) .* w;
    % wartosc w mianowniku
    mianownik = dot(N, N, 2) .* w;
    % wartosc calej calki
    alm = sum(licznik(:)) / sum(mianownik(:));
    % liczymy wspolczynnik
    A(m+2) = alm;
    a = a + k^2 * l * (l + 1) * (alm * conj(alm));
end
end