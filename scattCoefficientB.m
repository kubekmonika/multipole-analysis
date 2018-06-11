function [b, B] = scattCoefficientB(e, k, r, theta, phi)
%Oblicza wspolczynnik rozproszenia b(1).
%
%   SCATTCOEFFICIENTB(l, m, E, R, theta, phi)
%
%   e - rozklad pola (w postaci wektora)
%   k - wartosc wektora falowego
%   r - promien sfery
%   phi - katy azymutalne
%   theta - katy zenitalne

% l - stopien, l >= |m|
l = 1;
b = 0;
B = ones(1, 3);
e(isnan(e)) = 0;
for m = -l : 1 : l
    assert(l >= abs(m), 'Blad: |m| > l')
    % liczymy harmonike
    M = vsh(m, theta, phi, r, k, 'M');
    % zastepujemy NaN wartoscia 0
    M(isnan(M)) = 0;
    % wartosc w liczniku
    licznik = dot(M, e, 2) .* sin(theta);
    % wartosc w mianowniku
    mianownik = dot(M, M, 3) .* sin(theta);
    % wartosc calej calki
    blm = sum(licznik(:)) / sum(mianownik(:));
    % liczymy wspolczynnik
    B(m+2) = blm;
    b = b + k^2 * l * (l + 1) * (blm * conj(blm));
end
end