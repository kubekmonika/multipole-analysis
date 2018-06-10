function [a, A] = scattcoefficienta(l, k, e, r, theta, phi)
%Oblicza wspolczynnik rozproszenia a(l).
%
%   SCATTCOEFFICIENTA(l, m, E, R, theta, phi)
% 
%   l - stopien, l >= |m|
%   k - wartosc wektora falowego
%   e - rozklad pola (w postaci wektora)
%   r - promien sfery
%   phi - katy zenitalne
%   theta - katy azymutalne

a = 0;
A = ones(1, 3);
for m = -l : 1 : l
    alm = a_lm(l, m, k, e, r, theta, phi);
    A(m+2) = alm;
    a = a + k^2 * l * (l + 1) * (alm * conj(alm));
end
end

function alm = a_lm(l, m, k, e, r, theta, phi)
% Zwraca wartosc wspolczynnika rozproszenia a_lm
assert(l >= abs(m), 'Blad: |m| > l')
% VSH
% if m >= 0
%     N = vsh('N', l, m, theta, phi, R, k);
% else
%     N = (-1)^abs(m) * conj(vsh('N', l, abs(m), theta, phi, R, k));
% end
N = vshN(m, theta, phi, r, k);

% przeksztalcamy pole
E = sphreshapefield(e, length(theta), length(phi));
% zastepujemy NaN wartoscia 0
N(isnan(N)) = 0;
E(isnan(E)) = 0;
% licznik
licznik = dot(N, E, 3) .* sin(theta);
mianownik = dot(N, N, 3) .* sin(theta);
alm = sum(licznik(:)) / sum(mianownik(:));
end