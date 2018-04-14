function a = scattcoefficienta(n, k, e, R, theta, phi)
%Oblicza wspolczynnik rozproszenia a(n).
%
%   SCATTCOEFFICIENTA(n, m, E, R, theta, phi)
% 
%   n - stopien, n >= |m|
%   k - wartosc wektora falowego
%   e - rozklad pola (w postaci wektora)
%   R - promien sfery
%   theta - katy azymutalne
%   phi - katy zenitalne
% 
a = 0;
for m = -n : 1 : n
    anm = a_nm(m, n, k, e, R, theta, phi);
    a = a + k^2 * n * (n + 1) * (anm * conj(anm));
end

end

function anm = a_nm(m, n, k, e, R, theta, phi)
% Zwraca wartosc wspolczynnika rozproszenia a_nm
assert(n >= abs(m), 'Blad: |m| > n')
% VSH
N = vsh('N', n, m, theta, phi, R, k);
% przeksztalcamy pole
E = sphreshapefield(e, length(theta), length(phi));
% zastepujemy NaN wartoscia 0
N(isnan(N)) = 0;
E(isnan(E)) = 0;
% licznik 
licznik = dot(E, conj(N), 3) .* sin(theta);
mianownik = dot(N, conj(N), 3) .* sin(theta);
anm = sum(licznik(:)) / sum(mianownik(:));
end